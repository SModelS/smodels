#!/usr/bin/env python3

"""
.. module:: hs3Interface
   :synopsis: Code that delegates the computation of limits and likelihoods to
              pyhs3, using the HS3 statistical model format.

.. moduleauthor:: SModelS collaboration

This module is the HS3 counterpart of :mod:`smodels.statistics.pyhfInterface`.
It reads an HS3 workspace (a JSON dictionary), patches the SModelS signal
yields into the HistFactory channels, and uses :mod:`pyhs3` to evaluate the
likelihood.

The public interface deliberately mirrors :class:`PyhfData` and
:class:`PyhfUpperLimitComputer` where this is practical.  In particular,
``Hs3Data`` accepts the same arguments and region mappings as ``PyhfData``.
Region mappings may either use the ``"hs3"`` key, or fall back to the older
``"pyhf"`` key.
"""

import copy
import math
import os
import re
import sys
import time
import warnings
from typing import Dict, List, Optional, Text, Union

import numpy as np
from scipy import optimize
from scipy.stats import norm

from smodels.base.physicsUnits import UnitLumi, UnitXSec
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import (
    CLsfromNLL,
    NllEvalType,
    aposteriori,
    apriori,
    clsType,
    findRoot,
    observed,
)
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError

warnings.filterwarnings("ignore", r"invalid value encountered in log")
warnings.filterwarnings("ignore", r"divide by zero encountered in log")

try:
    import pyhs3
except ModuleNotFoundError:
    print("[SModelS:hs3Interface] pyhs3 import failed. Is the module installed?")
    sys.exit(-1)

try:
    import pytensor.tensor as pt
    from pyhs3.compile import function as pyhs3_function
    from pyhs3.context import Context
    from pyhs3.distributions import HistFactoryDistChannel
    from pytensor.graph.traversal import explicit_graph_inputs
except Exception as e:  # pragma: no cover - only triggered for broken installs
    print(f"[SModelS:hs3Interface] pyhs3/pytensor import failed: {e}")
    sys.exit(-1)


hs3info = {
    "ver": getattr(pyhs3, "__version__", "?"),
    "hasgreeted": False,
}

# Sets the maximum number of attempts for determining the bracketing interval
# for mu in getUpperLimitOnMu.
nattempts_max = 10

poiName = "mu_SIG"
signalUncertaintyName = "signalUncertainty"


def _splitBinName(name: str) -> tuple:
    """Split an SModelS/pyhf style channel name into channel and bin index.

    ``"SR_meff[2]"`` becomes ``("SR_meff", 2)``; a name without a bin suffix
    becomes ``(name, None)``.
    """
    match = re.match(r"^(.*)\[(\d+)\]$", name)
    if match:
        return match.group(1), int(match.group(2))
    return name, None


def _nbins(dist: dict) -> int:
    """Return the total number of bins of a HistFactory distribution."""
    nbins = 1
    for axis in dist.get("axes", []):
        nbins *= int(axis.get("nbins", 1))
    return nbins


def _regionChannelName(region: dict) -> Optional[str]:
    """Return the HS3 channel name for a region mapping.

    New database entries use ``"hs3"``.  For backwards compatibility the
    ``"pyhf"`` key is used as a fallback, since the channel names in an HS3
    workspace typically coincide with the pyhf ones.
    """
    for key in ("hs3", "pyhf", "smodels"):
        value = region.get(key)
        if value is not None:
            return value
    return None


def _constraintSpec(reference: str, distributions: Dict[str, dict]) -> dict:
    """Translate an HS3 constraint reference into a pyhs3-friendly spec."""
    if reference in ("Gauss", "Poisson", "LogNormal"):
        return {"type": reference, "sigma": 1.0, "aux": 0.0}

    dist = distributions.get(reference, {})
    distType = dist.get("type")
    sigma = dist.get("sigma", 1.0)
    spec = {"type": "Gauss", "sigma": 1.0, "aux": 0.0}
    if isinstance(sigma, str):
        spec["sigmaName"] = sigma
    else:
        spec["sigma"] = float(sigma)

    if distType == "gaussian_dist":
        spec["type"] = "Gauss"
        mean = dist.get("mean", 0.0)
        if isinstance(mean, str):
            spec["auxName"] = mean
        else:
            spec["aux"] = float(mean)
    elif distType == "poisson_dist":
        spec["type"] = "Poisson"
        spec["aux"] = 1.0
    elif distType in ("lognormal_dist", "log_normal_dist"):
        spec["type"] = "LogNormal"
        spec["aux"] = 1.0
    elif "poisson" in reference.lower():
        spec["type"] = "Poisson"
        spec["aux"] = 1.0
    elif "lognormal" in reference.lower() or "log_normal" in reference.lower():
        spec["type"] = "LogNormal"
        spec["aux"] = 1.0
    return spec


def _normalizeConstraints(js: dict) -> Dict[str, dict]:
    """Make constraint names palatable for pyhs3.

    The HS3 specification allows constraints to refer to a distribution by
    name.  pyhs3 currently only understands the special literals ``"Gauss"``,
    ``"Poisson"`` and ``"LogNormal"``.  ROOT exports, for instance, use names
    like ``"alpha_JET_...Constraint"`` which refer to a ``gaussian_dist`` in
    the workspace.

    This function rewrites such references in place and returns the constraint
    information (type, width, auxiliary data) needed to build the constraint
    terms with the correct normalisation.  pyhs3 hardcodes unit widths and
    zero auxiliary data for ``Gauss`` constraints, so the interface builds
    those terms itself.

    Constrained ``normfactor`` modifiers are not supported by pyhs3.  Their
    constraint is removed from the modifier, but the information is kept in
    the returned mapping so that it can still be added to the likelihood.
    """
    distributions = {
        d.get("name"): d
        for d in js.get("distributions", [])
        if d.get("name") is not None
    }
    constraintSpecs = {}
    droppedNormFactorConstraint = False
    for dist in js.get("distributions", []):
        if dist.get("type") != "histfactory_dist":
            continue
        for sample in dist.get("samples", []):
            for modifier in sample.get("modifiers", []):
                reference = modifier.get("constraint")
                if reference is None:
                    continue
                spec = _constraintSpec(reference, distributions)
                parameter = modifier.get("parameter")
                if parameter is not None:
                    constraintSpecs.setdefault(parameter, spec)
                if modifier.get("type") == "normfactor":
                    modifier.pop("constraint", None)
                    spec["isNormFactor"] = True
                    droppedNormFactorConstraint = True
                    continue
                modifier["constraint"] = spec["type"]
    if droppedNormFactorConstraint:
        logger.debug(
            "HS3 workspace contains constrained normfactor modifiers; "
            "pyhs3 does not support these, so the constraints are added "
            "to the likelihood by the HS3 interface instead."
        )
    return constraintSpecs


def _modifierParameters(dist: dict) -> List[str]:
    """Collect all parameter names used by the modifiers of a channel."""
    parameters = []
    for sample in dist.get("samples", []):
        for modifier in sample.get("modifiers", []):
            parameter = modifier.get("parameter")
            if parameter is not None:
                parameters.append(parameter)
            parameters.extend(
                name
                for name in modifier.get("parameters", [])
                if name is not None
            )
    return parameters


def _defaultParameterValue(name: str) -> float:
    """A reasonable nominal value for a parameter that is missing."""
    if name == poiName:
        return 1.0
    if name.startswith("gamma") or name == "Lumi":
        return 1.0
    return 0.0


def _defaultParameterBounds(name: str) -> tuple:
    """Reasonable domain bounds for a parameter that is missing."""
    if name == poiName:
        return (-5.0, 10.0)
    if name == signalUncertaintyName:
        return (-5.0, 5.0)
    if name.startswith("gamma") or name == "Lumi":
        return (0.0, 10.0)
    return (-5.0, 5.0)


class Hs3Data:
    """Holds data for use with pyhs3.

    :ivar nsignals: signal predictions dictionary of dictionaries, one for
        each json file, one entry per signal region
    :ivar inputJson: HS3 workspace as a dictionary
    :ivar regions: the SModelS regions
    :ivar globalInfo: None, or link to globalInfo (for debugging)
    :ivar jsonFileName: name of json file (optional)
    """

    def __init__(
        self,
        nsignals: Dict,
        inputJson: dict,
        regions: Optional[list] = None,
        includeCRs: bool = False,
        signalUncertainty: Optional[float] = None,
        globalInfo: Optional[object] = None,
        jsonFileName: Optional[str] = None,
        analysisName: Optional[str] = None,
    ):
        self.globalInfo = globalInfo
        self.nsignals = nsignals
        self.getTotalSignalYield()
        self.jsonFileName = jsonFileName
        self.errorFlag = False
        self.regions = regions
        self.inputJson = inputJson
        self.includeCRs = includeCRs
        self.signalUncertainty = signalUncertainty
        self.analysisName = analysisName
        self.combinations = None
        if jsonFileName is not None:
            self.combinations = os.path.splitext(os.path.basename(jsonFileName))[0]

        # information filled by getWSInfo()
        self.distributions = {}
        self.channelMap = {}
        self.channelsInfo = None

        self.checkConsistency()
        self.getWSInfo()

    def getTotalSignalYield(self) -> float:
        """The total yield in all signal regions."""
        total = 0.0
        for signal in self.nsignals.values():
            try:
                total += signal
            except TypeError:
                raise SModelSError(f"type {self.nsignals} is {type(signal)}")
        self.totalSignalYield = total
        return total

    def getTotalXSec(self) -> UnitXSec:
        """The total fiducial cross section for this computer."""
        return self.totalSignalYield / self.globalInfo.lumi

    def checkConsistency(self) -> None:
        """Check various inconsistencies of the Hs3Data attributes."""
        if not isinstance(self.nsignals, dict):
            logger.error("The 'nsignals' parameter must be a dictionary")
            self.errorFlag = True
            self.zeroSignalsFlag = True
            return

        self.zeroSignalsFlag = all([s == 0 for s in self.nsignals.values()])
        if self.regions is None:
            return

        nBinsJson = 0
        for region in self.regions:
            if region.get("type") == "SR" or (
                region.get("type") == "CR"
                and self.includeCRs
                and region.get("smodels") is not None
            ):
                nBinsJson += 1
        if nBinsJson != len(self.nsignals):
            logger.error(
                f"The number of signals ({len(self.nsignals)}) provided is "
                f"different from the number of signal bins for json "
                f"({nBinsJson}) for {self.regions}"
            )
            self.errorFlag = True

    def _inferRegions(self) -> None:
        """Infer HS3 channel names if the region mappings do not provide them."""
        if not self.distributions:
            self.errorFlag = True
            return

        keys = list(self.nsignals.keys())
        totalBins = sum(_nbins(dist) for dist in self.distributions.values())
        if totalBins != len(keys):
            totalBins = -1

        if self.regions is None:
            self.regions = []
            cursor = 0
            for name, dist in self.distributions.items():
                nbins = _nbins(dist)
                if totalBins > 0:
                    for i in range(nbins):
                        smodelsName = keys[cursor]
                        cursor += 1
                        hs3Name = name if nbins == 1 else f"{name}[{i}]"
                        self.regions.append(
                            {
                                "smodels": smodelsName,
                                "type": "SR",
                                "hs3": hs3Name,
                            }
                        )
                elif len(keys) == len(self.distributions):
                    self.regions.append(
                        {"smodels": keys[cursor], "type": "SR", "hs3": name}
                    )
                    cursor += 1
                else:
                    logger.error(
                        "Cannot infer HS3 region mappings: provide explicit "
                        "region mappings with an 'hs3' or 'pyhf' key."
                    )
                    self.errorFlag = True
                    return
            return

        # Regions are given, but none of them has a channel name.  Assign
        # channels in the order in which they appear in the workspace.
        if totalBins > 0 and len(self.regions) == totalBins:
            cursor = 0
            for name, dist in self.distributions.items():
                nbins = _nbins(dist)
                for i in range(nbins):
                    region = self.regions[cursor]
                    region["hs3"] = name if nbins == 1 else f"{name}[{i}]"
                    cursor += 1
            return
        if len(self.regions) == len(self.distributions):
            for region, name in zip(self.regions, self.distributions):
                region["hs3"] = name
            return

        logger.error(
            "Cannot infer HS3 region mappings: provide explicit region "
            "mappings with an 'hs3' or 'pyhf' key."
        )
        self.errorFlag = True

    def getWSInfo(self) -> None:
        """Extract the channel information from the HS3 workspace."""
        if self.errorFlag:
            return
        if not isinstance(self.inputJson, dict):
            logger.error(
                f"The 'inputJson' parameter must be of type dict not "
                f"{type(self.inputJson)}"
            )
            self.errorFlag = True
            return

        distributions = [
            dist
            for dist in self.inputJson.get("distributions", [])
            if dist.get("type") == "histfactory_dist"
        ]
        self.distributions = {dist["name"]: dist for dist in distributions}
        if not self.distributions:
            logger.error("No 'histfactory_dist' distributions found in HS3 workspace")
            self.errorFlag = True
            return

        if self.regions is None or not any(
            _regionChannelName(region) for region in self.regions
        ):
            self._inferRegions()
        if self.errorFlag:
            return

        # Collect all region mappings per channel.
        for region in self.regions:
            channelName = _regionChannelName(region)
            if channelName is None:
                logger.warning(f"Region {region} has no HS3/pyhf channel name")
                continue
            baseName, binIndex = _splitBinName(channelName)
            if baseName not in self.distributions:
                logger.warning(
                    f"HS3 channel {baseName} from region {region} not found in "
                    f"workspace"
                )
                continue
            entry = self.channelMap.setdefault(
                baseName,
                {
                    "name": baseName,
                    "dist": self.distributions[baseName],
                    "regions": [],
                    "included": False,
                },
            )
            entry["regions"].append((binIndex, region))
            included = region.get("type") == "SR" or (
                region.get("type") == "CR"
                and self.includeCRs
                and region.get("smodels") is not None
            )
            if included:
                entry["included"] = True

        # Turn the region mappings into one signal vector per channel.
        signalRegions = []
        otherRegions = []
        for name, entry in self.channelMap.items():
            dist = entry["dist"]
            nbins = _nbins(dist)
            signals = [0.0] * nbins
            assigned = [False] * nbins
            smodelsNames = []
            regions = entry["regions"]
            directRegions = [(i, r) for i, r in regions if i is None]
            binnedRegions = [(i, r) for i, r in regions if i is not None]
            binnedRegions.sort(key=lambda item: item[0])

            for binIndex, region in binnedRegions:
                smodelsName = region.get("smodels")
                if smodelsName is None:
                    continue
                if binIndex >= nbins:
                    logger.error(
                        f"Region {region} points to bin {binIndex} of channel "
                        f"{name}, which only has {nbins} bins"
                    )
                    self.errorFlag = True
                    continue
                signals[binIndex] = float(self.nsignals.get(smodelsName, 0.0))
                assigned[binIndex] = True
                smodelsNames.append(smodelsName)

            for binIndex, region in directRegions:
                smodelsName = region.get("smodels")
                if smodelsName is None:
                    continue
                if nbins == 1:
                    signals[0] = float(self.nsignals.get(smodelsName, 0.0))
                    assigned[0] = True
                    smodelsNames.append(smodelsName)
                elif len(directRegions) == nbins:
                    # One region per bin, in the order given by the database.
                    pos = directRegions.index((binIndex, region))
                    signals[pos] = float(self.nsignals.get(smodelsName, 0.0))
                    assigned[pos] = True
                    smodelsNames.append(smodelsName)
                else:
                    logger.error(
                        f"Json region {name} has {nbins} bins, but no bin "
                        f"indices are given in the region mappings"
                    )
                    self.errorFlag = True

            entry["nbins"] = nbins
            entry["signals"] = signals
            entry["smodelsNames"] = smodelsNames
            if entry["included"]:
                if not all(assigned):
                    logger.error(
                        f"Json region {name} has {nbins} bins, but only "
                        f"{sum(assigned)} are implemented!"
                    )
                    self.errorFlag = True
                signalRegions.append(
                    {
                        "name": name,
                        "smodelsName": ";".join(smodelsNames),
                        "signal": signals,
                        "size": nbins,
                    }
                )
            else:
                otherRegions.append({"name": name})

        self.channelsInfo = {
            "signalRegions": signalRegions,
            "otherRegions": otherRegions,
        }
        logger.debug(f"HS3 channels: {list(self.channelMap.keys())}")


class _Hs3ModelBundle:
    """A compiled pyhs3 model together with fitting helpers.

    pyhs3 does not provide minimisation or hypothesis-testing helpers, so this
    class compiles the log-probability (and its gradient) with pytensor and
    exposes the small set of numerical operations required by SModelS.
    """

    def __init__(
        self,
        model,
        workspace=None,
        workspaceJson: Optional[dict] = None,
        dataNames: Optional[dict] = None,
        poiName: str = poiName,
        constraintSpecs: Optional[Dict[str, dict]] = None,
        auxOverrides: Optional[Dict[str, float]] = None,
    ):
        self.model = model
        self.workspace = workspace
        self.workspaceJson = workspaceJson or {}
        self.dataNames = dataNames or {}
        self.poiName = poiName
        self.constraintSpecs = constraintSpecs or {}
        self.auxOverrides = auxOverrides or {}

        logProb = self._buildLogProb(model)
        inputs = []
        seen = set()
        for var in explicit_graph_inputs([logProb]):
            if var.name is not None and var.name not in seen:
                seen.add(var.name)
                inputs.append(var)
        self.inputs = inputs
        self.inputNames = [var.name for var in inputs]
        if poiName not in self.inputNames:
            raise SModelSError(
                f"POI '{poiName}' is not a free parameter of the HS3 model"
            )
        self.poiIndex = self.inputNames.index(poiName)

        self._nllFunction = pyhs3_function(
            inputs=inputs,
            outputs=-2.0 * logProb,
            on_unused_input="ignore",
            mode="FAST_RUN",
        )
        self._gradFunction = None
        if len(inputs) <= 12:
            # For small models the compiled gradient is fast and reduces the
            # number of likelihood evaluations.  For larger workspaces the
            # pytensor gradient graph becomes prohibitively expensive to
            # evaluate, so SciPy's numerical differentiation is used instead.
            try:
                # Gradient of the objective actually minimised, i.e. -2 log L.
                grad = pt.grad(-2.0 * pt.sum(logProb), wrt=inputs)
                self._gradFunction = pyhs3_function(
                    inputs=inputs,
                    outputs=grad,
                    on_unused_input="ignore",
                    mode="FAST_RUN",
                )
            except Exception as e:  # pragma: no cover - unusual graphs
                logger.debug(f"Could not compile HS3 gradient: {e}")

        freeParams = model.free_params
        nominalParams = model.nominal_params
        self.defaults = np.array(
            [
                float(freeParams.get(name, nominalParams.get(name, 0.0)))
                for name in self.inputNames
            ],
            dtype=float,
        )
        self.bounds = []
        for name in self.inputNames:
            lower, upper = model.domain.get(name, (None, None))
            defaultLower, defaultUpper = _defaultParameterBounds(name)
            if lower is None:
                lower = defaultLower
            if upper is None:
                upper = defaultUpper
            # A lower bound of exactly zero makes the expected rates vanish at
            # the boundary, so ``log(rate)`` is ``-inf`` and L-BFGS-B's line
            # search gives up immediately.  pyhf uses ``1e-10`` for its
            # constrained parameters; do the same here.
            if float(lower) == 0.0:
                lower = 1e-10
            self.bounds.append((float(lower), float(upper)))

        # Compile the expected-rate expressions.  These are the Poisson means,
        # before the Poisson probability is formed.
        self._expectedFunctions = {}
        context = Context(model.parameters)
        for distObject in model._distribution_objects:
            if not isinstance(distObject, HistFactoryDistChannel):
                continue
            name = distObject.name
            try:
                expression = distObject._compute_expected_rates(
                    context, distObject._get_total_bins()
                )
                expectedInputs = [
                    var
                    for var in explicit_graph_inputs([expression])
                    if var.name is not None
                ]
                expectedFunction = pyhs3_function(
                    inputs=expectedInputs,
                    outputs=expression,
                    on_unused_input="ignore",
                    mode="FAST_RUN",
                )
                self._expectedFunctions[name] = (
                    expectedFunction,
                    [var.name for var in expectedInputs],
                )
            except Exception as e:  # pragma: no cover - unusual graphs
                logger.debug(f"Could not compile expected rates for {name}: {e}")

    # ------------------------------------------------------------------
    # likelihood construction
    # ------------------------------------------------------------------
    def _resolveSigma(self, model, spec: dict) -> float:
        """Return the width of a constraint distribution."""
        sigmaName = spec.get("sigmaName")
        if sigmaName:
            return float(model.nominal_params.get(sigmaName, spec.get("sigma", 1.0)))
        return float(spec.get("sigma", 1.0))

    def _resolveAux(self, model, parameter: str, spec: dict) -> float:
        """Return the auxiliary data value for a constraint distribution.

        For Asimov data the auxiliary measurement is set to the profiled
        nuisance parameter, as in pyhf.
        """
        if parameter in self.auxOverrides:
            return float(self.auxOverrides[parameter])
        auxName = spec.get("auxName")
        if auxName:
            return float(model.nominal_params.get(auxName, spec.get("aux", 0.0)))
        return float(spec.get("aux", 0.0))

    def _gaussianConstraint(self, model, parameter: str, spec: dict):
        """Build the log-probability of a Gaussian (or log-normal) constraint.

        pyhs3 always uses unit widths and zero auxiliary data for its
        ``Gauss`` constraints.  This method uses the width and mean of the
        referenced HS3 constraint distribution instead.
        """
        aux = self._resolveAux(model, parameter, spec)
        sigma = abs(self._resolveSigma(model, spec))
        sigma = max(sigma, 1e-12)
        value = model.parameters[parameter]
        return (
            -0.5 * ((aux - value) / sigma) ** 2
            - pt.log(sigma)
            - 0.5 * math.log(2.0 * math.pi)
        )

    def _staterrorLiteConstraint(self, distObject, model):
        """Build the channel-level Barlow-Beeston lite constraint.

        This mirrors pyhs3's implementation, but allows the auxiliary data to
        be replaced by the profiled nuisance parameters when generating Asimov
        data.
        """
        from pyhs3.distributions.histfactory.modifiers import StatErrorModifier

        staterrorModifier = distObject._find_staterror_modifier()
        if staterrorModifier is None:
            return None
        gammaParams = staterrorModifier.parameters
        constraintType = staterrorModifier.constraint
        totalBins = distObject._get_total_bins()
        totalNominal = np.zeros(totalBins)
        totalVariance = np.zeros(totalBins)
        for sample in distObject.samples:
            hasStaterror = any(
                isinstance(modifier, StatErrorModifier)
                for modifier in sample.modifiers
            )
            if not hasStaterror:
                continue
            totalNominal += np.asarray(sample.data.contents, dtype=float)
            totalVariance += np.square(
                np.asarray(sample.data.errors, dtype=float)
            )
        totalSigma = np.sqrt(totalVariance)

        terms = []
        for index, parameter in enumerate(gammaParams):
            nominal = totalNominal[index]
            sigma = totalSigma[index]
            if nominal <= 0:
                continue
            gamma = model.parameters[parameter]
            if constraintType == "Poisson":
                tau = (nominal / sigma) ** 2
                aux = float(self.auxOverrides.get(parameter, 1.0)) * tau
                mean = gamma * tau
                terms.append(aux * pt.log(mean) - mean - pt.gammaln(aux + 1.0))
            else:
                relError = sigma / nominal
                aux = float(self.auxOverrides.get(parameter, 1.0))
                terms.append(
                    -0.5 * ((aux - gamma) / relError) ** 2
                    - pt.log(relError)
                    - 0.5 * math.log(2.0 * math.pi)
                )
        if not terms:
            return None
        if len(terms) == 1:
            return terms[0]
        return pt.add(*terms)

    def _buildLogProb(self, model):
        """Assemble the joint log-probability with corrected constraints.

        The main Poisson terms are taken from pyhs3.  The constraint terms are
        rebuilt so that non-unit constraint widths and non-zero auxiliary data
        are handled correctly, and so that the Barlow-Beeston lite constraints
        are included in the joint likelihood.
        """
        from pyhs3.distributions.histfactory.modifiers import StatErrorModifier

        context = Context(model.parameters)
        terms = []
        for poisson in model._hfdc_poisson.values():
            terms.append(pt.log(poisson))

        seen = set()
        for distObject in model._distribution_objects:
            if not isinstance(distObject, HistFactoryDistChannel):
                continue
            if distObject.name not in model._hfdc_poisson:
                continue
            lite = getattr(distObject, "barlow_beeston_method", "full") == "lite"
            for dedupKey, modifier, sampleData in distObject.constraint_specs():
                if lite and isinstance(modifier, StatErrorModifier):
                    # The combined channel-level constraint is added below.
                    continue
                if dedupKey is not None:
                    if dedupKey in seen:
                        continue
                    seen.add(dedupKey)
                    spec = self.constraintSpecs.get(dedupKey)
                    if spec is not None and spec.get("type") in ("Gauss", "LogNormal"):
                        terms.append(
                            self._gaussianConstraint(model, dedupKey, spec)
                        )
                        continue
                terms.append(
                    pt.log(modifier.make_constraint(context, sampleData))
                )
            if lite:
                liteConstraint = self._staterrorLiteConstraint(distObject, model)
                if liteConstraint is not None:
                    terms.append(liteConstraint)

        # Constrained normfactors have no modifier that pyhs3 could attach a
        # constraint to, so add their constraint terms here.
        for parameter, spec in self.constraintSpecs.items():
            if not spec.get("isNormFactor") or parameter in seen:
                continue
            if spec.get("type") in ("Gauss", "LogNormal"):
                seen.add(parameter)
                terms.append(self._gaussianConstraint(model, parameter, spec))

        if not terms:
            return pt.constant(0.0)
        if len(terms) == 1:
            return terms[0]
        return pt.add(*terms)

    def __call__(self, parameters: np.ndarray) -> float:
        """Return twice the negative log-likelihood."""
        return self.twiceNLL(parameters)

    def twiceNLL(self, parameters: np.ndarray) -> float:
        """Return twice the negative log-likelihood."""
        try:
            value = self._nllFunction(*np.asarray(parameters, dtype=float))
            value = float(np.asarray(value).reshape(-1)[0])
        except Exception:
            return float("inf")
        if not np.isfinite(value):
            return float("inf")
        return value

    def gradient(self, parameters: np.ndarray) -> Optional[np.ndarray]:
        """Return the gradient of twice the negative log-likelihood."""
        if self._gradFunction is None:
            return None
        try:
            grad = self._gradFunction(*np.asarray(parameters, dtype=float))
            grad = np.asarray(grad, dtype=float).reshape(-1)
        except Exception:
            return None
        if not np.all(np.isfinite(grad)):
            return None
        return grad

    def expectedRates(self, parameters: np.ndarray) -> Dict[str, np.ndarray]:
        """Return the expected rates for all channels at ``parameters``."""
        parameters = np.asarray(parameters, dtype=float)
        values = dict(zip(self.inputNames, parameters))
        result = {}
        for name, (function, inputNames) in self._expectedFunctions.items():
            args = [values[inputName] for inputName in inputNames]
            try:
                rates = np.asarray(function(*args), dtype=float).reshape(-1)
            except Exception:
                rates = np.full(len(self.dataNames.get(name, [])), np.nan)
            result[name] = rates
        return result

    def fit(
        self,
        fixedPoi: Optional[float] = None,
        poiBounds: Optional[tuple] = None,
        allowNegativeSignals: bool = False,
        initial: Optional[np.ndarray] = None,
        maxiter: int = 2000,
    ) -> tuple:
        """Minimise twice the negative log-likelihood.

        :param fixedPoi: if not None, fix the POI to this value and profile
            only the nuisance parameters.
        :param poiBounds: optional bounds for the POI.
        :param allowNegativeSignals: if True, widen the POI bounds to
            ``(-5, 10)``.
        :param initial: initial parameter values.
        :returns: ``(parameters, twice_nll, result)``
        """
        x0 = self.defaults.copy() if initial is None else np.asarray(initial, dtype=float)
        bounds = list(self.bounds)
        if poiBounds is not None:
            bounds[self.poiIndex] = tuple(poiBounds)
        elif allowNegativeSignals:
            bounds[self.poiIndex] = (-5.0, 10.0)
        else:
            bounds[self.poiIndex] = (0.0, 10.0)

        if fixedPoi is not None:
            x0[self.poiIndex] = fixedPoi
            freeIndices = [i for i in range(len(x0)) if i != self.poiIndex]
            if not freeIndices:
                return x0, self.twiceNLL(x0), None

            def objective(values):
                parameters = x0.copy()
                parameters[freeIndices] = values
                parameters[self.poiIndex] = fixedPoi
                return self.twiceNLL(parameters)

            def jacobian(values):
                parameters = x0.copy()
                parameters[freeIndices] = values
                parameters[self.poiIndex] = fixedPoi
                grad = self.gradient(parameters)
                if grad is None:
                    return None
                return grad[freeIndices]

            # Only pass ``jac`` when an analytic gradient is available.
            # A callable that returns ``None`` makes SciPy treat the gradient
            # as zero, which stops the fit at the initial point.
            result = optimize.minimize(
                objective,
                x0[freeIndices],
                jac=(jacobian if self._gradFunction is not None else None),
                method="L-BFGS-B",
                bounds=[bounds[i] for i in freeIndices],
                options={"maxiter": maxiter},
            )
            parameters = x0.copy()
            parameters[freeIndices] = result.x
            parameters[self.poiIndex] = fixedPoi
            return parameters, float(result.fun), result

        # Only pass ``jac`` when an analytic gradient is available.  Passing a
        # callable that returns ``None`` makes SciPy treat the gradient as
        # zero, so the optimiser stops at the initial point.
        result = optimize.minimize(
            self.twiceNLL,
            x0,
            jac=(self.gradient if self._gradFunction is not None else None),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": maxiter},
        )
        return result.x, float(result.fun), result

    def sigmaMu(self, result, parameters: np.ndarray) -> float:
        """Estimate the uncertainty on the POI from the Hessian.

        ``scipy``'s L-BFGS-B result only carries an *approximate* inverse
        Hessian.  If the fit starts at the minimum, that approximation is the
        identity matrix, which would incorrectly suggest ``sigma_mu = 1``.
        We therefore only trust it when it differs from the identity and fall
        back to a numerical Hessian otherwise (as ``pyhfInterface`` does).
        """
        if result is not None and hasattr(result, "hess_inv"):
            try:
                hessInv = result.hess_inv
                if hasattr(hessInv, "todense"):
                    hessInv = hessInv.todense()
                hessInv = np.asarray(hessInv, dtype=float)
                variance = hessInv[self.poiIndex, self.poiIndex]
                if np.isfinite(variance) and variance > 0:
                    sigma = float(np.sqrt(variance))
                    if abs(sigma - 1.0) > 1e-4:
                        return sigma
            except Exception:
                pass

        # Numerical fallback.  Use the full Hessian for small models, and a
        # one-dimensional approximation for larger ones.
        n = len(self.inputNames)
        epsilon = 1e-4
        try:
            if n <= 10:
                hessian = np.zeros((n, n))
                for i in range(n):
                    for j in range(i + 1):
                        stepI = epsilon * np.eye(n)[:, i]
                        stepJ = epsilon * np.eye(n)[:, j]
                        par11 = self.twiceNLL(parameters + stepI + stepJ)
                        par12 = self.twiceNLL(parameters - stepI + stepJ)
                        par21 = self.twiceNLL(parameters + stepI - stepJ)
                        par22 = self.twiceNLL(parameters - stepI - stepJ)
                        value = (par11 - par12 - par21 + par22) / (4 * epsilon**2)
                        hessian[i, j] = value
                        hessian[j, i] = value
                inverse = np.linalg.inv(hessian)
                variance = inverse[self.poiIndex, self.poiIndex]
            else:
                step = epsilon * np.eye(n)[:, self.poiIndex]
                par11 = self.twiceNLL(parameters + 2 * step)
                par12 = self.twiceNLL(parameters)
                par22 = self.twiceNLL(parameters - 2 * step)
                hessian = (par11 - 2 * par12 + par22) / (4 * epsilon**2)
                variance = 1.0 / hessian if hessian != 0 else float("inf")
            if np.isfinite(variance) and variance > 0:
                return float(np.sqrt(variance))
        except Exception:
            pass
        return 1.0


class Hs3UpperLimitComputer:
    """Compute upper limits using an HS3 workspace and pyhs3.

    The class mirrors :class:`~smodels.statistics.pyhfInterface.PyhfUpperLimitComputer`.
    Since pyhs3 does not ship fit or hypothesis-testing helpers, the model is
    compiled with pytensor and fitted with SciPy.  CLs values are obtained
    from the asymptotic formula in
    :func:`~smodels.statistics.basicStats.CLsfromNLL`.
    """

    def __init__(
        self,
        data: Hs3Data,
        cl: float = 0.95,
        lumi: Optional[UnitLumi] = None,
    ):
        self.data = data
        self.lumi = lumi
        self.name = data.jsonFileName
        if data.globalInfo is not None:
            self.name = f"{data.globalInfo.id}:{data.jsonFileName}"
        self.dataType = "hs3"
        self.allowNegativeSignals = False
        self.nsignals = copy.deepcopy(self.data.nsignals)
        self.inputJson = self.data.inputJson
        self.channelsInfo = self.data.channelsInfo
        self.includeCRs = self.data.includeCRs
        self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.cl = cl
        self.scale = 1.0
        self.sigma_mu = None
        self.alreadyBeenThere = False

        self._bundleCache = {}
        self._nllCache = {}
        self._nllMinCache = {}
        self._asimovCache = {}

        self.welcome()

    # ------------------------------------------------------------------
    # bookkeeping
    # ------------------------------------------------------------------
    def welcome(self) -> None:
        """Greet the world."""
        if hs3info["hasgreeted"]:
            return
        logger.info(
            f"HS3 interface, we are using pyhs3 v{hs3info['ver']} with "
            f"pytensor as backend."
        )
        hs3info["hasgreeted"] = True

    def _clearCaches(self) -> None:
        self._bundleCache.clear()
        self._nllCache.clear()
        self._nllMinCache.clear()
        self._asimovCache.clear()

    def _signalsKey(self) -> tuple:
        return tuple(
            (key, round(float(value), 8)) for key, value in sorted(self.nsignals.items())
        )

    def backup(self) -> None:
        self.bu_signal = copy.deepcopy(self.data.nsignals)

    def restore(self) -> None:
        if not hasattr(self, "bu_signal"):
            return
        self.data.nsignals = copy.deepcopy(self.bu_signal)
        del self.bu_signal

    def rescale(self, factor: float) -> None:
        """Rescale the signal predictions and the model.

        This mirrors the pyhf computer.  ``CLs(mu)`` interprets ``mu`` as an
        absolute signal strength, so the rescaling is absorbed by
        ``self.scale``.
        """
        for regionName in self.nsignals.keys():
            self.nsignals[regionName] = self.nsignals[regionName] * factor
        try:
            self.alreadyBeenThere = self.nsignals == self.nsignals_2
        except AttributeError:
            pass
        self.scale *= factor
        logger.debug(f"new signal scale: {self.scale}")
        self._clearCaches()
        try:
            self.nsignals_2 = self.nsignals_1.copy()
        except AttributeError:
            pass
        self.nsignals_1 = self.nsignals.copy()

    def _reset(self) -> None:
        """Reset the signals and the model to the original data."""
        self.nsignals = copy.deepcopy(self.data.nsignals)
        self.scale = 1.0
        self._clearCaches()

    def getTotalXSec(self) -> UnitXSec:
        return self.data.getTotalXSec()

    # ------------------------------------------------------------------
    # workspace/model construction
    # ------------------------------------------------------------------
    def _selectedChannels(self) -> List[str]:
        return [
            name
            for name, info in self.data.channelMap.items()
            if info.get("included", False)
        ]

    def _findDataForChannel(self, js: dict, channel: str) -> Optional[str]:
        """Find a data entry that matches the channel's axis names."""
        dist = self.data.distributions[channel]
        axisNames = tuple(axis.get("name") for axis in dist.get("axes", []))
        for datum in js.get("data", []):
            if datum.get("type") != "binned":
                continue
            datumAxes = tuple(axis.get("name") for axis in datum.get("axes", []))
            if datumAxes == axisNames:
                return datum.get("name")
        for datum in js.get("data", []):
            if datum.get("type") == "binned" and channel in datum.get("name", ""):
                return datum.get("name")
        return None

    def _chooseLikelihood(self, js: dict, selected: List[str]) -> dict:
        """Choose the likelihood that covers most selected channels."""
        selectedSet = set(selected)

        def coverage(likelihood):
            count = 0
            for distRef in likelihood.get("distributions", []):
                name = distRef if isinstance(distRef, str) else distRef.get("name")
                if name in selectedSet:
                    count += 1
            return count

        def observedScore(likelihood):
            text = str(likelihood.get("name", "")).lower()
            for dataRef in likelihood.get("data", []):
                name = dataRef if isinstance(dataRef, str) else dataRef.get("name", "")
                text += " " + str(name).lower()
            return 0 if "asimov" in text and "obs" not in text else 1

        likelihoods = js.get("likelihoods", [])
        if likelihoods:
            best = max(
                likelihoods,
                key=lambda likelihood: (coverage(likelihood), observedScore(likelihood)),
            )
            if coverage(best) > 0:
                return best
        return {"name": "smodels_likelihood", "distributions": [], "data": []}

    def _chooseAnalysis(self, js: dict, likelihood: dict) -> dict:
        """Choose an analysis that uses the selected likelihood."""
        analyses = [
            analysis
            for analysis in js.get("analyses", [])
            if analysis.get("likelihood") == likelihood.get("name")
        ]
        if not analyses:
            analyses = js.get("analyses", [])

        def observedScore(analysis):
            text = str(analysis.get("name", "")).lower()
            return 0 if "asimov" in text and "obs" not in text else 1

        if analyses:
            return max(analyses, key=observedScore)
        return {
            "name": "smodels_analysis",
            "likelihood": likelihood.get("name"),
            "domains": [],
            "init": None,
        }

    def _prepareJson(
        self,
        dataContents: Optional[Dict[str, List[float]]] = None,
        aprioriData: bool = False,
    ) -> tuple:
        """Create a small, patched HS3 workspace for the selected channels."""
        js = copy.deepcopy(self.inputJson)
        constraintSpecs = _normalizeConstraints(js)

        selected = self._selectedChannels()
        if not selected:
            raise SModelSError("No HS3 channels selected")

        # Keep only the selected HistFactory channels.  Constraint
        # distributions are not needed after _normalizeConstraints().
        distributions = {
            dist.get("name"): dist
            for dist in js.get("distributions", [])
            if dist.get("type") == "histfactory_dist"
            and dist.get("name") in selected
        }
        missing = [name for name in selected if name not in distributions]
        if missing:
            raise SModelSError(f"HS3 channels {missing} not found in workspace")
        js["distributions"] = [distributions[name] for name in selected]

        # Choose a likelihood and build the distribution/data pairs.
        likelihood = self._chooseLikelihood(js, selected)
        dataByName = {datum.get("name"): datum for datum in js.get("data", [])}
        pairs = []
        pairedChannels = set()
        for distRef, dataRef in zip(
            likelihood.get("distributions", []), likelihood.get("data", [])
        ):
            distName = distRef if isinstance(distRef, str) else distRef.get("name")
            dataName = dataRef if isinstance(dataRef, str) else dataRef.get("name")
            if distName in selected:
                if dataName not in dataByName:
                    dataName = self._findDataForChannel(js, distName)
                pairs.append((distName, dataName))
                pairedChannels.add(distName)
        for channel in selected:
            if channel in pairedChannels:
                continue
            dataName = self._findDataForChannel(js, channel)
            if dataName is None:
                raise SModelSError(
                    f"Could not find binned data for HS3 channel {channel}"
                )
            pairs.append((channel, dataName))

        # Filter data and the likelihood to the selected pairs.
        dataNames = []
        for _, dataName in pairs:
            if dataName not in dataNames:
                dataNames.append(dataName)
        js["data"] = [dataByName[name] for name in dataNames if name in dataByName]
        dataByName = {datum.get("name"): datum for datum in js["data"]}
        likelihood = dict(likelihood)
        likelihood["name"] = likelihood.get("name", "smodels_likelihood")
        likelihood["distributions"] = [distName for distName, _ in pairs]
        likelihood["data"] = [dataName for _, dataName in pairs]
        js["likelihoods"] = [likelihood]

        # Add the signal sample and (optionally) its uncertainty.
        for channel in selected:
            dist = distributions[channel]
            # ``rescale`` scales the model signal, while CLs() converts the
            # absolute mu back to a value relative to this scale.  This keeps
            # the physics unchanged and only improves the numerical range.
            signals = [
                float(x) * self.scale
                for x in self.data.channelMap[channel]["signals"]
            ]
            # SModelS supplies its own signal sample.  Drop an existing
            # SModelS-style signal sample if the workspace already has one,
            # to avoid counting the signal twice.
            samples = []
            for sample in dist.get("samples", []):
                if sample.get("name") == "bsm":
                    continue
                sampleParameters = _modifierParameters({"samples": [sample]})
                if poiName in sampleParameters:
                    continue
                samples.append(sample)
            modifiers = [
                {"name": poiName, "parameter": poiName, "type": "normfactor"}
            ]
            if (
                self.data.signalUncertainty is not None
                and self.data.signalUncertainty != 0
            ):
                uncertainty = float(self.data.signalUncertainty)
                modifiers.append(
                    {
                        "name": signalUncertaintyName,
                        "parameter": signalUncertaintyName,
                        "type": "histosys",
                        "data": {
                            "hi": {
                                "contents": [
                                    x * (1.0 + uncertainty) for x in signals
                                ]
                            },
                            "lo": {
                                "contents": [
                                    x * (1.0 - uncertainty) for x in signals
                                ]
                            },
                        },
                    }
                )
            samples.append(
                {
                    "name": "bsm",
                    "data": {"contents": signals, "errors": [0.0] * len(signals)},
                    "modifiers": modifiers,
                }
            )
            dist["samples"] = samples

        # Replace the data contents for expected or Asimov data.
        if aprioriData or dataContents is not None:
            contentsByDataName = {}
            for channel, dataName in pairs:
                datum = dataByName.get(dataName)
                if datum is None:
                    raise SModelSError(
                        f"HS3 data {dataName} for channel {channel} not found"
                    )
                if dataContents is not None:
                    contents = dataContents.get(channel)
                    if contents is None:
                        raise SModelSError(
                            f"No Asimov data provided for HS3 channel {channel}"
                        )
                else:
                    contents = np.zeros(_nbins(distributions[channel]))
                    for sample in distributions[channel].get("samples", []):
                        if sample.get("name") == "bsm":
                            continue
                        if poiName in _modifierParameters({"samples": [sample]}):
                            continue
                        contents += np.asarray(
                            sample.get("data", {}).get("contents", []), dtype=float
                        )
                contents = np.asarray(contents, dtype=float).reshape(-1)
                if len(contents) != len(datum.get("contents", [])):
                    raise SModelSError(
                        f"HS3 data for channel {channel} has "
                        f"{len(datum.get('contents', []))} bins, expected "
                        f"{len(contents)}"
                    )
                contentsByDataName[dataName] = contents
            for dataName, contents in contentsByDataName.items():
                dataByName[dataName]["contents"] = [float(x) for x in contents]

        # Choose the analysis and make sure it refers to our likelihood.
        analysis = self._chooseAnalysis(js, likelihood)
        analysis = dict(analysis)
        analysis["likelihood"] = likelihood["name"]
        js["analyses"] = [analysis]

        # Pick or create a parameter point.
        parameterPoints = {
            point.get("name"): point for point in js.get("parameter_points", [])
        }
        initName = analysis.get("init")
        if initName not in parameterPoints:
            if parameterPoints:
                initName = next(iter(parameterPoints))
            else:
                initName = "smodels_parameters"
                parameterPoints[initName] = {"name": initName, "parameters": []}
        parameterPoint = parameterPoints[initName]
        parameters = {
            parameter.get("name"): parameter
            for parameter in parameterPoint.get("parameters", [])
        }

        # Make sure every parameter used by the selected channels exists.
        required = [poiName]
        for channel in selected:
            required.extend(_modifierParameters(distributions[channel]))
        if signalUncertaintyName in required:
            required = [name for name in required if name != signalUncertaintyName]
            required.append(signalUncertaintyName)
        for name in required:
            if name not in parameters:
                parameters[name] = {
                    "name": name,
                    "value": _defaultParameterValue(name),
                }
        parameterPoint["parameters"] = list(parameters.values())
        js["parameter_points"] = [parameterPoint]
        analysis["init"] = initName

        # Pick or create a product domain and make sure it has all axes.
        domains = {
            domain.get("name"): domain for domain in js.get("domains", [])
        }
        analysisDomains = analysis.get("domains", [])
        targetDomain = None
        for domainName in analysisDomains:
            domain = domains.get(domainName)
            if domain is not None and domain.get("type") == "product_domain":
                targetDomain = domain
                break
        if targetDomain is None:
            for domain in domains.values():
                if domain.get("type") == "product_domain":
                    targetDomain = domain
                    break
        if targetDomain is None:
            targetDomain = {
                "name": "smodels_domain",
                "type": "product_domain",
                "axes": [],
            }
            domains[targetDomain["name"]] = targetDomain
        axes = {
            axis.get("name"): axis for axis in targetDomain.get("axes", [])
        }
        for name in required:
            if name in axes:
                continue
            lower, upper = _defaultParameterBounds(name)
            axes[name] = {
                "name": name,
                "min": lower,
                "max": upper,
                "nbins": 1000,
            }
        targetDomain["axes"] = list(axes.values())
        analysis["domains"] = [targetDomain["name"]]
        js["domains"] = [targetDomain]

        return js, constraintSpecs

    def _patchStatErrorData(self, workspace) -> None:
        """Allow pyhs3 to build models with Barlow-Beeston lite staterror.

        pyhs3's lite mode stores no ``data`` on the ``staterror`` modifier, but
        its model builder still calls ``make_constraint`` for every modifier.
        Supplying the per-sample uncertainties makes model construction
        succeed.  The resulting per-sample constraints are ignored in favour of
        the proper channel-level lite constraint (see ``_Hs3ModelBundle``).
        """
        from pyhs3.distributions.histfactory.modifiers import (
            StatErrorData,
            StatErrorModifier,
        )

        for distObject in workspace.distributions or []:
            if not isinstance(distObject, HistFactoryDistChannel):
                continue
            if getattr(distObject, "barlow_beeston_method", "lite") != "lite":
                continue
            for sample in distObject.samples:
                errors = list(getattr(sample.data, "errors", []) or [])
                for modifier in sample.modifiers:
                    if not isinstance(modifier, StatErrorModifier):
                        continue
                    if modifier.data is not None:
                        continue
                    nParams = len(modifier.parameters)
                    uncertainties = errors if len(errors) == nParams else [1.0] * nParams
                    modifier.data = StatErrorData(
                        uncertainties=[float(x) for x in uncertainties]
                    )

    def _buildBundle(
        self,
        dataContents: Optional[Dict[str, List[float]]] = None,
        aprioriData: bool = False,
        auxOverrides: Optional[Dict[str, float]] = None,
    ) -> Optional[_Hs3ModelBundle]:
        """Build a pyhs3 model for the current signals and data."""
        js = None
        try:
            js, constraintSpecs = self._prepareJson(
                dataContents=dataContents, aprioriData=aprioriData
            )
            workspace = pyhs3.Workspace(**js)
            self._patchStatErrorData(workspace)
            analysis = js["analyses"][0]
            model = workspace.model(
                analysis["name"],
                parameter_set=analysis.get("init"),
                progress=False,
            )
            dataNames = {
                distName: dataName
                for distName, dataName in zip(
                    js["likelihoods"][0]["distributions"],
                    js["likelihoods"][0]["data"],
                )
            }
            return _Hs3ModelBundle(
                model,
                workspace=workspace,
                workspaceJson=js,
                dataNames=dataNames,
                constraintSpecs=constraintSpecs,
                auxOverrides=auxOverrides,
            )
        except Exception as e:
            logger.error(f"Could not construct HS3 model: {e}")
            if logger.level <= 10:
                logger.debug("HS3 workspace was:\n%s", js)
            return None

    def _baseBundle(self, evaluationType: NllEvalType) -> Optional[_Hs3ModelBundle]:
        """Return the observed or a-priori bundle."""
        if evaluationType == apriori:
            mode = "apriori"
        else:
            mode = "observed"
        key = (mode, self.scale, self._signalsKey())
        if key in self._bundleCache:
            return self._bundleCache[key]
        bundle = self._buildBundle(aprioriData=(mode == "apriori"))
        self._bundleCache[key] = bundle
        return bundle

    def _asimovBundle(
        self, mu: float, evaluationType: NllEvalType
    ) -> Optional[_Hs3ModelBundle]:
        """Return a bundle whose data are the Asimov data for signal ``mu``."""
        eT = observed if evaluationType == aposteriori else evaluationType
        key = (round(float(mu), 8), eT, self.scale, self._signalsKey())
        if key in self._asimovCache:
            return self._asimovCache[key]

        base = self._baseBundle(eT)
        if base is None:
            return None
        try:
            parameters, _, _ = base.fit(fixedPoi=float(mu), poiBounds=(-5.0, 10.0))
        except Exception as e:
            logger.info(f"HS3 fixed-POI fit for Asimov data failed: {e}")
            parameters = base.defaults.copy()
            parameters[base.poiIndex] = float(mu)
        contents = base.expectedRates(parameters)
        # As in pyhf, the auxiliary measurements of the Asimov dataset are set
        # to the profiled nuisance parameters.
        auxOverrides = {
            name: float(value)
            for name, value in zip(base.inputNames, parameters)
        }
        bundle = self._buildBundle(
            dataContents=contents, auxOverrides=auxOverrides
        )
        self._asimovCache[key] = bundle
        return bundle

    def updateWorkspace(self, evaluationType: NllEvalType = observed):
        """Return the pyhs3 workspace used for the given evaluation type."""
        bundle = self._baseBundle(evaluationType)
        if bundle is None:
            return None
        return bundle.workspace

    def generateAsimovData(
        self, mu: Union[None, float] = 0.0, evaluationType: NllEvalType = observed
    ) -> Optional[_Hs3ModelBundle]:
        """Generate Asimov data for the model at signal strength ``mu``.

        If ``mu`` is None, the observed (or a-priori) bundle is returned
        unchanged.
        """
        if mu is None:
            return self._baseBundle(evaluationType)
        return self._asimovBundle(float(mu), evaluationType)

    # ------------------------------------------------------------------
    # likelihood interface
    # ------------------------------------------------------------------
    def nll(
        self,
        mu: float = 1.0,
        evaluationType: NllEvalType = observed,
        asimov: Optional[float] = None,
    ) -> Optional[float]:
        """Return the negative log-likelihood with the POI fixed to ``mu``.

        :param evaluationType: one of observed, apriori, aposteriori
        :param asimov: if not None, generate Asimov data for this signal
            strength and evaluate the likelihood on it.
        """
        key = (round(float(mu), 8), evaluationType, asimov, self.scale,
               self._signalsKey())
        if key in self._nllCache:
            return self._nllCache[key]

        if asimov is None:
            bundle = self._baseBundle(evaluationType)
        else:
            bundle = self._asimovBundle(float(asimov) / self.scale, evaluationType)
        if bundle is None:
            return None
        try:
            _, twiceNLL, _ = bundle.fit(
                fixedPoi=float(mu) / self.scale, poiBounds=(-5.0, 10.0)
            )
        except Exception as e:
            logger.info(f"HS3 fixed-POI fit failed for mu={mu:.3f}: {e}")
            return None
        result = twiceNLL / 2.0
        self._nllCache[key] = result
        return result

    def nll_min(
        self,
        evaluationType: NllEvalType = observed,
        allowNegativeSignals: bool = False,
    ) -> Optional[dict]:
        """Return the minimum negative log-likelihood and the fitted POI."""
        key = (evaluationType, allowNegativeSignals, self.scale, self._signalsKey())
        if key in self._nllMinCache:
            return copy.deepcopy(self._nllMinCache[key])

        bundle = self._baseBundle(evaluationType)
        if bundle is None:
            return None
        try:
            parameters, twiceNLL, result = bundle.fit(
                allowNegativeSignals=allowNegativeSignals
            )
        except Exception as e:
            logger.info(f"HS3 unconditional fit failed: {e}")
            return None
        muhat = float(parameters[bundle.poiIndex]) * self.scale
        sigmaMu = bundle.sigmaMu(result, parameters) * self.scale
        ret = {
            "muhat": muhat,
            "sigma_mu": sigmaMu,
            "nll_min": twiceNLL / 2.0,
        }
        self._nllMinCache[key] = ret
        return copy.deepcopy(ret)

    def clsFromNLLs(
        self, mu: float, evaluationType: NllEvalType, return_type: str
    ) -> Optional[float]:
        """Compute CLs from the NLL values, mirroring the pyhf computer."""
        ret = self.nll_min(evaluationType=evaluationType)
        if ret is None:
            return None
        nll_min = ret["nll_min"]
        muhat = ret["muhat"]
        nll = self.nll(mu=mu, evaluationType=evaluationType)
        eT = observed if evaluationType == aposteriori else evaluationType
        nllA = self.nll(mu=mu, evaluationType=eT, asimov=0.0)
        nll_minA = self.nll(mu=0.0, evaluationType=eT, asimov=0.0)
        if nll is None or nllA is None or nll_minA is None:
            return None
        if evaluationType == aposteriori:
            nll = nllA
            nll_min = nll_minA
        result = CLsfromNLL(
            nllA,
            nll_minA,
            nll,
            nll_min,
            muhat > mu,
            return_type,
            report_all=True,
        )
        return result["returns"]

    def CLs(
        self,
        mu: float,
        evaluationType: NllEvalType,
        return_type: Text = "CLs",
        nSigma: int = 0,
        reset_data: bool = True,
        **kwargs,
    ) -> Optional[float]:
        """Compute CLs for an absolute signal strength ``mu``."""
        if reset_data:
            self._reset()
        if "pmSigma" in kwargs:
            assert kwargs["pmSigma"] == 0, "no CLs with pmSigma for HS3"
        return self._CLs(mu / self.scale, evaluationType, return_type, nSigma)

    def _CLs(
        self,
        mu_rel: float,
        evaluationType: NllEvalType,
        return_type: Text = "CLs",
        nSigma: int = 0,
    ) -> Optional[float]:
        """Compute CLs for a signal strength relative to the current scale."""
        assert return_type in ["CLs-alpha", "alpha-CLs", "1-CLs", "CLs"], (
            f"Unknown return type: {return_type}."
        )
        if self.data.errorFlag or self.zeroSignalsFlag:
            return None
        if not np.isfinite(mu_rel):
            return None

        ret = self.nll_min(evaluationType=evaluationType)
        if ret is None:
            return None
        nll_min = ret["nll_min"]
        muhat = ret["muhat"]

        mu = mu_rel * self.scale
        nll = self.nll(mu=mu, evaluationType=evaluationType)
        eT = observed if evaluationType == aposteriori else evaluationType
        nllA = self.nll(mu=mu, evaluationType=eT, asimov=0.0)
        nll_minA = self.nll(mu=0.0, evaluationType=eT, asimov=0.0)
        if nll is None or nllA is None or nll_minA is None:
            return None
        if evaluationType == aposteriori:
            nll = nllA
            nll_min = nll_minA
        try:
            # Reproduce pyhf's asymptotic CLs (arXiv:1007.1727) in the
            # sqrt(qmu) test-statistic space.
            qmu = 0.0 if (nll < nll_min or muhat > mu) else 2.0 * (nll - nll_min)
            qA = max(0.0, 2.0 * (nllA - nll_minA))
            sqrtQmu = np.sqrt(qmu)
            sqrtQA = np.sqrt(qA)
            if nSigma != 0:
                # Expected band.  pyhf orders its expected set from +2sigma
                # to -2sigma, while SModelS uses the opposite convention for
                # the requested band, hence the flipped sign.
                ns = -nSigma
                clsb = 1.0 - norm.cdf(ns + sqrtQA)
                clb = norm.cdf(-ns)
            else:
                if sqrtQmu <= sqrtQA or sqrtQA == 0.0:
                    testStat = sqrtQmu - sqrtQA
                else:
                    testStat = (qmu - qA) / (2.0 * sqrtQA)
                clsb = 1.0 - norm.cdf(testStat + sqrtQA)
                clb = norm.cdf(-testStat)
            cl = clsb / clb if clb > 0.0 else 0.0
            return float(clsType(cl, return_type, self.cl))
        except Exception as e:
            logger.error(f"CLs computation failed for mu={mu_rel:.3f}: {e}")
            return None

    def getSign(
        self, lo_mu: float, hi_mu: float, evaluationType: NllEvalType
    ) -> float:
        """Return CLs(lo_mu) * CLs(hi_mu)."""
        lo = self._CLs(lo_mu, evaluationType, "CLs-alpha")
        hi = self._CLs(hi_mu, evaluationType, "CLs-alpha")
        if lo is None or hi is None:
            return float("nan")
        return lo * hi

    def getUpperLimitOnMu(
        self, evaluationType: NllEvalType = observed, nSigma: int = 0, **kwargs
    ) -> Optional[float]:
        """Compute the upper limit on the signal strength modifier.

        :param evaluationType: if apriori, use the expected SM background
            as observations.
        :param nSigma: compute the limit for this many standard deviations
            away from the central expected limit.
        :returns: the upper limit at ``self.cl`` confidence level.
        """
        if "pmSigma" in kwargs:
            kwargs.pop("pmSigma")
        if len(kwargs) > 0:
            logger.warning(f"HS3 computer will ignore {kwargs}")
        return self._getUpperLimitOnMu(evaluationType, nSigma)

    def _getUpperLimitOnMu(
        self, evaluationType: NllEvalType = observed, nSigma: int = 0
    ) -> Optional[float]:
        self._reset()
        if self.data.errorFlag:
            return None
        if self.zeroSignalsFlag:
            logger.debug(f"no signals were found for {self.data.jsonFileName}")
            return None

        startUL = time.time()
        lo_mu = 0.0
        hi_mu = 1.0

        f_lo = self.CLs(lo_mu, evaluationType, "alpha-CLs", nSigma,
                        reset_data=False)
        if f_lo is None:
            return None
        if f_lo >= 0.0:
            return 0.0

        f_hi = None
        for attempt in range(nattempts_max + 1):
            f_hi = self.CLs(hi_mu, evaluationType, "alpha-CLs", nSigma,
                            reset_data=False)
            if f_hi is None:
                hi_mu *= 2.0
                continue
            if f_hi > 0.0:
                break
            lo_mu = hi_mu
            hi_mu *= 2.0
        else:
            logger.warning("could not bracket the root for getUpperLimitOnMu")
            return None

        if f_hi is None or not np.isfinite(f_hi):
            return None
        try:
            ul = findRoot(
                self.CLs,
                lo_mu,
                hi_mu,
                args=(evaluationType, "alpha-CLs", nSigma, False),
                rtol=1e-3,
                xtol=1e-3,
            )
        except Exception as e:
            logger.warning(f"root finding for HS3 upper limit failed: {e}")
            return None

        endUL = time.time()
        logger.debug(
            f"getUpperLimitOnMu elapsed time : {endUL - startUL:1.4f} secs"
        )
        return float(ul)

    def totalYieldsFromSignals(self, poi_test: float) -> List[float]:
        """Return the total yields (signal + background) for debugging."""
        bundle = self._baseBundle(observed)
        if bundle is None:
            return []
        parameters = bundle.defaults.copy()
        parameters[bundle.poiIndex] = float(poi_test)
        expected = bundle.expectedRates(parameters)
        yields = []
        for channel in self._selectedChannels():
            yields.extend([float(x) for x in expected.get(channel, [])])
        return yields


__all__ = ["Hs3Data", "Hs3UpperLimitComputer"]

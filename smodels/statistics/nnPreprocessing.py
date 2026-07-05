#!/usr/bin/env python3
"""
.. module:: preprocessing_nnAdapter
   :synopsis: Preprocessing and inverse-preprocessing functions for the
   NNAdapter. Handles feature standardisation and nLL postprocessing,
   including uncertainty propagation through arbitrary transformation chains.

"""

__all__ = [
    "preprocess_features",
    "postprocess_nLLs",
    "postprocess_nLLs_errors",
]

import numpy as np
from typing import Optional


# ---------------------------------------------------------------------------
# Elementary transforms (forward and inverse)
# ---------------------------------------------------------------------------

def _log_with_negatives(x: np.ndarray) -> np.ndarray:
    """Signed log1p transform, works for negative values, exactly invertible."""
    return np.sign(x) * np.log1p(np.abs(x))


def _undo_log_with_negatives(x: np.ndarray) -> np.ndarray:
    """Exact inverse of _log_with_negatives."""
    return np.sign(x) * np.expm1(np.abs(x))


def _standardize(
    x: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
) -> np.ndarray:
    """Apply standardisation using pre-computed mean and std."""
    return (x - mean) / std


def _undo_standardize(
    x: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
) -> np.ndarray:
    """Invert standardisation."""
    return x * std + mean


# ---------------------------------------------------------------------------
# Feature preprocessing (forward, used at inference time)
# ---------------------------------------------------------------------------

def _get_fn(fn_str: str):
    """Return a named scalar transform."""
    match fn_str:
        case "log":
            return np.log
        case "exp":
            return np.exp
        case "sqrt":
            return np.sqrt
        case "inverse":
            return lambda x: 1.0 / x
        case "log_w_negatives":
            return _log_with_negatives
        case _:
            raise ValueError(f"Unknown transform '{fn_str}'.")


def preprocess_features(
    features_raw: np.ndarray,
    trafos: Optional[dict] = None,
    mean: Optional[np.ndarray] = None,
    std: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """Apply feature preprocessing using saved mean/std from training.

    :param features_raw: raw input features, shape (batch, n_features)
    :param trafos: dict mapping set-name to list of transform strings,
        e.g. ``{"fvs_standardized": ["log_w_negatives", "standardization"]}``
    :param mean: per-feature mean saved at training time (required when
        "standardization" is in trafos)
    :param std: per-feature std saved at training time (required when
        "standardization" is in trafos)
    :returns: (scaled_features, mean, std)
    """
    assert np.isfinite(features_raw).all(), "Non-finite values in features_raw"

    feature_sets = {}
    if trafos:
        for set_name, fn_list in trafos.items():
            transformed = features_raw.copy()
            for fn_str in fn_list:
                if fn_str == "standardization":
                    assert mean is not None and std is not None, \
                        "mean and std must be provided for standardization at inference time"
                    transformed = _standardize(transformed, mean, std)
                else:
                    transformed = _get_fn(fn_str)(transformed)
            feature_sets[set_name] = transformed

    scaled = list ( map ( float, feature_sets["fvs_standardized"] ) )
    # scaled = np.concatenate(list(feature_sets.values()), axis=1) if feature_sets else features_raw
    return scaled, mean, std


# ---------------------------------------------------------------------------
# nLL postprocessing (inverse transforms)
# ---------------------------------------------------------------------------

def _get_inv_fn(fn_str: str):
    """Return the inverse of a named scalar transform."""
    match fn_str:
        case "log":
            return np.exp
        case "exp":
            return np.log
        case "sqrt":
            return lambda x: x ** 2
        case "inverse":
            return lambda x: 1.0 / x
        case "log_w_negatives":
            return _undo_log_with_negatives
        case _:
            raise ValueError(f"No known inverse for transform '{fn_str}'.")


def postprocess_nLLs(
    nLLs: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    trafos: Optional[list] = None,
) -> np.ndarray:
    """Undo nLL preprocessing in reverse order.

    :param nLLs: preprocessed nLL deltas output by the NN, shape (4,)
    :param mean: per-output mean saved at training time
    :param std: per-output std saved at training time
    :param trafos: list of transform strings applied during training,
        e.g. ``["log_w_negatives", "standardization"]``
    :returns: unpreprocessed nLL deltas, same shape as nLLs
    """
    if trafos:
        for fn_str in reversed(trafos):
            if fn_str == "standardization":
                nLLs = _undo_standardize(nLLs, mean, std)
            else:
                nLLs = _get_inv_fn(fn_str)(nLLs)
    assert np.isfinite(nLLs).all(), "Non-finite values after postprocess_nLLs"
    return nLLs


def postprocess_nLLs_errors(
    errors: np.ndarray,
    nLLs_prepd: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    trafos: Optional[list] = None,
    eps: float = 1e-5,
) -> np.ndarray:
    """Propagate heteroskedastic errors through the inverse nLL preprocessing.

    Uses central-difference numerical differentiation to compute
    ``|d(postprocess)/d(nLL)| * sigma`` for each output independently.

    :param errors: NN-predicted uncertainties on the preprocessed deltas,
        shape (4,), these are the errors on nLLs_prepd
    :param nLLs_prepd: the central (mean) preprocessed nLL deltas, shape (4,)
    :param mean: per-output mean saved at training time
    :param std: per-output std saved at training time
    :param trafos: same transform list used in postprocess_nLLs
    :param eps: step size for numerical differentiation
    :returns: propagated uncertainties on the unpreprocessed nLL deltas,
        shape (4,)
    """
    f_plus  = postprocess_nLLs(nLLs_prepd + eps, mean, std, trafos)
    f_minus = postprocess_nLLs(nLLs_prepd - eps, mean, std, trafos)
    deriv = np.abs(f_plus - f_minus) / (2.0 * eps)
    return deriv * errors

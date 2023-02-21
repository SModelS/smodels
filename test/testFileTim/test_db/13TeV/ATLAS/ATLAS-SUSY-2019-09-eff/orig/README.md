# HistFactory Likelihoods for the ATLAS SUSY Electroweak
The HistFactory likelihoods of the three-lepton searches are serialised in JSON.
The likelihoods can be provided as input or manipulated with [pyhf](https://iris-hep.org/projects/pyhf.html).

For each analysis channel (on-shell, off-shell) a background-only workspace is provided 
(bkg_onshell.json, bkg_offshell.json). Information relative to the signal models
is contained in files called 'patchsets'. Patchsets are provided for
each analysis channel and for each SUSY Electroweak scanario considered:
wino/bino (+) for the on-shell channel and wino/bino (+), wino/bino (-) and higgsino for
the off-shell channel. 

### Performing a fit

In order to perform a fit, a background workspace has to be applied a patch.
For example, if the desired signal point is defined by m($\chi^{\pm} _{1}$) = 250 GeV and m($\chi^{0} _{1}$) = 50 GeV, the likelihood including the on-shell regions is built with:
```
pyhf patchset apply --name CN_WZ_250_50_harmonised_onshell bkg_onshell.json onshell_winobino_plus_patchset.json --output-file CN_WZ_250_50_onshell_from_patch.json
```
The string following the --name option is the patch name. 
All patches are called follwoing the pattern:
```
CN_WZ_{chargino mass}_{neutralino mass}_harmonised_{channel = onshell or offshell}
```
The adjective 'harmonised' appearing in the patch name remarks the fact that the onshell and offshell likelihoods share the same conventions for the naming of parameters.

The exclusion fit can be run with:
```
pyhf cls CN_WZ_250_50_onshell_from_patch.json
```

### Combined fits

For the wino/bino (+) scenario, the signal points with  m($\chi^{\pm} _{1}$)  -  m($\chi^{0} _{1}$) between 78 and 108 GeV are in common between the on-shell and the off-shell channels, allowing for a combined fit.

To permorm a combined fit:

* create the workspaces using the patchsets:
```
pyhf patchset apply --name CN_WZ_200_100_harmonised_onshell bkg_onshell.json onshell_winobino_plus_patchset.json --output-file CN_WZ_200_100_onshell_from_patch.json
```
```
pyhf patchset apply --name CN_WZ_200_100_harmonised_offshell bkg_offshell.json offshell_winobino_plus_patchset.json --output-file CN_WZ_200_100_offshell_from_patch.json
```
* combine the workspaces
```
pyhf combine --join outer CN_WZ_200_100_onshell_from_patch.json CN_WZ_200_100_offshell_from_patch.json --output-file combined_CN_WZ_200_100.json
```
* run the fit
```
pyhf cls combined_CN_WZ_200_100.json 
```

### Combination with the 2L compressed workspaces

The likelihoods of the [2L compressed analysis](https://doi.org/10.1103/PhysRevD.101.052005) can be downloaded at https://www.hepdata.net/record/ins1767649.
To perform a combination the parameters have to be renamed accoring to the mapping specified in `compressed_mapping.py`.

### Additional Notes

`pyhf` v0.6.2+ is recommended to use these patchsets. See [scikit-hep/pyhf#1488](https://github.com/scikit-hep/pyhf/pull/1488) for more details.

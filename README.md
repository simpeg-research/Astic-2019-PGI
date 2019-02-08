**[summary](#summary) | [contents](#contents) | [usage](#usage) | [running the notebooks](#running-the-notebooks) | [issues](#issues) | [citation](#citation) | [license](#license)**

# Open source software for petrophysically and geologically guided geophysical inversion

[![Build Status](https://travis-ci.org/simpeg-research/Astic-2019-PGI.svg?branch=master)](https://travis-ci.org/simpeg-research/Astic-2019-PGI)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/simpeg-research/Astic-2019-PGI/master?filepath=index.ipynb)
[![Azure](https://notebooks.azure.com/launch.png)](https://notebooks.azure.com/import/gh/simpeg-research/Astic-2019-PGI)
[![Zenodo](https://zenodo.org/badge/124603211.svg)](https://zenodo.org/badge/latestdoi/124603211)
[![License](https://img.shields.io/github/license/simpeg-research/Astic-2019-PGI.svg)](https://github.com/simpeg-research/Astic-2019-PGI/blob/master/LICENSE)
[![SimPEG](https://img.shields.io/badge/powered%20by-SimPEG-blue.svg)](http://simpeg.xyz)

Notebooks and python scripts to reproduce the figures shown in Astic & Oldenburg (2019).

<img src="figures/currents.png" width=40% align="middle">

## Summary

We propose a new framework for incorporating petrophysical and geological information into voxel-based geophysical inversion. By quantitatively linking these data into a single framework, we recover a final inverted model that reproduces the observed, or desired, petrophysical and geological features while fitting the geophysical data. In the hope of facilitating this exploration and promoting reproducibility of geophysical simulations and inversions, we make the scripts of those examples publicly available thanks to the open source software package, SimPEG. This allows researchers to interrogate all of the components and to facilitate the exploration of our new inversion strategies. We highlight different capabilities of our methodology by inverting magnetotelluric and DC resistivity data in 1D and 2D respectively. Finally we apply our framework to inverting airborne over frequency domain data, acquired in Australia, for the detection and characterization of saline contamination of freshwater.

## Contents

There are 4 notebooks in this repository:

- [1_TEM_VerticalConductor_2D_forward.ipynb](/notebooks/1_TEM_VerticalConductor_2D_forward.ipynb) : runs a forward simulation of an airborne electromagnetic simulation over a conductive plate. This notebook was used to generate figures 1-4 in the abstract
- [2_TEM_VerticalConductor_1D_stitched_inversion.ipynb](/notebooks/2_TEM_VerticalConductor_1D_stitched_inversion.ipynb) : Using the forward simulated data from the previous notebook, we run 1D inversions over the plate (Figure 5 in the abstract).
- [3_TEM_VerticalConductor_2D_inversion_load.ipynb](/notebooks/3_TEM_VerticalConductor_2D_inversion_load.ipynb) : This notebook loads the 2D inversion results over the plate (Figure 6 in the abstract). The 2D inversion was run using the script [2dinv_smooth.py](/notebooks/2d_inv_smooth/2dinv_smooth.py).
- [4_TEM_VerticalConductor_parametric_inversion_load.ipynb](/notebooks/4_TEM_VerticalConductor_parametric_inversion_load.ipynb) : This notebook loads the 2D parametric inversion inversion results (Figure 7 in the abstract). The 2D parametric inversion was run using the script [2dinv_parametric.py](/notebooks/2d_inv_parametric/2d_inv_parametric.py) .

In addition, there are two notebooks used for demos in the workshop [3D EM Modelling and Inversion with Open Source Resources](https://courses.geosci.xyz/aem2018):

- [TEM_VerticalConductor_2D_forward.ipynb](/demo_notebooks/TEM_VerticalConductor_2D_forward.ipynb) : runs a forward simulation of an airborne electromagnetic simulation over a conductive plate. Similar to that in the notebooks directory.
- [TDEM_1D_inversion.ipynb](/demo_notebooks/TDEM_1D_inversion.ipynb): In this notebook, we run a 1D inversion for a single airborne time domain EM sounding

## Usage

Dependencies are specified in [requirements.txt](/requirements.txt)

```
pip install -r requirements.txt
```

To run the notebooks locally, you will need to have python installed,
preferably through [anaconda](https://www.anaconda.com/download/) .

You can then clone this repository. From a command line, run

```
git clone https://github.com/simpeg-research/Astic-2019-PGI.git
```

Then `cd` into the `Astic-2019-PGI` directory:

```
cd Astic-2019-PGI
```

To setup your software environment, we recommend you use the provided conda environment

```
conda env create -f environment.yml
conda activate aem-environment
```


alternatively, you can install dependencies through pypi

```
pip install -r requirements.txt
```

You can then launch Jupyter

```
jupyter notebook
```

Jupyter will then launch in your web-browser.

## Running the notebooks

Each cell of code can be run with `shift + enter` or you can run the entire notebook by selecting `cell`, `Run All` in the toolbar.

For more information on running Jupyter notebooks, see the [Jupyter Documentation](https://jupyter.readthedocs.io/en/latest/)

## Issues

Please [make an issue](https://github.com/simpeg-research/Astic_2019_PGI/issues) if you encounter any problems while trying to run the notebooks.

## Citations

If you build upon or use these examples in your work, please cite:

Astic, T. & Oldenburg, D. W. (2018). Petrophysically guided geophysical inversion using a dynamic Gaussian mixture model prior. In SEG Technical Program Expanded Abstracts 2018 (pp. 2312-2316).


```
@inbook{Astic2018,
author = {Thibaut Astic and Douglas W. Oldenburg},
title = {Petrophysically guided geophysical inversion using a dynamic Gaussian mixture model prior},
booktitle = {SEG Technical Program Expanded Abstracts 2018},
chapter = {},
pages = {2312-2316},
year = {2018},
doi = {10.1190/segam2018-2995155.1},
URL = {https://library.seg.org/doi/abs/10.1190/segam2018-2995155.1},
eprint = {https://library.seg.org/doi/pdf/10.1190/segam2018-2995155.1}
}
```

These examples were built following the structures developed by Heagy et al.(2018)

Heagy, L. J., Kang, S., Cockett, R., & Oldenburg, D. W. (2018). Open source software for simulations and inversions of airborne electromagnetic data. In 7th International Workshop on Airborne Electromagnetics (pp. 1–5).

```
@inproceedings{Heagy2018,
author = {Heagy, Lindsey J and Kang, Seogi and Cockett, Rowan and Oldenburg, Douglas W.},
booktitle = {7th International Workshop on Airborne Electromagnetics},
keywords = {finite volume,frequency domain,inversion,open source software,time domain},
pages = {1--5},
title = {{Open source software for simulations and inversions of airborne electromagnetic data}},
year = {2018}
}
```


## License
These notebooks are licensed under the [MIT License](/LICENSE) which allows academic and commercial re-use and adaptation of this work.

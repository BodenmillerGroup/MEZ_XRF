# MEZ-XRF publication repository
Multi-element Z-tag X-ray fluorescence (MEZ-XRF) is a multiplex bioimaging microscopy method. It is based on X-ray fluorescent mapping of elements targeted to different features of interest in tissue via affinity reagents (e.g., antibodies) conjugated to different (Z) elements (Z-tags). The method was developed and published in [*Multi-element Z-tag imaging by X-ray fluorescence microscopy for next-generation multiplex imaging*](https://www.nature.com/articles/s41592-023-01977-x).

This repository contains Python code and Jupyter notebooks for processing and analysis of the MEZ-XRF (and related imaging mass cytometry, IMC) raw data presented in the linked publication. This includes code to generate the panels in the manuscript Figures.

## Raw data
MEZ-XRF, IMC and microscopy raw data analysed in [*Multi-element Z-tag imaging by X-ray fluorescence microscopy for next-generation multiplex imaging*](https://www.nature.com/articles/s41592-023-01977-x) are available at [https://doi.org/10.5281/zenodo.7949102](https://doi.org/10.5281/zenodo.7949102).

XRF data were collected at the high energy beamline ID15A of the [European Synchrotron Research Facility](https://www.esrf.fr/home/UsersAndScience/Experiments/StructMaterials/ID15A.html). IMC data were collected on Hyperion instruments and microscopy data on an Olympus slidescanner in the Bodenmiller lab. See the publication for details.

## Repo layout
The repository contains the following folders:

| Setup | Contains a conda `environment.yml` for setting up a `MEZ-XRF` Python environment to run analysis |
| Code | Python code for running analysis |
| Notebooks | Notebooks for processing data |

*Disclaimer: i ~~learnt~~ to code through this project, so expect  ameteurish decisions and style :).

## Code rationale
### high_plex_hdf.py
To simplify downstream analysis, multichannel images (2D XRF data and related imaging mass cytometry data) were repackaged to a `high_plex_hdf.py` class. The `hph` class can be written to disk as a `.hdf` file.

The `hph` class stores:
- *images*, can store multiple versions including raw, normalised or filtered
- *masks*, mainly single-cell segmentation masks
- *channel metadata*, stored as a Pandas dataframe, useful to tie multiple names and observations to channels
- *anndata*, single-cell measurements directly related to the above masks and images

The `hph` class has associated functions including single-cell segmentation with [Deepcell](https://github.com/vanvalenlab/deepcell-tf/tree/master) and handy plotting functions used in multiple notebooks (e.g., to link high dimension plot cluster colour schemes to colour matched 2D single-cell segmentation mask plots).

### hph_adata_analysis.py
`hph_adata_analysis.py` handles multiple `hph.hdf` files for grouped analysis and ploting (e.g., to automatically match grayscale levels for single channel plots of different images for intensity comparison).

## Notebooks
Notebooks were used to process and analyse raw data. A brief description of each notebook follows.

### Data processing notebooks
- [1_XRF_Notebook](link)
  - Stitches serial XRF scans together (useful for interrupted scans that were restarted)
- [2_XRF_Notebook](link)
  - Deconvolutes XRF scans based on a PyMCA configuration file (supplied with raw data)
- [3_XRF_Notebook](link)
  - Normalises XRF channel intensities to beam intensity (which drops ~2 % every hour before beam top up), and repackages XRF data to `hph.h5` files

### XRF data analysis notebooks
4. [4_XRF_Notebook](link)
  - **(Figure X)**
  - Measures XRF limits of detection from a gelatin standard series
5. [5_XRF_Notebook](link)
  - **(Figure X)**
  - Segments single-cells from epithelial cell line images imaged with the silicon drift diode (SDD) detector
6. [5_XRF_Notebook](link)
  - **(Figure X)**
  - Cytometry of epithelial cell line single-cells imaged with the silicon drift diode (SDD) detector
7. [6_XRF_Notebook](link)
  - **(Figure X)**
  - Segments single-cells from breast cancer tissue microarray samples imaged with the silicon drift diode (SDD) detector
8. [7_XRF_Notebook](link)
  - **(Figure X)**
  - lorem ipsum

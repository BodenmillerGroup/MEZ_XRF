# MEZ-XRF publication repository
Multi-element Z-tag X-ray fluorescence (MEZ-XRF) is a multiplex imaging method. The method was developed and published in [*Multi-element Z-tag imaging by X-ray fluorescence microscopy for next-generation multiplex imaging*](www.paper.com).

This repository contains Python code and Jupyter notebooks for processing and analysis of the MEZ-XRF (and related imaging mass cytometry, IMC) raw data presented in the publication. This includes code to generate the panels in the manuscript Figures.

## Raw data
The MEZ-XRF and IMC raw data analysed in [*Multi-element Z-tag imaging by X-ray fluorescence microscopy for next-generation multiplex imaging*](www.paper.com) is available at [zenodo doi](www.zenodo.com).

XRF data were collected at the high energy beamline ID15A of the [European Synchrotron Research Facility](www.ESRF.com). IMC data were collected on Hyperion instruments in the Bodenmiller lab. See the publication for details.

## Repo layout
The repository contains the following folders:

| Setup | Contains a conda `environment.yml` for setting up a `MEZ-XRF` Python environment to run analysis |

| Code | Python code for running analysis |

| Notebooks | Notebooks for processing data (see Nobelow for breakdown) |

## Code rationale
### high_plex_hdf
To simplify downstream analysis, multichannel images (2D XRF data and related imaging mass cytometry data) were repackaged to a [`high_plex_hdf.py`](link) class. The `hph` class can be written to disk as a `.hdf` file.

The `hph` class stores:
- *images* (multiple versions, i.e. raw, normalised, filtered)
- *masks* (mainly single-cell segmentation masks)
- *channel metadata* (as a Pandas dataframe, useful to tie multiple names to channels)
- [*anndata*](scanpy) single-cell measurements  

`hph` class functions include single-cell segmentation with [Deepcell](link) and handy plotting functions (e.g. to link high dimension plot cluster colour schemes to colour matched 2D single-cell segmentation mask plots).

[`hph_adata_analysis.py`](link) adds a layer on top of `hph`, to handle multiple `hph` files for grouped analysis and ploting (e.g. to automatically match grayscale levels for single channel plots of different images for intensity comparison).

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

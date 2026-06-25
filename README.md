# Effective numerosity drives adaptation: map-specific recalibration of number-selective cortex revealed by perceptual grouping
The data analysis code for the paper of "Effective numerosity drives adaptation: map-specific recalibration of number-selective cortex revealed by perceptual grouping"

# Analysis code for: *Effective numerosity drives adaptation: map-specific recalibration of number-selective cortex revealed by perceptual grouping*

This repository contains the analysis and visualization code, together with anonymized processed data, for the manuscript:

**Effective numerosity drives adaptation: map-specific recalibration of number-selective cortex revealed by perceptual grouping**

## Overview

This repository provides the code used to analyze the behavioral and neuroimaging data and to generate the figures and tables reported in the manuscript. It includes scripts for the main behavioral experiment, the additional behavioral experiment, statistical analyses, and figure/table visualization for both the main text and the Supplementary Information.

To facilitate transparency and reproducibility, the repository also includes organized behavioral data and anonymized processed neuroimaging-derived data that can be shared without revealing participant identity.

## Repository Structure

```text
code/
├── data/
│   ├── T_diff_LME_withLinearAndLogDV.mat
│   └── behavioral data/
│       ├── additional behavioral experiment data/
│       ├── sub01_adaptation condition.mat
│       ├── sub01_addition.mat
│       ├── sub01_no adpatation.mat
│       ├── sub01_PSE of 40 connected.mat
│       ├── ...
├── behavioral_experiment.m
├── Figure_2A_C.m
├── Figure_3A_B.m
├── Figure_3C.m
├── Figure_4B_TableS1.m
├── Figure_S1.m
├── Table1_TableS4.m
└── remaining_Figures_and_Tables.m
```

## File Description

### Analysis scripts
* `behavioral_experiment.m`
  Main analysis script for the behavioral experiment.

* `Figure_2A_C.m`
  Code for generating Figure 2A–C.

* `Figure_3A_B.m`
  Code for generating Figure 3A-B.

* `Figure_3C.m`
  Code for generating Figure 3C.

* `Figure_4B_TableS1.m`
  Code for generating Figure 4B and Table S1.

* `Figure_S1.m`
  Code for generating Supplementary Figure S1.

* `Table1_TableS4.m`
  Code for generating Supplementary Tables 1 and S4.

* `remaining_Figures_Tables.m`
  Code for generating the remaining figures and tables not covered by the scripts listed above.

### Data files

* `data/T_diff_LME_withLinearAndLogDV.mat`
  Anonymized processed neuroimaging-derived data used for statistical modeling and figure/table generation.

* `data/behavioral data/`
  Behavioral data for the main and additional experiments.

Within the behavioral data folder, the file types have the following meanings:

* `*_no adpatation.mat`
  Behavioral data from the **no adaptation** condition.

* `*_adaptation condition.mat`
  Behavioral data collected under the different **adaptation conditions**.

* `*_PSE of 40 connected.mat`
  Data used to estimate the participant’s **subjective perceived numerosity** for the display containing **40 connected dots**.

* `*_addition.mat`
  Data from the **additional behavioral experiment**, collected in a later follow-up session.

* `additional behavioral experiment data/`
  Folder containing organized data related to the additional behavioral experiment.

## Notes on Behavioral Dataset Organization

The behavioral dataset includes data from the main experiment and the follow-up additional behavioral experiment.

### Main behavioral experiment

The main behavioral analyses reported in the manuscript were primarily based on **participants 1–13**, because at that stage the additional behavioral experiment had not yet been conducted.

### Additional behavioral experiment

The additional behavioral experiment was conducted later and includes data from **participants 1–26**.

For some participants, there is **no separate `*_addition.mat` file**. In those cases, this indicates that the participant completed the full experimental protocol (main experiment plus additional behavioral experiment), and the relevant follow-up data were retained within the corresponding `*_adaptation condition.mat` file rather than saved as a separate addition file.

## Data Availability

This repository includes:

* anonymized processed neuroimaging-derived data
* behavioral data from the formal experiment
* behavioral data from the additional behavioral experiment

These shared data have been organized so that they do **not reveal participant identity**.

The repository does **not** include raw participant-specific neuroimaging data in native individual space. Raw data are not publicly shared because of privacy considerations, informed consent limitations, and ethical restrictions associated with participant-specific imaging data.

## Reproducibility

This repository supports reproducibility of the reported analyses by providing:

* the analysis scripts used in the manuscript
* the scripts used to generate figures and tables
* the processed/anonymized data required for these analyses

Because raw neuroimaging data in individual space are not included, the repository does not provide full reproduction from raw acquisition onward. However, it does provide the computational workflow underlying the reported statistical results and visualizations.

## Software Requirements

The analyses were conducted in **MATLAB**.

Please use a recent MATLAB release and ensure that all required toolboxes are installed. Before running the scripts, add the `code/` folder and its subfolders to the MATLAB path.

If necessary, local file paths in the scripts should be adjusted to match the directory structure on your machine.

## How to Use

1. Clone or download this repository.
2. Open MATLAB.
3. Add the `code/` directory and all subfolders to the MATLAB path.
4. Run the relevant script depending on the figure, table, or analysis to be reproduced.

Typical usage:

* Run `behavioral_experiment.m` for the main behavioral analyses.
* Run `additional_behavioral_experiment.m` for the follow-up behavioral analyses.
* Run `Figure2A_C.m` for Figure 2.
* Run `Figure_3A_B.m`, and `Figure_3C.m` for Figure 3.
* Run `Figure_4B_TableS1.m` for Figure 4B and Table S1.
* Run `Figure_S1.m` for Supplementary Figure S1.
* Run `Table1_TableS4.m` for Supplementary Tables 1 and Table S4.
* Run `remaining_Figures_Tables.m` for the remaining manuscript outputs.

## Mapping Between Scripts and Manuscript Outputs

| Manuscript output                         | Script                               |
| ----------------------------------------- | ------------------------------------ |
| Main behavioral analyses                  | `behavioral_experiment.m`            |
| Additional behavioral experiment analyses | `additional_behavioral_experiment.m` |
| Figure 2A–C                               | `Figure_2A_C.m`                      |
| Figure 3A-B                               | `Figure_3A_B.m`                      |
| Figure 3C                                 | `Figure_3C.m`                        |
| Figure 4B                                 | `Figure_4B_TableS1.m`                |
| Table S1                                  | `Figure_4B_TableS1.m`                |
| Figure S1                                 | `Figure_S1.m`                        |
| Tables 1                                  | `Table1_TableS4.m`                   |
| Tables S4                                 | `Table1_TableS4.m`                   |
| Other figures/tables                      | `remaining_Figures_and_Tables.m`     |

## Citation

If you use this code or data, please cite the associated manuscript:

> [ ].

If the manuscript is still under review, you may use:

> [ ].

## License

This repository is provided for academic research use only.

If you wish to apply an open-source license, this section can be replaced with the appropriate license text.

## Contact

For questions regarding the code, data organization, or analyses, please contact:

* Yuxuan Cai
* School of Psychology, South China Normal University
* y.cai@m.scnu.edu.cn


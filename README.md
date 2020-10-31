# covid-npi
Analysis code for:

**A viral perspective on worldwide non-pharmaceutical interventions against COVID-19**

Jean-Philippe Rasigade, Anaïs Barray, Julie Teresa Shapiro, Charlène Coquisart, Yoann Vigouroux, Antonin Bal, Grégory Destras, Philippe Vanhems, Bruno Lina, Laurence Josset, Thierry Wirth.

[medRXiv preprint manuscript](https://www.medrxiv.org/content/10.1101/2020.08.24.20180927v2)

## Requirements
R software (tested on Windows 10 with v3.6.1).

Packages
* ape
* data.table
* lubridate
* survival
* coxme
* stringr
* beeswarm
* plotrix
* wordcloud

A workstation with at least 24 CPU cores is recommended to run scripts involving mixed-effect models. Use the special flag **rebuild = FALSE** to bypass these computations.

## Usage
Scripts should be run individually in the following order:

### Pre-treatment scripts

These scripts perform alignment operations, especially between country identifiers; and extraction of branch length data.

* S1_dataprep_countries.R
* S2_dataprep_tree.R
* S3_dataprep_branch_lengths.R
* S4_dataprep_npis.R

### Analysis scripts

These scripts generate the results and figures of the study.

* A1_descriptive_stats_by_country.R
* A2_branch_length_distribution.R
* A3_divergence_rates_per_country.R
* A4_npi_timing.R
* A5_npi_univariate.R
* A6_npi_multivariate.R
* A7_SIR_models.R

### Flowchart of selection steps
* Z1_flowchart.R
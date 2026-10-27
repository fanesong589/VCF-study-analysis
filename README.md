Analysis Code for Vertebral Compression Fracture Reporting and Subsequent Outcomes

This repository contains the statistical analysis code for a multicenter retrospective cohort study of vertebral compression fracture reporting on abdominopelvic CT and subsequent fractures and clinical care. Patient-level study data and the commercial AI screening software are not included.

Repository Contents

README.md: Study overview and instructions.

analysis.R: Statistical analysis code.

data_dictionary.csv: Required variables, definitions, units, and coding.

All three files are located in the repository root.

Analysis Overview

The code covers descriptive statistics, primary and sensitivity analyses, and pooled adjusted probabilities and risk differences.

Primary analyses use Firth’s penalized-likelihood logistic regression for two outcomes: an additional moderate or severe vertebral compression fracture and subsequent clinical care within 90 days. Models are fitted for the pooled cohort and separately for each center.

Primary models include age, sex, body mass index, moderate and severe fracture counts, and prior fracture-related procedure status. Pooled models additionally include center. Full-adjustment sensitivity models additionally include hypertension, diabetes mellitus, and chronic kidney disease. A second sensitivity analysis excludes patients with prior fracture-related procedures and omits procedure status from the pooled primary model.

Adjusted odds ratios, profile penalized-likelihood 95% confidence intervals, and raw two-sided penalized likelihood-ratio P values are calculated. No multiple-comparison adjustment is applied.

For the pooled cohort, adjusted probabilities and risk differences are estimated by marginal standardization of conventional maximum-likelihood logistic regression models with the primary covariates. Reporting status is set to reported and unreported for each patient, and predicted probabilities are averaged across the analytic sample. Risk differences are calculated as reported minus unreported.

Percentile 95% confidence intervals are obtained from 2,000 patient-level bootstrap samples per outcome, stratified by center and reporting status. Each bootstrap sample undergoes model refitting and standardization. The random seed is 42.

Software Requirements

The analysis was validated using:

R 4.6.1

logistf 1.26.1

Install the validated logistf version and its dependencies as follows:

install.packages("remotes")
remotes::install_version(
  "logistf",
  version = "1.26.1",
  repos = "https://cloud.r-project.org"
)

Input Data

Prepare a CSV file with one row per patient and the variables specified in data_dictionary.csv. Keep patient-level data outside the public repository.

Center 1 represents Gangnam and Center 2 represents Sinchon. Center 2 is the reference category. Reporting status is coded as 0 for unreported and 1 for reported. Original procedure codes 1 and 3 indicate prior fracture-related procedures, whereas codes 0 and 2 indicate their absence.

The clinical care composite is positive when either a related clinic visit or initiation of a qualifying medication is recorded within 90 days. Patients meeting both criteria are counted once. Qualifying medications follow Supplementary Table S1 of the manuscript, including calcium and vitamin D preparations.

The code starts from the curated analytic cohort after eligibility assessment, image review, and electronic medical record extraction. It does not perform these clinical assessments. All 286 study patients had complete data for the reported models, and no imputation was performed.

Running the Analysis

Open analysis.R in R or RStudio and follow its input and output instructions. Use an authorized analytic dataset that conforms to data_dictionary.csv.

To reproduce the original bootstrap estimates precisely, preserve the original patient row order, factor levels, random-number settings, and software versions.

Interpretation and Data Availability

Adjusted probabilities describe the selected analytic sample. They should not be interpreted as population incidence estimates, exact five-year cumulative incidence, or causal effects.

Patient-level data are not publicly distributed in this repository. Access requests should follow the manuscript’s data sharing statement and institutional requirements. Reproduction of the study estimates requires authorized access to the underlying study data.

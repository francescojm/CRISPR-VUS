# CRISPR-VUS

## Overview

This repository contains the analysis workflow for the CRISPR-VUS project.

The repository includes a Google Colab notebook and a set of R files. The Colab notebook is the main file used to run the analysis. The R files contain functions that are sourced and used directly inside the notebook.

---

## Data retrieval

This section will describe how to retrieve the data required to run the analysis.

To be completed.

Information to add later:

- Data sources.
- Download links or access instructions.
- Expected folder structure.

---

## Notebook organization

The notebook is divided into 19 sections.

Sections 1 to 4 are required setup sections. They must be run every time before running any later section of the notebook.

Sections 5 to 19 are independent from each other after sections 1 to 4 have been run, provided that the required intermediate results already exist.

For example, to run section 15, it should be sufficient to run:

```text
1 → 2 → 3 → 4 → 15
```

However, during the first complete run of the notebook, it is mandatory to run sections 1 to 8 in order. Sections 1 to 8 compute and aggregate the main results required by the rest of the workflow.

Therefore:

- For any notebook run, sections 1 to 4 must always be executed first.
- For the first full run, sections 1 to 8 must be executed.
- After the first full run, downstream sections can be run independently, as long as sections 1 to 4 are run first and the required intermediate results are available.

---

## Notebook sections

| Section | Title | Notes |
|---|---|---|
| 1 | User settings and analysis parameters | Required for every run |
| 2 | Data loading | Required for every run |
| 3 | Preprocessing | Required for every run |
| 4 | Tissue selection | Required for every run |
| 5 | DAM tissue-specific analysis | Required for the first complete run |
| 6 | Integration with IntOGen driver information | Required for the first complete run |
| 7 | Drug response analysis | Required for the first complete run |
| 8 | Post-processing, first summaries, and driver enrichment analyses | Required for the first complete run |
| 9 | Global summaries and visualization | Independent after sections 1 to 4 and required intermediate results |
| 10 | Reactome enrichment | Independent after sections 1 to 4 and required intermediate results |
| 11 | Cross-tissue DAM analysis | Independent after sections 1 to 4 and required intermediate results |
| 12 | STRING analysis | Independent after sections 1 to 4 and required intermediate results |
| 13 | Co-occurrence analysis | Independent after sections 1 to 4 and required intermediate results |
| 14 | DR validation | Independent after sections 1 to 4 and required intermediate results |
| 15 | DAM position annotation | Independent after sections 1 to 4 and required intermediate results |
| 16 | Add SIFT and PolyPhen scores | Independent after sections 1 to 4 and required intermediate results |
| 17 | IntOGen and COSMIC patient analysis | Independent after sections 1 to 4 and required intermediate results |
| 18 | Summary of patient data | Independent after sections 1 to 4 and required intermediate results |
| 19 | DAM patient actionability and tractability analysis | Independent after sections 1 to 4 and required intermediate results |

---

## How to run the notebook

### First complete run

At the first complete run, execute sections 1 to 8 in order:

```text
1 → 2 → 3 → 4 → 5 → 6 → 7 → 8
```

These sections perform the initial setup and compute and aggregate the main results.

### Running a specific section after the first complete run

After the first complete run has been completed, sections 9 to 19 can be run independently.

Before running any of these sections, always run sections 1 to 4.

Example: to run section 15, execute:

```text
1 → 2 → 3 → 4 → 15
```

This assumes that the main results computed and aggregated by sections 1 to 8 are already available.

---

## Dependencies

To be completed.

This section will include the software and package requirements needed to run the notebook and the R scripts.

Information to add later:

- R version.
- Required R packages.
- Colab-specific setup instructions.



## Citation

To be added before publication.
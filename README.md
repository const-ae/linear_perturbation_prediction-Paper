# Code repository for *Deep learning-based predictions of gene perturbation effects do not yet outperform simple linear methods*

This repository contains the code to reproduce our analysis.

- The **notebooks** folder contains the R scripts used for the analysis and to make the figures
  - [Double Perturbation Analysis](https://htmlpreview.github.io/?https://github.com/const-ae/linear_perturbation_prediction-Paper/blob/main/notebooks/double_perturbation_analysis.html)
  - [Single Perturbation Analysis](https://htmlpreview.github.io/?https://github.com/const-ae/linear_perturbation_prediction-Paper/blob/main/notebooks/single_perturbation_analysis.html)
- The **benchmark** folder contains the scripts to reproduce the benchmark results
  - The **benchmark/src** contains individual scripts to run each method
  - The **benchmark/conda_environments** and **benchmark/renv** contain the details about the software versions
  - The **benchmark/submission** contains the script to launch the scripts using my [custom](https://github.com/const-ae/MyWorkflowManager) workflow manager

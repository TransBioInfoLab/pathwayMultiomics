### pathwayMultiomics: An R package for efficient integrative analysis of multi-omics datasets with matched or un-matched samples

Gabriel J. Odom, Antonio Colaprico, Tiago C. Silva, X. Steven Chen, Lily Wang       
2021-08-27  

#### For user manuals, please see <https://transbioinfolab.github.io/pathwayMultiomics/articles/intro_to_minimax.html>

Recent advances in technology have made multi-omics datasets increasingly available to researchers. To leverage the wealth of information in multi-omics data, a number of integrative analysis strategies have been proposed recently. However, effectively extracting biological insights from these large, complex datasets remains challenging. In particular, matched samples with multiple types of omics data measured on each sample are often required for multi-omics analysis tools, which can significantly reduce the sample size. Another challenge is that analysis techniques such as dimension reductions, which extract association signals in high dimensional datasets by estimating a few variables that explain most of the variations in the samples, are typically applied to whole-genome data, which can be computationally demanding. 

Here we present pathwayMultiomics, a pathway-based approach for integrative analysis of multi-omics data with categorical, continuous, or survival outcome variables. The input of pathwayMultiomics is pathway P-values for individual omics data types, which are then integrated using a novel statistic, the MiniMax statistic, to prioritize pathways dysregulated in multiple types of omics datasets. Importantly, pathwayMultiomics is computationally efficient and does not require matched samples in multi-omics data. We performed a comprehensive simulation study to show that pathwayMultiomics significantly outperformed currently available multi-omics tools with improved power and well-controlled false-positive rates. In addition, we also analyzed real multi-omics datasets to show that pathwayMultiomics was able to recover known biology by nominating biologically meaningful pathways in complex diseases such as Alzheimer???s disease.

To install this package, you will need the `devtools::` package and either [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Windows) or [Xcode](https://developer.apple.com/xcode/) (Mac) to build packages from GitHub. Then, install this package via

```
library(devtools)
install_github("TransBioInfoLab/pathwayMultiomics")
```

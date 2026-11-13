options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("ggplot2", "vegan", "ggpubr", "dplyr", "ggrepel", "ggpp", "factoextra", "knitr", "GOplot", "ape", "stringr", "remotes", "BiocManager")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq", ask = FALSE, update = FALSE)
if (!requireNamespace("fantaxtic", quietly = TRUE)) remotes::install_github("gmteunisse/fantaxtic", upgrade = "never")

# Lobaria pulmonaria microbiome across a European climate gradient

R scripts and processed datasets accompanying **Bogale et al. (2024), “The microbiome of the lichen Lobaria pulmonaria varies according to climate on a subcontinental scale”**, published in *Environmental Microbiology Reports*, 16(3), e13289.

[Read the article](https://doi.org/10.1111/1758-2229.13289) · [Raw proteomics data: PXD049079](https://www.ebi.ac.uk/pride/archive/projects/PXD049079)

## Study

The study investigates bacterial microbiome composition and protein function in 42 *L. pulmonaria* thalli from 14 sites across five European countries. Analyses connect metaproteome variation with long-term climate, short-term weather, sampling season, and fungal and algal genotypes.

Climatic regions: Sub-Atlantic Lowland (SAL), Alpine (ALP), Hercynian Montane (HM), and Subcontinental Lowland (SCL).

## Contents

| Path | Contents |
| --- | --- |
| `analysis/lichen_data_analysis_original.R` | Original supplied script, preserved byte for byte |
| `analysis/lichen_data_analysis.R` | Portable copy with minimal execution preparation |
| `analysis/run_analysis.R` | Runner that resolves paths and saves plot pages and session information |
| `data/` | Eleven supplied CSV files, preserved byte for byte |
| `results/volcanoPlot_original.svg` | Original supplied plot |
| `docs/DATA_DICTIONARY.md` | Data descriptions and analysis coverage |
| `docs/PROVENANCE.md` | Source and exact preparation changes |
| `CITATION.cff` | Citation metadata for the associated paper |

## Run

Install R and a compatible set of packages, then open a terminal in the repository root:

```sh
Rscript scripts/install_dependencies.R
Rscript analysis/run_analysis.R
```

The runner writes `results/analysis_plots.pdf` and `results/sessionInfo.txt`. Package installation requires internet access. Historical package versions were not supplied; the installer does not reconstruct the original computational environment. For interactive work, set the working directory to `data/` and source `analysis/lichen_data_analysis.R` using its full path.

## Analysis coverage

The supplied script covers climate clustering, taxonomic composition, NMDS with environmental fitting, visualization of saved PERMANOVA results, plots of saved limma results, GO annotations, and genotype ordinations. PERMANOVA and limma tables are provided as precomputed results; the original fitting code is not included. Raw spectrum processing, imputation benchmarking, genotype clustering/Mantel tests, the sampling map, and the eggNOG functional-composition panel are outside this script's coverage. This is a companion code and data archive, not a complete raw-data-to-manuscript workflow.

## Validation status

Dataset dimensions, sample identifiers, numerical values, taxonomic table alignment, significance counts, and original-file checksums were checked in the preparation environment. R was unavailable there: package installation and full figure generation have **not** been executed. The included GitHub Actions workflow checks R syntax and data integrity; it does not run the full analysis. No historical `renv.lock` or package versions are claimed.

## Citation and attribution

Please cite the article above when using these study materials. Full author metadata are in `CITATION.cff`.

The article is distributed under Creative Commons Attribution, as stated in the publication. No separate software or dataset license was supplied in the source archive; this repository does not assign a new license to those materials. See `docs/PROVENANCE.md`.

Maintainer: [Anteneh Tamirat Bogale](https://github.com/AntenehTamiratBogale).

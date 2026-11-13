# Processed datasets

CSV files are read using the source script’s encodings; retain original bytes and accented site names. Counts below exclude headers and include identifier columns.

| File | Rows | Columns | Description |
| --- | ---: | ---: | --- |
| `DEPs_limma_output.csv` | 73 | 10 | 73 selected differential protein records from the limma output. |
| `GO_Terms_level_1.csv` | 292 | 9 | Four GO category blocks for the 73 differential proteins. |
| `algae_genotypes.csv` | 46 | 16 | Algal microsatellite fragment lengths; country and site plus 14 loci. |
| `all_protein_limma_output.csv` | 217 | 11 | Precomputed limma results for 217 proteins; includes labels and activity classes. |
| `data_env_prot.csv` | 42 | 237 | 42 samples, 20 metadata/environment columns, followed by 217 protein abundance columns. |
| `data_permanova_output.csv` | 12 | 3 | Precomputed PERMANOVA R-squared components, including residuals. |
| `fungi_genotypes.csv` | 46 | 10 | Fungal microsatellite fragment lengths; country and site plus eight loci. |
| `meta_data.csv` | 14 | 2 | Climatic region assignment for 14 sampling sites. |
| `oxidoreductase_GOterms.csv` | 31 | 6 | Binary oxidoreductase annotations and logFC for 31 proteins. |
| `tax_abundance_data.csv` | 58 | 15 | Taxon abundance table: 58 rows by 14 sites, plus row identifier. |
| `tax_class_data.csv` | 58 | 3 | Phylum and family for the same 58 abundance rows, plus row identifier. |

## Environmental fields

`sites_names` identifies individual thalli; `sites` identifies locations. `SSJ` is sampling season as Julian day; `T2W`, `H2W`, and `P2W` represent two-week temperature, humidity, and precipitation cover. `AMT`, `MDR`, `TS`, `AP`, `PDQ`, and `PS` are the six bioclimatic variables. Interpret units and source scaling using the publication and original environmental data sources. Columns 21–237 contain the processed protein abundance matrix used by the script.

## Interpretation

The limma table records the contrast string `GroupSAL- GroupALP`; retain the source direction when interpreting logFC. `B_statistics` and `neg_log10_p_val` are separate supplied fields. The taxonomic abundance table contains metaproteomic abundances, not amplicon ASV counts. The random tree constructed by the script is a plotting data structure and is not an inferred phylogeny. Genotype tables contain 46 records and must not be matched by row position to the 42 metaproteome samples. Saved result tables are inputs, not newly fitted results.

Protein display names vary slightly between the two limma tables; differential records are checked using their numerical logFC and adjusted P values.

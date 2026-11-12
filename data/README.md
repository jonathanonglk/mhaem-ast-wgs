# MIC phenotype files

Ordinal MIC rank encodings used as pyseer input by `env1_pipeline.sh`
(SNP GWAS) and `env2_pangenome.sh` (pan-genome GWAS).

Format: tab-separated, header `Samples<TAB>Trait`, one row per isolate.
60 rows: 56 clinical isolates from unique patients, plus replicates.

Values are **ordinal ranks, not µg/mL**. Rank increases with MIC, one step
per doubling dilution of the SLOMYCO panel. CLSI M24S non-susceptibility
breakpoints fall at:

| Drug | Non-susceptible | Rank |
|---|---|---|
| CIP | ≥ 2 µg/mL  | ≥ 5 |
| AMI | ≥ 32 µg/mL | ≥ 6 |
| RIF | ≥ 2 µg/mL  | ≥ 5 |

For CIP, ranks 4–7 correspond to 1, 2, 4 and 8 µg/mL. Full rank-to-MIC
mappings are in Supplementary Tables 1 and 3 of the manuscript.

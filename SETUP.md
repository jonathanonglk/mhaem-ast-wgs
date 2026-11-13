# Environment setup

Two conda environments are needed. They are kept separate because the
pan-genome tools and the variant-calling tools pin incompatible versions of
their shared dependencies.

| | `env_1` | `env_2` |
|---|---|---|
| Script | `env1_pipeline.sh` | `env2_pangenome.sh` |
| Does | trimming, QC, Snippy variant calling, snippy-core, IQ-TREE, SNP GWAS, targeted resistance genes, intein status | assembly, annotation, pan-genome, pan-genome GWAS, Figure 4, verification checks |
| Run order | first | second (uses `env_1` outputs) |

Install [Miniforge](https://github.com/conda-forge/miniforge) (or Miniconda)
first, then create the environments below.

---

## Option 1 (preferred): recreate from the exported environment files

These reproduce the exact package versions used for the published analysis.

```bash
conda env create -n env_1 -f env_mhaem_wgs.yml
conda env create -n env_2 -f env_pangenome.yml
```

---

## Option 2: create from package lists

Use this if the YAML files are unavailable. Solver results can differ, so
check the versions afterwards.

### env_1

```bash
conda create -y -n env_1 -c conda-forge -c bioconda python=3.10
conda install -y -n env_1 -c conda-forge -c bioconda \
    snippy=4.6.0 fastp fastqc multiqc iqtree \
    pyseer=1.4.0 samtools=1.23 bcftools=1.23.1 \
    dendropy pandas numpy scipy matplotlib
```

Snippy supplies BWA, FreeBayes, SnpEff and SAMtools. The versions used for the
published analysis were BWA-MEM 0.7.18, FreeBayes 1.3.6, SnpEff 5.0e,
SAMtools 1.23, BCFtools 1.23.1, IQ-TREE 3.1.3, pyseer 1.4.0, fastp 1.3.3,
FastQC 0.12.1 and MultiQC 1.35.

### env_2

```bash
conda create -y -n env_2 -c conda-forge -c bioconda python=3.8
conda install -y -n env_2 -c conda-forge -c bioconda \
    spades=4.3.0 quast=5.3.0 prokka=1.15.6 panaroo=1.8.0 \
    iqtree=3.1.3 pyseer=1.4.2 samtools bcftools \
    biopython=1.83 pandas scipy numpy matplotlib openpyxl
```

The versions used for the published analysis were SPAdes 4.3.0, QUAST 5.3.0,
Prokka 1.15.6, Panaroo 1.8.0, IQ-TREE 3.1.3, pyseer 1.4.2, Biopython 1.83,
pandas 2.0.3, scipy 1.7.0 and matplotlib 3.7.3.

---

## Check the environments

```bash
conda activate env_1
snippy --version && bwa 2>&1 | grep -i version && freebayes --version
samtools --version | head -1 && bcftools --version | head -1
iqtree --version | head -1 && pyseer --version
fastp --version && fastqc --version && multiqc --version
python -c "import pandas, numpy, scipy, matplotlib, dendropy; print('python deps OK')"

conda activate env_2
spades.py --version && quast.py --version | head -1 && prokka --version
panaroo --version && iqtree --version | head -1 && pyseer --version
python -c "import Bio, pandas, scipy, matplotlib, openpyxl; print('python deps OK')"
```

---

## Input files

Both scripts take their inputs as `KEY=VALUE` arguments (or environment
variables), so nothing needs to be edited inside the scripts.

| Input | Setting | Format |
|---|---|---|
| Raw paired FASTQ files | `RAW_READS` (env_1) | one directory containing `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`; Illumina-style `_R1_001.fastq.gz` names are renamed automatically |
| Reference genome | `REF_SRC` (env_1), `REF_GBFF` (env_2) | annotated GenBank file, e.g. `GCF_000340435.2.gbff` |
| MIC rank files | `MIC_AMI`, `MIC_CIP`, `MIC_RIF` (both) | tab-separated `sample<TAB>ordinal MIC rank`, with a header row; replicates may be present and are filtered out |
| env_1 outputs | `MHAEM_WORKDIR` (env_2) | the `WORKDIR` used for env_1 |

Run `bash env1_pipeline.sh --help` or `bash env2_pangenome.sh --help` for the
full list of settings and their defaults.

---

## Running

```bash
# 1. reference-based pipeline and SNP GWAS
conda activate env_1
bash env1_pipeline.sh \
    WORKDIR=~/mhaem_wgs \
    RAW_READS=~/data/fastq \
    REF_SRC=~/data/GCF_000340435.2.gbff \
    MIC_AMI=~/data/AMI.txt \
    MIC_CIP=~/data/CIP.txt \
    MIC_RIF=~/data/RIF.txt \
    THREADS=8

# 2. pan-genome analysis, Figure 4 and verification checks
conda activate env_2
bash env2_pangenome.sh \
    MHAEM_WORKDIR=~/mhaem_wgs \
    OUTDIR=~/pangenome \
    MIC_AMI=~/data/AMI.txt \
    MIC_CIP=~/data/CIP.txt \
    MIC_RIF=~/data/RIF.txt \
    THREADS_PER_JOB=8
```

Both scripts skip any step whose output already exists, so they can be
re-run safely after an interruption. To force a step to re-run, delete its
output first.

To run only some of the env_2 verification checks:

```bash
bash env2_pangenome.sh ADHOC_SECTIONS="F G"
```

---

## Verification scripts

These are independent of the pipeline: they re-derive published numbers from
the pipeline's own output files, so a successful run should reproduce them.
Run them after `env2_pangenome.sh` has finished. They need only `python3`
(scipy is used if present, otherwise an exact fallback is used instead).

| Script | Checks |
|---|---|
| `check_qc.sh` | fastp / FastQC / MultiQC outputs: files present, read-pair integrity, per-sample and summary QC metrics |
| `check_gwas_table5.sh` | every column of **Table 5** (features tested, unique patterns, Bonferroni threshold, lowest p, nominal hits, significant hits, h², lambda) for all six drug x approach combinations, and diffs them against the published values |
| `check_cosegregation.sh` | whether the top SNP signals share their presence/absence pattern with variants elsewhere in the genome (lineage marker) or are confined to one locus |
| `compare_abc_snp_mic.sh` | MIC distributions for isolates carrying the reference vs alternative allele at the B586_RS18800 variants |
| `test_kinship_absorption.sh` | for clusters that track the gyrA intein lineage, compares pyseer's uncorrected `filter-pvalue` with the kinship-corrected `lrt-pvalue` |

```bash
conda activate env_1          # any environment with python3 will do

bash check_qc.sh WORKDIR=~/mhaem_wgs

bash check_gwas_table5.sh WORKDIR=~/mhaem_wgs PANOUT=~/pangenome

bash check_cosegregation.sh WORKDIR=~/mhaem_wgs DRUG=CIP

bash compare_abc_snp_mic.sh \
    MIC=~/data/CIP.txt \
    VCF=~/mhaem_wgs/core/core.vcf \
    SNIPPY=~/mhaem_wgs/snippy

bash test_kinship_absorption.sh WORKDIR=~/mhaem_wgs PANDIR=~/pangenome/pyseer
```

`check_gwas_table5.sh` prints `OK` or `MISMATCH` per row. All six rows should
read `OK`:

| Approach | Drug | Features | Patterns | Bonferroni | Lowest p | Nominal | Sig | h2 | lambda |
|---|---|---|---|---|---|---|---|---|---|
| SNP (Snippy) | AMI | 24,071 | 821 | 6.09e-05 | 1.10e-03 | 166 | 0 | 0.38 | 1.82 |
| SNP (Snippy) | CIP | 24,294 | 826 | 6.05e-05 | 8.24e-03 | 254 | 0 | 0.44 | 3.14 |
| SNP (Snippy) | RIF | 24,294 | 826 | 6.05e-05 | 1.39e-03 | 74 | 0 | 0.50 | 1.28 |
| Pan-genome | AMI | 1,539 | 328 | 1.52e-04 | 4.57e-03 | 7 | 0 | 0.35 | 2.91 |
| Pan-genome | CIP | 1,562 | 334 | 1.50e-04 | 7.29e-03 | 21 | 0 | 0.61 | 2.28 |
| Pan-genome | RIF | 1,562 | 334 | 1.50e-04 | 1.14e-02 | 22 | 0 | 0.65 | 1.83 |

h2 is read from the pyseer run log, not the results table. If the log was not
kept, that column reports `NA` and is skipped; rerun pyseer with
`2>&1 | tee <DRUG>.log` to capture it.

Other published values these scripts should reproduce:

- 61 libraries: 56 clinical isolates from unique patients, 4 repeat isolates and ATCC 29548; replicates and the QC strain are excluded from the clinical summaries
- gyrA intein: 420 aa, codons 131-550, reference positions 7,656-8,915 (1,260 bp), present in 14 clinical isolates and the reference
- CIP non-susceptible 49/56, RIF non-susceptible 25/56, AMI non-susceptible 40/56
- The three B586_RS18800 variants (4,021,393 / 4,021,694 / 4,021,771): 16 isolates carry the reference allele, all CIP non-susceptible; those 16 are the 14 intein-positive isolates plus CT-WN16CD4677 and CT-WN17CD1188

---

## Phenotype file encoding

The `MIC_*` files hold **ordinal MIC ranks**, not concentrations in ug/mL.
Rank increases with MIC. For ciprofloxacin, rank 5 corresponds to 2 ug/mL,
which is the CLSI M24S non-susceptibility breakpoint; ranks 4, 6 and 7 are
1, 4 and 8 ug/mL. For amikacin, rank 6 corresponds to 32 ug/mL.

`compare_abc_snp_mic.sh` therefore applies its non-susceptibility cut on the
rank scale (`NS_CUT=5` for CIP, `6` for AMI by default). Passing a cut in
ug/mL will classify every isolate as non-susceptible.

---

## Notes

- `env1_pipeline.sh` patches the installed Snippy script to remove its
  `vt normalize` step, which fails under WSL. The patch is applied to the copy
  of Snippy inside `env_1` only.
- `env2_pangenome.sh` requires `env_1` to have been run: it reads the trimmed
  reads, the Snippy output, the core alignment and the intein status table.
- Part C of `env2_pangenome.sh` runs with `set +e`, so a check that cannot find
  its input reports and the script continues.
- Resources: SPAdes is the limiting step. Allow roughly 16 GB RAM per parallel
  job and several hours for 60+ isolates on 4–8 threads.

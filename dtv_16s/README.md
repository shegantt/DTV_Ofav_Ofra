# DTV_16S

16S rRNA (V4) bacterial community analyses for *Orbicella faveolata* and *O. franksi* under stable vs. diel thermal variation (DTV) priming and heat challenge. Raw reads are on NCBI SRA (PRJNA1442122).

Run the folders in order: `01` builds the cleaned ASV table, `02` analyzes alpha diversity, and `03` analyzes community composition and makes the bacterial figures. Both `02` and `03` start from `01.make.asv.table/ps.ready.rds`, and `03` also reads the diversity table written by `02`.

Note: Clonemate genets KA (= KG) and VB (= VD) are removed before all analyses in `02` and `03`, leaving 31 samples.

## 01.make.asv.table — ASV inference and cleaning

Files are listed roughly in the order they are made/used in the script.

| File | Description |
|---|---|
| `dtv_16s_revised.Rmd` / `.html` | Script: primer check, DADA2 ASV inference, SILVA taxonomy, contaminant and low-abundance filtering, and a final chunk that prints the read and ASV numbers reported in the manuscript. |
| `dtv.16s.seqtab.rds` | DADA2 sequence table right after merging paired reads, before length filtering and chimera removal (41 samples × 2,939 ASVs; 39 coral samples + 2 negative controls). |
| `dtv.16s.readstats.csv` | Reads per coral sample at each DADA2 step (input, filtered, denoised, merged, tabled, non-chimeric). The last two rows are the mean and SD of input reads (50,073 ± 30,016). |
| `dtv.16s.seqtab.nochim.RDS` | Sequence table after keeping 244–264 bp ASVs and removing chimeras (41 samples × 2,869 ASVs); columns are full sequences. |
| `dtv_16s_seqtab.nochim.fasta` | FASTA of the 2,869 non-chimeric ASV sequences, named `sq1`, `sq2`, … |
| `dtv.16s.taxa.withids.rds` | SILVA v138.1 taxonomy for each ASV, with rows renamed to `sq` IDs and the raw sequence kept as the last column. |
| `dtv.16s.seqtab.nochim.withids.rds` | Same as `dtv.16s.seqtab.nochim.RDS`, with columns renamed to `sq` IDs. |
| `dtv_16s_sample_info_fixed.csv` | Sample metadata (host species, treatment, timepoint, genotype, fragment), with clonemate genotypes merged in the `genotype` column (`KA_KG`, `VB_VD`) and the original genet in `genotype_original`. |
| `dtv.ps.cleaner.rds` | Phyloseq object after removing chloroplasts (237), mitochondria (22), non-bacterial ASVs (31), decontam contaminants (5) and negative controls (2,574 ASVs, 39 samples). |
| `ps.ready.rds` | **Final phyloseq object** used downstream (39 samples). It drops ASVs found in only one sample or with < 49 total reads, zeroes per-sample counts < 10 (465 ASVs left; 12,689 ± 10,464 reads per sample) and fills taxonomy gaps from the nearest assigned rank. |

## 02.diversity — alpha diversity

| File | Description |
|---|---|
| `dtv.16s.div_revised.Rmd` / `.html` | Script: Shannon and inverse Simpson diversity on non-rarefied and rarefied data, with linear mixed models (`lmer`, priming × heat, random genotype) and pairwise `emmeans` for each host species. |
| `dtv.16s.alphadiv.csv` | Per-sample Shannon, inverse Simpson, observed ASVs and evenness (non-rarefied, 31 samples) with sample metadata. Written by this script and read by `03` for Fig. 3D. |
| `dtv.ps.rare.3100.rds` | `ps.ready.rds` (clonemates removed) rarefied to 3,149 reads per sample (30 samples; VE18-1 dropped with 1,881 reads), used to check that the diversity results aren't driven by read depth. Saved once, because rarefying is random. |

## 03.comm.comp — community composition

| File | Description |
|---|---|
| `dtv.16s.community_revised.Rmd` / `.html` | Script: PCoA of Bray-Curtis dissimilarity on relative abundances, negative binomial GLMMs (`glmmTMB`) on ASV counts with likelihood ratio tests (the stats reported in the manuscript), phylum/family summaries and stacked bar plots. |
| `fig3CD_bact.pdf` / `.png` | Two-panel figure: PCoA of Bray-Curtis dissimilarity and mean Shannon diversity (manuscript Fig. 3C–D, before p-values and final layout were added). |
| `gg.ms.bac.revised.pdf` / `.png` | Relative abundance of the top phyla and families for each sample (manuscript Fig. S8, before final styling). |
| `ready.otutab.csv` | ASV count table exported from `ps.ready.rds` (39 samples × 465 ASVs). |

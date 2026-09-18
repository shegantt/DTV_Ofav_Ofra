This code is associated with the statistical analysis of O.faveolata and O. franksi that were maintained in lab systems for one year before being acclimated to 2.4°C DTV oscillations and then exposed to a heat challenge experiment. 

Physiology in the form of Fv/Fm, red channel intensity, and buoyant weights were measured every 30 days throughout the experiment, concluding in 4 timepoints.

Additional physiological and gene expression measures were only sampled after thermal treatment priming and at the end of the thermal experiment. These included measures of photosynthesis and dark respiration rates, DNA for ITS2 and 16S characterization, and RNA for transcriptomic analysis.

The code and files here encompass the following analyses and files.

Metadata file: THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv

Environmental Characterization: -------------------------------------------------------------------------
  Tank water quality measures over the experiment                     -> Tank_WQ-Sheet1.csv
  Raw temperature data for control and DTV tanks over the experiment  -> hobo_all.csv
  R code for analyzing temperature data                               -> temp_analysis.Rmd

Physiological analysis: ----------------------------------------------------------------------------------
  R code for analysis of all physiological measures across timepoints -> DTV_Phys_Expt.R

Gene Expression Sequence Processing: ---------------------------------------------------------------------
  Code for high-performance computing processing of RNA sequences    -> GE_MappingToCounts_README.txt

Gene Expression Analysis: --------------------------------------------------------------------------------

  Host files:
  Code for calling single-nucleotide polymorphisms for clone ID      -> DTV_SNPs_SymID_README.txt
  Isoform to gene translation file for sequence identification       -> Ofav_Made_iso2geneName_16Dec2024.tab
  Isoform to GO term identification for sequence identification      -> Ofav_Made_iso2go_nrify_16Dec2024.tab
  R code for Coral Host RNA sequence analysis                        -> Host_GE_DeSeqPCA&Venn.R
  R code for Coral Host WGCNA analysis                               -> Host_GE_WGCNA.R

  Sybiodiniaceae files:
  R code for characterizing Symbiodiniaceae communities from RNA     -> RNA-based_SymProportions.R
  R code for Symbiodiniaceae RNA sequence analysis                   -> Sym_GE_PCA&Processing&Plasticity&Venn.R

ITS2 - Symbiodiniaceae Community Characterization:
  Raw Symbiodiniaceae ITS2 reads                                      -> dtv_its2_rawreads.xslx
  R code for Symbiodiniaceae analysis                                 -> dtv_its2_revised.Rmd
  Metadata for Symbiodiniaceae samples                                -> dtv_its2_sampleinfo_fixed.csv
  Symportal output - Taxonomy                                         -> dtv_symportal_taxa.csv
  Symportal output ITS2 type profiles                                 -> dtv_symportal_type_profiles.csv

16S Bacterial Community Characterization:
  R code for 16S analysis                                             -> 
  16S taxonomy identification                                         ->
  16S ASV counts                                                      ->
  16S ASV Metadata                                                    ->




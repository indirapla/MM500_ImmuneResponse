#https://bioconductor.org/packages/release/bioc/vignettes/easier/inst/doc/easier_user_manual.html#1_Introduction

#BiocManager::install("easier")

library("easier")
library(AnnotationDbi)
library(org.Hs.eg.db)
library('tidyverse')
#============================================================================

# protein expression data (1)
RNA_tpm_seg <- as.data.frame(readxl::read_excel("data/Cohort2_rna.xlsx")) # Open the data
rownames(RNA_tpm_seg) <- RNA_tpm_seg$Gene
RNA_tpm_seg <- RNA_tpm_seg[,2:ncol(RNA_tpm_seg)]

# gene_protID = clusterProfiler::bitr(prot.exp1_protID, fromType="UNIPROT", 
#                                            toType="ALIAS", 
#                                            OrgDb="hgu95av2.db") #  convert biological IDs to geneSymbol gene ID.
# 
# prot.exp1$UNIPROT <- rownames(prot.exp1)
# 
# prot.exp1a <- plyr::join_all(list(prot.exp1,gene_protID), by='UNIPROT',match='first')
# 
RNA_tpm_seg1 <- na.omit(RNA_tpm_seg)


# RNA expression data (2)
RNA_counts_MM <- #assays(dataset_mariathasan)[["counts"]]

# cohort cancer type 
cancer_type_MM <- 'SKCM'   # skin cutaneous melanoma (SKCM)

#================ANALISIS=======================================================
# 3.2 Compute hallmarks of immune response
hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", 
                                  "IFNy", "Ayers_expIS", "Tcell_inflamed",
                                  "RIR", "TLS")

immune_response_scores_MM <- compute_scores_immune_response(RNA_tpm = RNA_tpm_seg1, 
                                                         selected_scores = hallmarks_of_immune_response)
head(immune_response_scores_MM)
#write.csv2(immune_response_scores_MM,'output/RNA_2do_immune_response_scores_MM.csv')

# 3.3 Compute quantitative descriptors of the TME

cell_fractions_MM <- compute_cell_fractions(RNA_tpm = RNA_tpm_seg1)
head(cell_fractions_MM)

#write.csv2(cell_fractions_MM,'output/RNA_2do_cell_fractions_MM.csv')


# By applying DoRothEA (Garcia-Alonso et al. 2019) method to TPM data from
# RNA-seq, the activity of 118 transcription factor (TF) can be inferred as follows:

tf_activities_MM <- compute_TF_activity(RNA_tpm = RNA_tpm_seg1)

head(tf_activities_MM[,1:5])

#write.csv2(tf_activities_MM,'output/RNA_2do_TFs_activities_MM.csv')


# Ligand-receptor (LR)
lrpair_weights_MM <- compute_LR_pairs(RNA_tpm = RNA_tpm_seg1,
                                   cancer_type = "pancan")
head(lrpair_weights_MM[,1:5])

#write.csv2(lrpair_weights_MM,'output/RNA_2do_lrpair_weights_MM.csv')

#Using the ligand-receptor weights as input, 169 cell-cell interaction scores can be derived as in the chunk below.

ccpair_scores_MM <- compute_CC_pairs(lrpairs = lrpair_weights_MM, 
                                  cancer_type = "pancan")

#> CC pairs computed
head(ccpair_scores_MM[,1:5])

#write.csv2(ccpair_scores_MM,'output/RNA_2do_ccpair_scores_MM.csv')

# 3.4 Obtain patients’ predictions of immune response

predictions_MM <- predict_immune_response(immunecells = cell_fractions_MM,
                                          #pathways = pathway_activities_MM,
                                          tfs = tf_activities_MM,
                                          lrpairs = lrpair_weights_MM,
                                          ccpairs = ccpair_scores_MM,
                                          cancer_type = cancer_type_MM, 
                                          verbose = TRUE)


# 3.6 score of likelihood of immune response

output_eval_no_resp_MM <- assess_immune_response(predictions_immune_response = predictions_MM,
                                                 RNA_tpm = RNA_tpm_seg1,
                                                 # TMB_values = TMB,
                                                 easier_with_TMB = "weighted_average",
                                                 weight_penalty = 0.5)

                                                 # TMB_values = TMB,
                                                 # easier_with_TMB = "weighted_average",
                                                 # weight_penalty = 0.5)

# Figure 1 output: Boxplot of patients’ easier score showing its distribution across the 10 different tasks.
# Figure 2 output: Scatterplot of patients’ prediction when combining easier score with tumor mutational 
#                  burden showing its distribution across the 10 different tasks.

output_eval_no_resp_MM[[1]]
output_eval_no_resp_MM[[2]]


# 3.6.1 Retrieve easier scores of immune response

easier_derived_scores_MM <- retrieve_easier_score(predictions_immune_response = predictions_MM,
                                                  # TMB_values = TMB,
                                                  easier_with_TMB = c("weighted_average", 
                                                                       "penalized_score"),
                                                  weight_penalty = 0.5)

head(easier_derived_scores_MM)
#write.csv(easier_derived_scores_MM,'output/easier_derived_scores_RNA_119MM.csv')

# 3.7 Interpret response to immunotherapy through systems biomarkers

output_biomarkers <- explore_biomarkers(immunecells = cell_fractions_MM,
                                        #pathways = pathway_activities_MM,
                                        lrpairs = lrpair_weights_MM,
                                        tfs = tf_activities_MM,
                                        ccpairs = ccpair_scores_MM,
                                        cancer_type = cancer_type_MM)
                                        #patient_response = patient_ICBresponse)

write.csv2(ccpair_scores_MM,'output/RNA_2do_ccpair_scores_MM.csv')
#=======================================================================


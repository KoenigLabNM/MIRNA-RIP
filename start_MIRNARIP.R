library(MIRNARIP)

###############################################################################
# CHANGE paths and file names to fit your data
###############################################################################

#global paths
basic_output_path = "/home/user/MIRNARIP/output/"
data_path = "/home/user/MIRNARIP/data/"

##############################################################################

#Gene and miRNA expression data files
gene_expression_file = paste(data_path,"RNAseqV2_expression_matrix_update_symbols.csv",sep="")
mirna_expression_file = paste(data_path,"miRNA_seq_expression_matrix.csv",sep="")

##############################################################################

#miRNA - target gene mapping file (i.e. from TarBase or other ressources)
mirna_target_gene_file = paste(data_path,"miRNA_target_gene_TarBase.csv",sep="")

# a mapping file of mature miRNA IDs to pre-miRNA IDs
# i.e. "hsa-miR-146a-5p"   "hsa-mir-146a"
mirna_mature_pre_mirna_mapping_file = paste(data_path,"total_mature_pre_mirna_mapping.csv",sep="")

# a mapping file of experimental miRNA IDs (i.e. the rownames in your miRNA expression matrix) and pre-mirRNA IDs
# i.e. "hsa-mir-486"   "hsa-mir-486-1"
exp_mirna_mapping_file = paste(data_path,"total_experimental_pre_mirna_mapping.csv",sep="")

# a file of miRNA candidates
mirna_candidates_file = paste(data_path,"mirna_candidates_example.txt",sep="")

##############################################################################
#Parameters

#Num cores
num_cores_wf_1 = 40
num_cores_wf_4_pwl = 40
num_cores_wf_5_lin = 40
num_cores_wf_8 = 40

#cutoff for correlation of miRNA and gene expression 
#default = 0
neg_cor_cutoff = 0

##Parameter for piecewise linear and linear model
#Number of cross-validations to perform
num_fold_cross_validation = 5

#constraints for PWL model
beta_m = 10
big_m = 1000

#miRNA expression quantiles
#default: 20% step quantiles
quantiles = c(0.0,0.2,0.4,0.6,0.8,1.0)
#number of miRNA quantiles i.e. how many steps in the model
num_mirna_expression_bin = length(quantiles)

#Parameter for GO enrichment
min_target_gene_num = 5
go_node_size = 5
top_go_sign_cutoff = 0.05
good_prediction_cutoff = 0.25

################################################################################
# Internal file name definitions. DO NOT CHANGE!
################################################################################

#Output file names
mirna_expression_z_score_file_name = "mirna_expression_intersect_z_score_matrix.csv"
gene_expression_z_score_file_name = "gene_expression_intersect_z_score_matrix.csv"
intersect_samples_file_name = "intersect_samples.txt"
mirna_gene_cor_file_name = "mirna_all_gene_correlation_matrix.csv"
exp_mirna_with_targets_and_data_file_name = "total_experimental_mirna_with_targets_and_data.txt"
combined_model_cor_summary_file_name = "mirna_combined_model_cor_summary.csv"

#File suffixes
mirna_quantile_file_suffix = "_mirna_expression_bin_matrix.csv"
model_prediction_file_suffix = "_prediction_correlation.csv"
model_result_paramter_file_suffix = "_result_parameter_matrix.csv"
cor_summary_file_suffix = "_prediction_cor_summary.csv"
global_cor_summary_file_suffix = "global_cor_summary.csv"
enrichment_file_suffix = "_enrichment.csv"

##############################################################################
##Fixed output path definitions

mirna_bin_dir_name = "mirna_bin/"
mirna_bin_path = paste(basic_output_path,mirna_bin_dir_name,sep="")

piecewise_linear_cor_dir_name = "model/mirna/piecewise_linear/cor/"
piecewise_linear_beta_parameter_dir_name = "model/mirna/piecewise_linear/beta/"

linear_cor_dir_name = "model/mirna/linear/cor/"
linear_beta_dir_name = "model/mirna/linear/beta/"

piecewise_linear_cor_summary_dir = "model/mirna/piecewise_linear/cor_summary/"
linear_cor_summary_dir = "model/mirna/linear/cor_summary/"

piecewise_linear_results_dir = "model/mirna/piecewise_linear/results/"
linear_results_dir = "model/mirna/linear/results/"

combined_model_dir = "model_combination/"
combined_model_enrichment_dir = "enrichment/GO/good_prediction_genes_background_all_target_genes/"

piecewise_linear_cor_path = paste(basic_output_path,piecewise_linear_cor_dir_name,sep="")
piecewise_linear_beta_path = paste(basic_output_path,piecewise_linear_beta_parameter_dir_name,sep="")

linear_cor_path = paste(basic_output_path,linear_cor_dir_name,sep="")
linear_beta_path = paste(basic_output_path,linear_beta_dir_name,sep="")

piecewise_linear_cor_summary_path = paste(basic_output_path,piecewise_linear_cor_summary_dir,sep="")
linear_cor_summary_path = paste(basic_output_path,linear_cor_summary_dir,sep="")

piecewise_linear_results_path = paste(basic_output_path,piecewise_linear_results_dir,sep="")
linear_results_path = paste(basic_output_path,linear_results_dir,sep="")

combined_model_path = paste(basic_output_path,combined_model_dir,sep="")
combined_model_enrichment_path = paste(basic_output_path,combined_model_dir,combined_model_enrichment_dir,sep="")

#files
mirna_expression_z_score_file = paste(basic_output_path,mirna_expression_z_score_file_name,sep="")
gene_expression_z_score_file = paste(basic_output_path,gene_expression_z_score_file_name,sep="")
intersect_samples_file = paste(basic_output_path,intersect_samples_file_name,sep="")
mirna_gene_cor_file = paste(basic_output_path,mirna_gene_cor_file_name,sep="")
exp_mirna_with_targets_and_data_file = paste(basic_output_path,exp_mirna_with_targets_and_data_file_name,sep="")
piecewise_linear_global_cor_summary_file = paste(piecewise_linear_results_path,global_cor_summary_file_suffix,sep="")
linear_global_cor_summary_file = paste(linear_results_path,global_cor_summary_file_suffix,sep="")
combined_model_cor_summary_file = paste(combined_model_path,combined_model_cor_summary_file_name,sep="")

##############################################################################

#create directories if they don't exist

if(!dir.exists(mirna_bin_path)){
  dir.create(mirna_bin_path)
}

if(!dir.exists(piecewise_linear_cor_path)){
  dir.create(piecewise_linear_cor_path, recursive=TRUE)
}

if(!dir.exists(piecewise_linear_beta_path)){
  dir.create(piecewise_linear_beta_path, recursive=TRUE)
}

if(!dir.exists(linear_cor_path)){
  dir.create(linear_cor_path, recursive=TRUE)
}

if(!dir.exists(linear_beta_path)){
  dir.create(linear_beta_path, recursive=TRUE)
}

if(!dir.exists(piecewise_linear_cor_summary_path)){
  dir.create(piecewise_linear_cor_summary_path, recursive=TRUE)
}

if(!dir.exists(linear_cor_summary_path)){
  dir.create(linear_cor_summary_path, recursive=TRUE)
}

if(!dir.exists(combined_model_path)){
  dir.create(combined_model_path, recursive=TRUE)
}

if(!dir.exists(combined_model_enrichment_path)){
  dir.create(combined_model_enrichment_path, recursive=TRUE)
}

if(!dir.exists(piecewise_linear_results_path)){
  dir.create(piecewise_linear_results_path, recursive=TRUE)
}

if(!dir.exists(linear_results_path)){
  dir.create(linear_results_path, recursive=TRUE)
}


##############################################################################
#Workflow

#1)

print("Step 1: start computing correlation...")
compute_mirna_gene_expression_correlation(mirna_expression_file,gene_expression_file,intersect_samples_file,num_cores_wf_1, mirna_gene_cor_file)
print("Step 1 finished")
print("----------------")

########################

#2)

print("Step 2: start computing z_scores...")
print("miRNA expression")
compute_z_score_matrix(mirna_expression_file, intersect_samples_file, mirna_expression_z_score_file)
print("gene expression")
compute_z_score_matrix(gene_expression_file, intersect_samples_file, gene_expression_z_score_file)
print("Step 2 finished")
print("----------------")

########################

#3)

print("Step 3: start computing SOS type 2 values...")
compute_mirna_expression_quantiles(mirna_expression_z_score_file, intersect_samples_file, mirna_bin_path, mirna_quantile_file_suffix, exp_mirna_with_targets_and_data_file, quantiles)
print("Step 3 finished")
print("----------------")

########################

#4)

print("Step 4: run piecewise linear model...")
run_piecewise_linear_model(num_cores_wf_4_pwl, num_fold_cross_validation, num_mirna_expression_bin, neg_cor_cutoff, beta_m, big_m, gene_expression_z_score_file, intersect_samples_file, mirna_bin_path, mirna_quantile_file_suffix, mirna_gene_cor_file, exp_mirna_mapping_file, mirna_mature_pre_mirna_mapping_file, mirna_candidates_file, mirna_target_gene_file, piecewise_linear_cor_path, piecewise_linear_beta_path, model_prediction_file_suffix, model_result_paramter_file_suffix)
print("Step 4 finished")
print("----------------")

########################

#5)

print("Step 5: run linear model...")
run_linear_model(num_cores_wf_5_lin, num_fold_cross_validation, neg_cor_cutoff, big_m, gene_expression_z_score_file, mirna_expression_z_score_file, intersect_samples_file, mirna_gene_cor_file, exp_mirna_mapping_file, mirna_mature_pre_mirna_mapping_file, mirna_candidates_file, mirna_target_gene_file, linear_cor_path, linear_beta_path, model_prediction_file_suffix, model_result_paramter_file_suffix)
print("Step 5 finished")
print("----------------")

########################

#6)

print("Step 6: summarize correlation files...")
print("Piecewise linear model:")
summarize_cor_files_by_mirna(piecewise_linear_cor_path, model_prediction_file_suffix, piecewise_linear_cor_summary_path, cor_summary_file_suffix, piecewise_linear_global_cor_summary_file)
print("Linear model:")
summarize_cor_files_by_mirna(linear_cor_path, model_prediction_file_suffix, linear_cor_summary_path, cor_summary_file_suffix, linear_global_cor_summary_file)
print("Step 6 finished")
print("----------------")

########################

#7)

print("Step 7: create combined model correlation files...")
create_combined_model_correlation_file(linear_global_cor_summary_file, piecewise_linear_global_cor_summary_file, combined_model_cor_summary_file)
print("Step 7 finished")
print("----------------")

########################

#8)

print("Step 8: perform GO enrichment analysis...")
perform_model_predictions_GO_enrichment_analysis(num_cores_wf_8, mirna_gene_cor_file, exp_mirna_mapping_file, mirna_mature_pre_mirna_mapping_file, mirna_target_gene_file, combined_model_cor_summary_file, combined_model_enrichment_path, enrichment_file_suffix, min_target_gene_num, go_node_size, top_go_sign_cutoff, neg_cor_cutoff, good_prediction_cutoff, mirna_candidates_file)
print("Step 8 finished")
print("----------------")

  
  
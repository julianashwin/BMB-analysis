standardise <- function(x){
  x <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  return(x)
}






cor_mat_plot <- function(df, qual_vars){
  # Create correlation matrix
  cormat_hh <- rcorr(as.matrix(df[,qual_vars]))
  # Coefs
  cormat_hh_coef <- cormat_hh$r
  cormat_hh_coef[upper.tri(cormat_hh_coef)] <- NA
  diag(cormat_hh_coef) <- NA
  cormat_hh_coef <- tibble(melt(cormat_hh_coef, na.rm = T, value.name = "coef")) %>%
    filter(Var1 %in% qual_vars & Var2 %in% qual_vars) %>%
    mutate(coef = round(coef,2))
  # P-values
  cormat_hh_pval <- cormat_hh$P
  cormat_hh_pval[upper.tri(cormat_hh_pval)] <- NA
  diag(cormat_hh_pval) <- NA
  cormat_hh_pval <- melt(cormat_hh_pval, na.rm = T, value.name = "pval") %>%
    filter(Var1 %in% qual_vars & Var2 %in% qual_vars) %>%
    mutate(alpha = case_when(pval < 0.1 ~ "sig", TRUE ~ "not")) %>%
    mutate(stars = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**", pval < 0.1 ~ "*", TRUE ~ ""))
  # Merge back 
  cormat_hh <- left_join(cormat_hh_coef, cormat_hh_pval, by = c("Var1", "Var2"))
  # Plot matrix
  corplot_hh <- ggplot(cormat_hh, aes(x = (Var2), y = fct_rev(Var1), fill = coef)) + theme_minimal() + 
    geom_tile(color = "white") + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) + 
    labs(x = "", y = "", fill = "Correlation") + 
    geom_text(aes(label = str_c(coef,stars), alpha = alpha), color = "black", size = 2.3) + 
    scale_alpha_manual(values = c("sig" = 1, "not" = 0.5), guide = 'none') +
    ggtitle("")
  
  return(corplot_hh)
  
}










import_data <- function(import_dir, experiment){
  
  if (experiment == "4"){
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_1/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_1/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_1/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE) 
  } else if (experiment == "4q"){
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_2q/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_2q/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_2q/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE) 
  } else if (experiment == "4qa"){
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_2qa/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_2qa/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_2qa/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE) 
  } else if (experiment == "5"){
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_1/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r2_sent_df <- read.csv(paste0(import_dir, "experiment_1/r2/enhanced_sent_data.csv"), stringsAsFactors = FALSE)
    r2_sent_df <- aggregate(r2_sent_df[,c(paste0(qual_vars, "_pred_proba"))], by = 
                              list(uid = r2_sent_df$uid), FUN = mean)
    names(r2_sent_df) <- str_replace(names(r2_sent_df), "_pred_proba", "")
    acts_preds <- names(r2_hh_df)[which(str_detect(names(r2_hh_df), "_act|_pred"))]
    r2_hh_df <- merge(r2_hh_df[,c("uid", "annotation_status", acts_preds)], r2_sent_df, by = "uid")
    
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_1/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_sent_df <- read.csv(paste0(import_dir, "experiment_1/r3/enhanced_sent_data.csv"), stringsAsFactors = FALSE)
    r3_sent_df <- aggregate(r3_sent_df[,c(paste0(qual_vars_r3, "_pred_proba"))], by = 
                              list(uid = r3_sent_df$uid), FUN = mean)
    names(r3_sent_df) <- str_replace(names(r3_sent_df), "_pred_proba", "")
    acts_preds <- names(r3_hh_df)[which(str_detect(names(r3_hh_df), "_act|_pred"))]
    r3_hh_df <- merge(r3_hh_df[,c("uid", "annotation_status", acts_preds)], r3_sent_df, by = "uid")
    
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_1/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE) 
    pool_sent_df <- read.csv(paste0(import_dir, "experiment_1/r2_r3/enhanced_sent_data.csv"), stringsAsFactors = FALSE) 
    pool_sent_df <- aggregate(pool_sent_df[,c(paste0(qual_vars, "_pred_proba"))], by = 
                                list(uid = pool_sent_df$uid), FUN = mean)
    names(pool_sent_df) <- str_replace(names(pool_sent_df), "_pred_proba", "")
    acts_preds <- names(pool_hh_df)[which(str_detect(names(pool_hh_df), "_act|_pred"))]
    pool_hh_df <- merge(pool_hh_df[,c("uid", "annotation_status", acts_preds)], pool_sent_df, by = "uid")
  } else if (experiment == "5q"){
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_2q/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r2_sent_df <- read.csv(paste0(import_dir, "experiment_2q/r2/enhanced_sent_data.csv"), stringsAsFactors = FALSE)
    r2_sent_df[is.na(r2_sent_df)] <- 0
    r2_sent_df <- aggregate(r2_sent_df[,c(paste0(qual_vars, "_pred_proba"))], by = 
                              list(uid = r2_sent_df$uid), FUN = mean)
    names(r2_sent_df) <- str_replace(names(r2_sent_df), "_pred_proba", "")
    acts_preds <- names(r2_hh_df)[which(str_detect(names(r2_hh_df), "_act|_pred"))]
    r2_hh_df <- merge(r2_hh_df[,c("uid", "annotation_status", acts_preds)], r2_sent_df, by = "uid")
    
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_2q/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_sent_df <- read.csv(paste0(import_dir, "experiment_2q/r3/enhanced_sent_data.csv"), stringsAsFactors = FALSE)
    r3_sent_df[is.na(r3_sent_df)] <- 0
    r3_sent_df <- aggregate(r3_sent_df[,c(paste0(qual_vars_r3, "_pred_proba"))], by = 
                              list(uid = r3_sent_df$uid), FUN = mean)
    names(r3_sent_df) <- str_replace(names(r3_sent_df), "_pred_proba", "")
    acts_preds <- names(r3_hh_df)[which(str_detect(names(r3_hh_df), "_act|_pred"))]
    r3_hh_df <- merge(r3_hh_df[,c("uid", "annotation_status", acts_preds)], r3_sent_df, by = "uid")
    
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_2q/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE) 
    pool_sent_df <- read.csv(paste0(import_dir, "experiment_2q/r2_r3/enhanced_sent_data.csv"), stringsAsFactors = FALSE) 
    pool_sent_df[is.na(pool_sent_df)] <- 0
    pool_sent_df <- aggregate(pool_sent_df[,c(paste0(qual_vars, "_pred_proba"))], by = 
                                list(uid = pool_sent_df$uid), FUN = mean)
    names(pool_sent_df) <- str_replace(names(pool_sent_df), "_pred_proba", "")
    acts_preds <- names(pool_hh_df)[which(str_detect(names(pool_hh_df), "_act|_pred"))]
    pool_hh_df <- merge(pool_hh_df[,c("uid", "annotation_status", acts_preds)], pool_sent_df, by = "uid")
  } else if (experiment == "5qa"){
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_2qa/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r2_sent_df <- read.csv(paste0(import_dir, "experiment_2qa/r2/enhanced_sent_data.csv"), stringsAsFactors = FALSE)
    r2_sent_df[is.na(r2_sent_df)] <- 0
    r2_sent_df <- aggregate(r2_sent_df[,c(paste0(qual_vars, "_pred_proba"))], by = 
                              list(uid = r2_sent_df$uid), FUN = mean)
    names(r2_sent_df) <- str_replace(names(r2_sent_df), "_pred_proba", "")
    acts_preds <- names(r2_hh_df)[which(str_detect(names(r2_hh_df), "_act|_pred"))]
    r2_hh_df <- merge(r2_hh_df[,c("uid", "annotation_status", acts_preds)], r2_sent_df, by = "uid")
    
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_2qa/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_sent_df <- read.csv(paste0(import_dir, "experiment_2qa/r3/enhanced_sent_data.csv"), stringsAsFactors = FALSE)
    r3_sent_df[is.na(r3_sent_df)] <- 0
    r3_sent_df <- aggregate(r3_sent_df[,c(paste0(qual_vars_r3, "_pred_proba"))], by = 
                              list(uid = r3_sent_df$uid), FUN = mean)
    names(r3_sent_df) <- str_replace(names(r3_sent_df), "_pred_proba", "")
    acts_preds <- names(r3_hh_df)[which(str_detect(names(r3_hh_df), "_act|_pred"))]
    r3_hh_df <- merge(r3_hh_df[,c("uid", "annotation_status", acts_preds)], r3_sent_df, by = "uid")
    
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_2qa/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE) 
    pool_sent_df <- read.csv(paste0(import_dir, "experiment_2qa/r2_r3/enhanced_sent_data.csv"), stringsAsFactors = FALSE) 
    pool_sent_df[is.na(pool_sent_df)] <- 0
    pool_sent_df <- aggregate(pool_sent_df[,c(paste0(qual_vars, "_pred_proba"))], by = 
                                list(uid = pool_sent_df$uid), FUN = mean)
    names(pool_sent_df) <- str_replace(names(pool_sent_df), "_pred_proba", "")
    acts_preds <- names(pool_hh_df)[which(str_detect(names(pool_hh_df), "_act|_pred"))]
    pool_hh_df <- merge(pool_hh_df[,c("uid", "annotation_status", acts_preds)], pool_sent_df, by = "uid")
  } else {
    r2_hh_df <- read.csv(paste0(import_dir, "experiment_",experiment,"/r2/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    r3_hh_df <- read.csv(paste0(import_dir, "experiment_",experiment,"/r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
    pool_hh_df <- read.csv(paste0(import_dir, "experiment_",experiment,"/r2_r3/enhanced_data_uid.csv"), stringsAsFactors = FALSE)
  }
  pool_hh_df$uid_round <- pool_hh_df$uid
  # Add variables to r2 and r3 where they weren't in the tree
  r3_hh_df$job_secular <- r3_hh_df$job_specific + r3_hh_df$home_based_work
  r3_hh_df$job_secular_act <- r3_hh_df$job_specific_act + r3_hh_df$home_based_work_act
  r3_hh_df$job_secular_pred <- r3_hh_df$job_specific_pred + r3_hh_df$home_based_work_pred
  
  # Create panel where possible
  panel_hh_df <- merge(r2_hh_df, r3_hh_df, by = "uid")
  names(panel_hh_df) <- str_replace_all(names(panel_hh_df), "\\.x", "_r2")
  names(panel_hh_df) <- str_replace_all(names(panel_hh_df), "\\.y", "_r3")
  # Merge with quant
  r2_hh_df <- merge(r2_hh_df, quant_df[,c("uid", quant_vars)], by = "uid")
  r3_hh_df <- merge(r3_hh_df, quant_df[,c("uid", quant_vars)], by = "uid")
  panel_hh_df <- merge(panel_hh_df, quant_df[,c("uid", quant_vars)], by = "uid")
  pool_hh_df$uid <- str_sub(pool_hh_df$uid,1,10)
  pool_hh_df <- merge(pool_hh_df, quant_df[,c("uid", quant_vars)], by = "uid")
  # Fix up pooled data
  pool_hh_df$eld_sex_r2[which(str_detect(pool_hh_df$uid_round, "R3"))] <- pool_hh_df$eld_sex_r3[which(str_detect(pool_hh_df$uid_round, "R3"))]
  pool_hh_df$eld_age_r2[which(str_detect(pool_hh_df$uid_round, "R3"))] <- pool_hh_df$eld_age_r3[which(str_detect(pool_hh_df$uid_round, "R3"))]
  pool_hh_df$num_child_r2[which(str_detect(pool_hh_df$uid_round, "R3"))] <- pool_hh_df$num_child_r3[which(str_detect(pool_hh_df$uid_round, "R3"))]
  pool_hh_df$parent_eduyears_r2[which(str_detect(pool_hh_df$uid_round, "R3"))] <- pool_hh_df$parent_eduyears_r3[which(str_detect(pool_hh_df$uid_round, "R3"))]
  pool_hh_df$parent_reledu_r2[which(str_detect(pool_hh_df$uid_round, "R3"))] <- pool_hh_df$parent_reledu_r3[which(str_detect(pool_hh_df$uid_round, "R3"))]
  pool_hh_df$int_trauma_exp_r2[which(str_detect(pool_hh_df$uid_round, "R3"))] <- pool_hh_df$int_trauma_exp_r3[which(str_detect(pool_hh_df$uid_round, "R3"))]
  # Extract annotated data
  r2_hh_annot <- r2_hh_df[which(r2_hh_df[,"annotation_status"]!="unannotated"),]
  r3_hh_annot <- r3_hh_df[which(r3_hh_df[,"annotation_status"]!="unannotated"),]
  pool_hh_annot <- pool_hh_df[which(pool_hh_df[,"annotation_status"]!="unannotated"),]
  
  
  # Remove r2 subscripts for r2 and pooled data
  names(r2_hh_annot) <- str_replace(names(r2_hh_annot), "_r2", "")
  names(r2_hh_df) <- str_replace(names(r2_hh_df), "_r2", "")
  names(pool_hh_annot) <- str_replace(names(pool_hh_annot), "_r2", "")
  names(pool_hh_df) <- str_replace(names(pool_hh_df), "_r2", "")
  # Remove r3 subscripts for r3 data
  names(r3_hh_annot) <- str_replace(names(r3_hh_annot), "_r3", "")
  names(r3_hh_df) <- str_replace(names(r3_hh_df), "_r3", "")
  
  
  # For experient 4 we use maximums rather than means
  if (str_detect(experiment, "4")){
    for (qual_var in qual_vars){
      # R2
      r2_hh_annot[,qual_var] <- as.numeric(r2_hh_annot[,qual_var] > 0)
      r2_hh_annot[,paste0(qual_var, "_pred")] <- as.numeric(r2_hh_annot[,paste0(qual_var, "_pred")] > 0)
      r2_hh_annot[,paste0(qual_var, "_act")] <- as.numeric(r2_hh_annot[,paste0(qual_var, "_act")] > 0)
      r2_hh_df[,qual_var] <- as.numeric(r2_hh_df[,qual_var] > 0)
      r2_hh_df[,paste0(qual_var, "_pred")] <- as.numeric(r2_hh_df[,paste0(qual_var, "_pred")] > 0)
      r2_hh_df[,paste0(qual_var, "_act")] <- as.numeric(r2_hh_df[,paste0(qual_var, "_act")] > 0)
      # R3
      r3_hh_annot[,qual_var] <- as.numeric(r3_hh_annot[,qual_var] > 0)
      r3_hh_annot[,paste0(qual_var, "_pred")] <- as.numeric(r3_hh_annot[,paste0(qual_var, "_pred")] > 0)
      r3_hh_annot[,paste0(qual_var, "_act")] <- as.numeric(r3_hh_annot[,paste0(qual_var, "_act")] > 0)
      r3_hh_df[,qual_var] <- as.numeric(r3_hh_df[,qual_var] > 0)
      r3_hh_df[,paste0(qual_var, "_pred")] <- as.numeric(r3_hh_df[,paste0(qual_var, "_pred")] > 0)
      r3_hh_df[,paste0(qual_var, "_act")] <- as.numeric(r3_hh_df[,paste0(qual_var, "_act")] > 0)
      # Pool
      pool_hh_annot[,qual_var] <- as.numeric(pool_hh_annot[,qual_var] > 0)
      pool_hh_annot[,paste0(qual_var, "_pred")] <- as.numeric(pool_hh_annot[,paste0(qual_var, "_pred")] > 0)
      pool_hh_annot[,paste0(qual_var, "_act")] <- as.numeric(pool_hh_annot[,paste0(qual_var, "_act")] > 0)
      pool_hh_df[,qual_var] <- as.numeric(pool_hh_df[,qual_var] > 0)
      pool_hh_df[,paste0(qual_var, "_pred")] <- as.numeric(pool_hh_df[,paste0(qual_var, "_pred")] > 0)
      pool_hh_df[,paste0(qual_var, "_act")] <- as.numeric(pool_hh_df[,paste0(qual_var, "_act")] > 0)
      # Panel 
      panel_hh_df[,paste0(qual_var,"_r2")] <- as.numeric(panel_hh_df[,paste0(qual_var,"_r2")] > 0)
      panel_hh_df[,paste0(qual_var,"_r3")] <- as.numeric(panel_hh_df[,paste0(qual_var,"_r3")] > 0)
    }
  }
  
  
  out_data <- list(r2_hh_df = r2_hh_df, r2_hh_annot = r2_hh_annot, 
                   r3_hh_df = r3_hh_df, r3_hh_annot = r3_hh_annot, 
                   pool_hh_df = pool_hh_df, pool_hh_annot = pool_hh_annot,
                   panel_hh_df = panel_hh_df)
  return(out_data)
  
}


rename_pool_quant <- function(pool_df){
  pool_df$eld_sex_r2[which(str_detect(pool_df$uid_round, "R3"))] <- pool_df$eld_sex_r3[which(str_detect(pool_df$uid_round, "R3"))]
  pool_df$eld_age_r2[which(str_detect(pool_df$uid_round, "R3"))] <- pool_df$eld_age_r3[which(str_detect(pool_df$uid_round, "R3"))]
  pool_df$num_child_r2[which(str_detect(pool_df$uid_round, "R3"))] <- pool_df$num_child_r3[which(str_detect(pool_df$uid_round, "R3"))]
  pool_df$parent_eduyears_r2[which(str_detect(pool_df$uid_round, "R3"))] <- pool_df$parent_eduyears_r3[which(str_detect(pool_df$uid_round, "R3"))]
  pool_df$parent_reledu_r2[which(str_detect(pool_df$uid_round, "R3"))] <- pool_df$parent_reledu_r3[which(str_detect(pool_df$uid_round, "R3"))]
  pool_df$int_trauma_exp_r2[which(str_detect(pool_df$uid_round, "R3"))] <- pool_df$int_trauma_exp_r3[which(str_detect(pool_df$uid_round, "R3"))]
  names(pool_df) <- str_replace(names(pool_df), "_r2", "")
  return(pool_df)
}


combine_w_quant <- function(enh_df, quant_df, round = "r2_r3"){
  
  enh_df$uid_round <- enh_df$uid
  enh_df$uid <- str_sub(enh_df$uid,1,10)
  enh_df <- merge(enh_df, quant_df, by = "uid")
  
  if (round == "r2_r3"){
    enh_df$eld_sex_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$eld_sex_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$eld_age_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$eld_age_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$num_child_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$num_child_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$parent_eduyears_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$parent_eduyears_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$parent_reledu_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$parent_reledu_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$int_trauma_exp_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$int_trauma_exp_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$int_sex_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$int_sex_r3[which(str_detect(enh_df$uid_round, "R3"))]
    enh_df$int_age_r2[which(str_detect(enh_df$uid_round, "R3"))] <- enh_df$int_age_r3[which(str_detect(enh_df$uid_round, "R3"))]
    names(enh_df) <- str_replace(names(enh_df), "_r2", "")
    
    enh_df$hh_in_r2 <- enh_df$hh_in
    
  } else if (round == "r2"){
    names(enh_df) <- str_replace(names(enh_df), "_r2", "")
    enh_df$hh_in_r2 <- enh_df$hh_in
  } else if (round == "r3"){
    names(enh_df) <- str_replace(names(enh_df), "_r3", "")
    enh_df$hh_in_r3 <- enh_df$hh_in
  } else {
    print("No valid round given")
  }
  
  return(enh_df)
  
}




compute_se <- function(input_df, annot, use_oob = FALSE, dec_place = 4){
  error_row <- data.frame(annot = annot)
  error_row[,c("Nh", "Nm", "mu_h", "mu_m", "sig2_h", "sig2_m", "sig2_eps", 
               "se_h", "se_m", "se_enh")] <- 0
  input_df <- input_df[which(!is.na(input_df[,annot])),]
  med_df <- filter(enh_df, bootstrap_run == "average")
  input_df <- filter(enh_df, bootstrap_run != "average")
  # Nh and Nm
  error_row$Nh <- length(unique(input_df$uid_round[which(input_df$annotation_status != "unannotated")]))
  error_row$Nm <- length(unique(input_df$uid_round[which(input_df$annotation_status == "unannotated")]))
  N_sum <- error_row$Nh + error_row$Nm
  # Separate out annotated and unannotated
  input_act_df <- input_df[which(input_df$annotation_status != "unannotated"),]
  input_pred_df <- input_df[which(input_df$annotation_status == "unannotated"),]
  # Variance of bootstrapped predictions
  error_row$sig2_m <- var(input_pred_df[,paste0(annot,"_pred")], na.rm = TRUE)
  # Variance of bootstrapped residuals
  if (use_oob){
    btstrp_residuals <- input_act_df[,paste0(annot,"_act")] - input_act_df[,paste0(annot,"_pred")]
    test_obs <- which(input_act_df$split=="test")
    btstrp_residuals <- 
      input_act_df[test_obs,paste0(annot,"_act")] - 
      input_act_df[test_obs,paste0(annot,"_pred")]
    tempvar <- var(input_act_df[test_obs,paste0(annot,"_act")], na.rm = TRUE)
  } else {
    btstrp_residuals <- med_df[,paste0(annot,"_act")] - med_df[,paste0(annot,"_pred")]
  }
  error_row$sig2_eps <- var(btstrp_residuals, na.rm=TRUE)
  # Variance of observed annotations in data 
  error_row$sig2_h <- var(input_act_df[,paste0(annot,"_act")], na.rm =T)
  # Mean in human sample
  error_row$mu_h <- mean(input_act_df[,paste0(annot,"_act")], na.rm =T)
  # Mean in machine sample
  error_row$mu_m <- mean(input_pred_df[,paste0(annot,"_pred")], na.rm =T)
  # standard error in human sample
  error_row$se_h <- sqrt(error_row$sig2_h/error_row$Nh)
  # standard error in machine sample
  error_row$se_m <- sqrt((error_row$sig2_m + error_row$sig2_eps)/error_row$Nm)
  # standard error in enhanced sample  
  error_row$se_enh <- sqrt((error_row$Nh*error_row$sig2_h + 
                              error_row$Nm*(error_row$sig2_m + error_row$sig2_eps))/
                             (N_sum^2))
  error_row[,names(error_row)[which(str_detect(names(error_row), "_"))]] <- 
    round(error_row[,names(error_row)[which(str_detect(names(error_row), "_"))]], digits = dec_place)
  
  return(error_row)
}






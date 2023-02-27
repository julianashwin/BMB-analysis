setwd("/Users/julianashwin/Documents/GitHub/BMB-analysis/")
rm(list=ls())

"
BMB data
"
library(tidyverse)
library(reshape2)
library(readxl)

source("src/analysis_fns.R")


bmb_df <- as_tibble(read_xlsx("data/BMB_annot_clean.xlsx", sheet = "BMB_annotations"))

# Separate df for annotated and unannotated docs
bmb_unannot <- filter(bmb_df, annotation_status == "unannotated")
bmb_annot <- filter(bmb_df, annotation_status == "annotated")%>%
  mutate(strict_direct_exp = replace_na(strict_direct_exp, 0), 
         broad_direct_exp = replace_na(broad_direct_exp, 0), 
         hearsay = replace_na(hearsay, 0))
write.csv(bmb_annot, "data/bmb_annotated.csv", row.names = F)
write.csv(bmb_unannot, "data/bmb_unannotated.csv", row.names = F)
#bmb_annot <- read.csv("data/bmb_annotated.csv")

# Import and merge quant data 
quant_df <- as_tibble(read.csv("data/hh_quant.csv"))
quant_df <- quant_df[,which(!str_detect(names(quant_df), "_r2"))]
names(quant_df) <- str_replace_all(names(quant_df), "_r3", "")

bmb_df <- left_join(bmb_df, quant_df, by = "uid")
bmb_annot <- left_join(bmb_annot, quant_df, by = "uid")


### Separate out the codes for each section
# A codes are about opinions/thoughts
A_codes <- names(bmb_df)[str_detect(names(bmb_df), "^A[0-9]+")]
# B codes are about positive experiences
B_codes <- names(bmb_df)[str_detect(names(bmb_df), "^B[0-9]+")]
# C codes are about negative experiences
C_codes <- names(bmb_df)[str_detect(names(bmb_df), "^C[0-9]+")]
# D codes are about government actions
D_codes <- names(bmb_df)[str_detect(names(bmb_df), "^D[0-9]+")]
# All qual codes
qual_vars <- c(A_codes, B_codes, C_codes, D_codes)
# Define the quant vars to use
quant_vars <- c("refugee", "int_sex", "int_age", "int_eduyears", "int_reledu",
                "int_trauma_exp", "num_child", "hh_cons_capita", "hh_income")


"
Look at sparsity of annotations in each section
"
sparsity_table <- crossing(codes= c("A", "B", "C", "D", "broad_direct_exp", "strict_direct_exp", "hearsay"), 
                           section = c("A", "B", "C", "D"), 
                           n_pos = NA, n_codes = NA, n_lines = NA, sparsity = NA)
for (ii in 1:nrow(sparsity_table)){
  sec <- sparsity_table$section[ii]
  codes <- sparsity_table$codes[ii]
  if (nchar(codes) == 1){
    code_names <- names(bmb_annot)[str_detect(names(bmb_annot), str_c("^",codes,"[0-9]+"))]
  } else {
    code_names <- codes
  }
  # Count all codes and lines in that section
  sparsity_table$n_pos[ii] <- sum(filter(bmb_annot, section == sec)[,code_names], na.rm = T)
  sparsity_table$n_codes[ii] <- length(code_names)
  sparsity_table$n_lines[ii] <- nrow(filter(bmb_annot, section == sec))
  # Calculate sparsity
  sparsity_table$sparsity[ii] <-sparsity_table$n_pos[ii]/
    (sparsity_table$n_codes[ii]*sparsity_table$n_lines[ii])
}

## Codes occus more in their own Section, but also quite frequently in other Sections as well
sparsity_table %>%
  mutate(codes = case_when(codes == "broad_direct_exp" ~ "Broad",
                           codes == "strict_direct_exp" ~ "Strict",
                           codes == "hearsay" ~ "Hearsay", TRUE ~ codes)) %>%
  mutate(codes = factor(codes, levels = c("A", "B", "C", "D", "Broad", "Strict", "Hearsay"), ordered = T)) %>%
  ggplot(aes(x = codes, y = sparsity)) + theme_bw() + 
  geom_point(aes(color = section)) + 
  labs(x = "Codes", y = "Average Code Sparsity", color = "Section")
ggsave("figures/human_only/sparsity_codes_sections.pdf", width = 5, height = 3)


meds_hh <- bmb_annot %>% 
  select(-section, -annotation_status, -Q_bangla, -A_bangla, -Q_en, -A_en) %>%
  group_by(uid) %>%
  summarise(across(where(is.numeric), mean),
            across(where(is.character), median))
meds_section <- bmb_annot %>% 
  select(-annotation_status, -Q_bangla, -A_bangla, -Q_en, -A_en) %>%
  group_by(uid, section) %>%
  summarise(across(where(is.numeric), mean),
            across(where(is.character), median))

# Keep only annotations in the "correct" section. 
meds_secA <- meds_section %>% filter(section == "A") %>%
  select(uid, A1_no_opinion:A16_religious_fraternity)
meds_secB <- meds_section %>% filter(section == "B") %>%
  select(uid, B1_no_benefits:B7_provided_food)
meds_secC <- meds_section %>% filter(section == "C") %>%
  select(uid, C1_no_problems:C17_overcrowded)
meds_secD <- meds_section %>% filter(section == "D") %>%
  select(uid, D1_no_actions:D13_increased_mobility)
meds_section <- meds_hh %>%
  select(uid, strict_direct_exp, broad_direct_exp, hearsay,
         refugee:int_trauma_witness) %>%
  left_join(meds_secA, by = "uid") %>%
  left_join(meds_secB, by = "uid") %>%
  left_join(meds_secC, by = "uid") %>%
  left_join(meds_secD, by = "uid")

rm(meds_secA, meds_secB, meds_secC, meds_secD)


"
Correlation matrix
"
library(Hmisc)

cor_mat_plot(meds_hh[,qual_vars], qual_vars) + 
  ggtitle("Whole Interview")
ggsave("figures/human_only/cormat_whole.pdf", width = 18, height = 12)

cor_mat_plot(meds_section[,qual_vars], qual_vars) + 
  ggtitle("'Correct' Section")
ggsave("figures/human_only/cormat_section.pdf", width = 18, height = 12)


"
Regressions
"
library(stargazer)
## What explains the A codes?
overall_string <-""
for (annot in A_codes){
  # Whole interview 
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, quant_vars), collapse = "+")))
  model1 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c(C_codes, quant_vars), collapse = "+")))
  model2 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, C_codes, quant_vars), collapse = "+")))
  model3 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars), collapse = "+")))
  model4 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars, "distancekmfromcamps"), collapse = "+")))
  model5 <- lm(fm, meds_hh)
  
  # Section only
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, quant_vars), collapse = "+")))
  model6 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c(C_codes, quant_vars), collapse = "+")))
  model7 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, C_codes, quant_vars), collapse = "+")))
  model8 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars), collapse = "+")))
  model9 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars, "distancekmfromcamps"), collapse = "+")))
  model10 <- lm(fm, meds_section)
  
  
  models <-list(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)
  stargazer(models,
            title = str_replace_all(annot, "_","."), df = F, table.placement = "H",
            column.sep.width = "4pt", font.size= "tiny",
            out=paste0("figures/human_only/regs/",annot,".tex"))
  rm(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)
  
  overall_string <- paste0(overall_string, "\\input{regs/", annot, ".tex}", " ")
}
writeLines(overall_string)





## What explains the B and C codes?
overall_string <-""
for (annot in c(B_codes, C_codes)){
  # Whole interview 
  fm <- as.formula(str_c(annot, "~", paste(c(quant_vars), collapse = "+")))
  model1 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             quant_vars), collapse = "+")))
  model2 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             quant_vars, "distancekmfromcamps"), collapse = "+")))
  model3 <- lm(fm, meds_hh)
  
  # Section only
  fm <- as.formula(str_c(annot, "~", paste(c(quant_vars), collapse = "+")))
  model4 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             quant_vars), collapse = "+")))
  model5 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             quant_vars,  "distancekmfromcamps"), collapse = "+")))
  model6 <- lm(fm, meds_section)
  
  
  models <-list(model1, model2, model3, model4, model5, model6)
  stargazer(models,
            title = str_replace_all(annot, "_","."), df = F, table.placement = "H",
            column.sep.width = "4pt", font.size= "tiny",
            out=paste0("figures/human_only/regs/",annot,".tex"))
  rm(model1, model2, model3, model4, model5, model6)
  
  overall_string <- paste0(overall_string, "\\input{regs/", annot, ".tex}", " ")
}
writeLines(overall_string)



## What explains the D codes?
overall_string <-""
for (annot in D_codes){
  # Whole interview 
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, quant_vars), collapse = "+")))
  model1 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c(C_codes, quant_vars), collapse = "+")))
  model2 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, C_codes, quant_vars), collapse = "+")))
  model3 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars), collapse = "+")))
  model4 <- lm(fm, meds_hh)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars, "distancekmfromcamps"), collapse = "+")))
  model5 <- lm(fm, meds_hh)
  
  # Section only
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, quant_vars), collapse = "+")))
  model6 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c(C_codes, quant_vars), collapse = "+")))
  model7 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c(B_codes, C_codes, quant_vars), collapse = "+")))
  model8 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars), collapse = "+")))
  model9 <- lm(fm, meds_section)
  fm <- as.formula(str_c(annot, "~", paste(c("strict_direct_exp", "broad_direct_exp", "hearsay", 
                                             B_codes, C_codes, quant_vars, "distancekmfromcamps"), collapse = "+")))
  model10 <- lm(fm, meds_section)
  
  
  models <-list(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)
  stargazer(models,
            title = str_replace_all(annot, "_","."), df = F, table.placement = "H",
            column.sep.width = "4pt", font.size= "tiny",
            out=paste0("figures/human_only/regs/",annot,".tex"))
  rm(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)
  
  overall_string <- paste0(overall_string, "\\input{regs/", annot, ".tex}", " ")
}
writeLines(overall_string)

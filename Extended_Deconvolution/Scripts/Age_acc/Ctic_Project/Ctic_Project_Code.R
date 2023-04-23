
#Ctic/Horvath/Hannum/PhenoAge predction
load("/Users/chenjiqing/Public/bladder cancer/Extended Deconvolution/Script/Age_acc/Ctic_Project/elastic_net_model1.RDATA")

# bladder cnacer beta
beta <- readRDS("/Users/chenjiqing/Public/bladder cancer/raw data/AllBC_Combat_Betas2.rds") # 833378 CpGs

library(dplyr)

#CticAge <- model %>% predict(t(beta)) 
#beta is the processed methylation beta matrix, but after QC some CpGs were filtered out. So we cannot estimate.
#CticAge2 <- as.data.frame(CticAge)

# load pheno
load("/Users/chenjiqing/Public/bladder cancer/Basic Analysis/data_clean2.RData")
rm(phenotype1,Recu)
# weird data
weird = phenotype2[which(phenotype2$ddstat=="dead" & phenotype2$RecurYN=="yes" &
                           phenotype2$survreftmd<phenotype2$recsurgtmd),]
# remove weird data
AllBC <- phenotype2 %>% filter(Sample_Name != weird$Sample_Name)
# cell proportion
load("/Users/chenjiqing/Public/bladder cancer/Extended Deconvolution/Bladder_Extended_cell_proportion2.RData")

AllBC_immune <- AllBC %>%
  mutate(Age = refage, Sex = sex, Tumor_grade = grade2, 
         Smoking_status = smkstat2, BCG = ImmunoTx) %>%
  dplyr::select(Sample_Name,Age,Sex,Tumor_grade,Smoking_status,BCG,TenDead, TenYearSurv, TenYearRFS, TenRFS, muscinv) %>%
  left_join(Extended2) %>%
  mutate(T_cell = CD4mem + CD4nv + CD8mem + CD8nv + Treg, # not show this one
         Lymph = Bmem + Bnv + NK + T_cell) %>%
  mutate(NLR = ifelse(Lymph == 0,Neu/0.01, Neu/Lymph)) %>%
  dplyr::select(-T_cell,-Lymph)
AllBC_immune_win <- AllBC_immune %>%
  mutate(Bas = ifelse(Bas <= quantile(AllBC_immune$Bas, probs = 0.98),
                      Bas,quantile(AllBC_immune$Bas, probs = 0.98)),
         Bmem = ifelse(Bmem <= quantile(AllBC_immune$Bmem, probs = 0.98),
                       Bmem,quantile(AllBC_immune$Bmem, probs = 0.98)),
         Bnv = ifelse(Bnv <= quantile(AllBC_immune$Bnv, probs = 0.98),
                      Bnv,quantile(AllBC_immune$Bnv, probs = 0.98)),
         CD4mem = ifelse(CD4mem <= quantile(AllBC_immune$CD4mem, probs = 0.98),
                         CD4mem,quantile(AllBC_immune$CD4mem, probs = 0.98)),
         CD4nv = ifelse(CD4nv <= quantile(AllBC_immune$CD4nv, probs = 0.98),
                        CD4nv,quantile(AllBC_immune$CD4nv, probs = 0.98)),
         CD8mem = ifelse(CD8mem <= quantile(AllBC_immune$CD8mem, probs = 0.98),
                         CD8mem,quantile(AllBC_immune$CD8mem, probs = 0.98)),
         CD8nv = ifelse(CD8nv <= quantile(AllBC_immune$CD8nv, probs = 0.98),
                        CD8nv,quantile(AllBC_immune$CD8nv, probs = 0.98)),
         Eos = ifelse(Eos <= quantile(AllBC_immune$Eos, probs = 0.98),
                      Eos,quantile(AllBC_immune$Eos, probs = 0.98)),
         Mono = ifelse(Mono <= quantile(AllBC_immune$Mono, probs = 0.98),
                       Mono,quantile(AllBC_immune$Mono, probs = 0.98)),
         Neu = ifelse(Neu <= quantile(AllBC_immune$Neu, probs = 0.02),
                      quantile(AllBC_immune$Neu, probs = 0.02),Neu),
         NK = ifelse(NK <= quantile(AllBC_immune$NK, probs = 0.98),
                     NK,quantile(AllBC_immune$NK, probs = 0.98)),
         Treg = ifelse(Treg <= quantile(AllBC_immune$Treg, probs = 0.98),
                       Treg,quantile(AllBC_immune$Treg, probs = 0.98)),
         NLR = ifelse(NLR <= quantile(AllBC_immune$NLR, probs = 0.98),
                      NLR,quantile(AllBC_immune$NLR, probs = 0.98)))

colnames(AllBC_immune_win)[1]<-"Sample_ID"

#CticAge2 <- CticAge2 %>% rownames_to_column(., var = "Sample_ID")

library(ENmix)

Horvath_Age<-methyAge(beta,type="hovath",normalize=FALSE,nCores=1)
Hannum_Age<-methyAge(beta,type="hannum",normalize=FALSE,nCores=1)
Pheno_Age<-methyAge(beta,type="phenoAge",normalize=FALSE,nCores=1)
DNAmAge<-inner_join(Horvath_Age,Hannum_Age,by = "SampleID")
DNAmAge<-inner_join(DNAmAge,Pheno_Age,by = "SampleID")
colnames(DNAmAge)[1]<-"Sample_ID"
DNAmAge[,2]<-as.numeric(DNAmAge[,2])
DNAmAge[,3]<-as.numeric(DNAmAge[,3])
DNAmAge[,4]<-as.numeric(DNAmAge[,4])
#DNAmAge<-inner_join(DNAmAge,CticAge2,by = "Sample_ID")

AllBC_immune_Age <- AllBC_immune_win %>% 
                    left_join(DNAmAge) %>%
                    mutate(hovath_acc = mAge_Hovath - Age,
                           hannum_acc = mAge_Hannum -Age,
                           pheno_acc = PhenoAge - Age)
save(AllBC_immune_Age, file = "/Users/chenjiqing/Public/bladder cancer/Extended Deconvolution/Script/Age_acc/Ctic_Project/AllBC_Age_acc2.RData")

# Age Acc distribution
load("/Users/chenjiqing/Public/bladder cancer/Extended Deconvolution/Script/Age_acc/Ctic_Project/AllBC_Age_acc2.RData")

NMIBC_immune_Age1 <- AllBC_immune_Age %>% 
                     filter(muscinv == "no") %>%
                     select(hovath_acc, hannum_acc, pheno_acc) %>%
                     gather(key = Acc_type, value = Acc)

Distr1 <- ggplot(NMIBC_immune_Age1, aes(x = Acc, by = Acc_type, 
                                        color = Acc_type, ..scaled..)) + 
          geom_density(alpha = 0.5) + 
          labs(x = "Age Acceleration (years)", y = "",colour = "") +
          theme_classic() +
          scale_color_nejm() +
          theme(legend.position = "bottom",
                text=element_text(size=12, face = "bold"),
                plot.margin = unit(c(6,0,6,0), "pt"))
# Age & Methylation Age distribution
NMIBC_immune_Age2 <- AllBC_immune_Age %>%
                     filter(muscinv == "no") %>%
                     select(Age, mAge_Hovath, mAge_Hannum, PhenoAge) %>%
                     gather(key = Age_type, value = mAge)

Distr2 <- ggplot(NMIBC_immune_Age2, aes(x = mAge, by = Age_type, 
                                        color = Age_type, ..scaled..)) + 
          geom_density(alpha = 0.5) + 
          labs(x = "Age (years)", y = "",colour = "") +
          theme_classic() +
          scale_color_nejm() +
          theme(legend.position = "bottom",
          text=element_text(size=12, face = "bold"),
          plot.margin = unit(c(6,0,6,0), "pt"))

Title = ggdraw() + draw_label("Age Distribution", x = 0, y = 1, vjust = 1, hjust = 0, size = 14, fontface = 'bold')
Comb_Plot = plot_grid(Title,Distr1, Distr2, ncol = 1, rel_heights = c(0.07,1, 1))
ggsave(filename = paste0("/Users/chenjiqing/Public/bladder cancer/Extended Deconvolution/Plot/Age and Age_Acc distribution.pdf"), plot = Comb_Plot, height = 9, width = 7.5)

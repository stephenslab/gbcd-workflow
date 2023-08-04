### this script performs survival analysis on PDAC bulk RNA-seq data using shared GEPs from gbcd and literature-derived subtype signatures (Fig. 4)
setwd("pdac")

### load in the required packages
library(Matrix)
library(readxl)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Seurat)
library(survival)
library(My.stepwise)
library(survminer)



########################################## load in the gene signatures for shared GEPs identified by gbcd ##########################################
### load in gbcd fit
fit.gbcd <- readRDS("pdac_gbcd_paper.rds")

### specify the indices of shared GEPs
k.idx1 <- c(3, 30, 7, 4, 23, 12, 5, 21, 54, 15, 16, 11, 19, 28)

### take the top driving genes for each shared GEP identified by gbcd
gene.lists <- list(NULL)
for(k in 1:length(k.idx1)){
  gene.lists[[k]] <- rownames(fit.gbcd$F$lfc)[fit.gbcd$F$lfc[, k.idx1[k]] > pmax(quantile(fit.gbcd$F$lfc[, k.idx1[k]], 1-200/nrow(fit.gbcd$F$lfc)), log2(1.5))]
}
names(gene.lists) <- paste0("GEP", 1:length(k.idx1))



########################################### load in the literature-derived subtype signatures ################################################
### load in signature gene sets for the subtype programs identified from Raghavan et al.
gene.sig <- read_xlsx("raghavan.xlsx", sheet=1, skip=3)
basal.r <- gene.sig$Gene[1:200]
gene.sig <- read_xlsx("raghavan.xlsx", sheet=2, skip=3)
classical.r <- gene.sig$Gene[1:200]

### load in the gene signatures for the subtype programs identified from Moffit et al.
gene.sig <- read_excel("moffit.xlsx")
basal.m <- gene.sig$symbol[order(gene.sig$F6_BasalLike, decreasing = TRUE)][1:200]
classical.m <- gene.sig$symbol[order(gene.sig$F8_Classical, decreasing = TRUE)][1:200]

### load in the gene signatures for the subtype programs identified from Chan-Seng-Yue et al.
gene.sig <- read_excel("chan_seng_yue.xlsx", skip=1, sheet=4)
classical1.c <- gene.sig$`Sig. 1 genes`[1:200]
basal1.c <- gene.sig$`Sig. 2 genes`[1:200]
classical2.c <- gene.sig$`Sig. 6 genes`[1:200]
basal2.c <- gene.sig$`Sig. 10 genes`[1:200]

### load in the gene signatures for the subtype programs identified from Hwang et al.
gene.sig <- read_excel("hwang.xlsx", skip=3, sheet=2)
classical.h <- gene.sig$`Classical-like`[1:200]
basaloid <- gene.sig$Basaloid[1:200]
squamoid <- gene.sig$Squamoid[1:200]
mesenchymal <- gene.sig$Mesenchymal[1:200]

### prepare the gene lists
gene.classical <- list(NULL)
gene.classical[[1]] <- classical.m
gene.classical[[2]] <- classical1.c
gene.classical[[3]] <- classical2.c
gene.classical[[4]] <- classical.r
gene.classical[[5]] <- classical.h
names(gene.classical) <- c("Classical_Moffit", "Classical_A_CK", "Classical_B_CK", "Classical_Raghavan", "Classical_Hwang")

gene.basal <- list(NULL)
gene.basal[[1]] <- basal.m
gene.basal[[2]] <- basal1.c
gene.basal[[3]] <- basal2.c
gene.basal[[4]] <- basal.r
gene.basal[[5]] <- basaloid
gene.basal[[6]] <- squamoid
gene.basal[[7]] <- mesenchymal
names(gene.basal) <- c("Basal_Moffit", "Basal_A_CK", "Basal_B_CK", "Basal_Raghavan", "Basaloid_Hwang", "Squamoid_Hwang", "Mesenchymal_Hwang")

### combine the gene lists
gene.lists <- c(gene.lists, gene.classical, gene.basal)



##################################### compute scores for all gene signatures in each patient from TCGA cohort ########################################
### load in data 
data1 <- readRDS("survival/tcga.rds")

### calculate the patient scores for gene signatures
dat <- data1$data.norm
tcga <- CreateSeuratObject(counts = t(dat), project = "tcga")
for(k in 1:length(gene.lists)){
  tcga <- AddModuleScore(object = tcga, features = list(gene.lists[[k]]), name = names(gene.lists)[k])
}
score1 <- tcga@meta.data[, -c(1:3)]



##################################### compute scores for all gene signatures in each patient from CPTAC cohort ########################################
### load in data
data2 <- readRDS("survival/cptac.rds")

### calculate the patient scores for gene signatures
dat <- data2$data.norm
cptac <- CreateSeuratObject(counts = t(dat), project = "cptac")
for(k in 1:length(gene.lists)){
  cptac <- AddModuleScore(object = cptac, features = list(gene.lists[[k]]), name = names(gene.lists)[k])
}
score2 <- cptac@meta.data[, -c(1:3)]



##################################### compute scores for all gene signatures in each patient from QCMG cohort ########################################
### load in data
data3 <- readRDS("survival/qcmg.rds")

### calculate the patient scores for gene signatures
dat <- data3$data.norm
genes <- intersect(colnames(dat), rownames(fit.gbcd$F$lfc))
dat <- dat[, genes]
qcmg <- CreateSeuratObject(counts = t(dat), project = "qcmg")
for(k in 1:length(gene.lists)){
  qcmg <- AddModuleScore(object = qcmg, features = list(gene.lists[[k]]), name = names(gene.lists)[k])
}
score3 <- qcmg@meta.data[, -c(1:3)]



############################### compute scores for all gene signatures in each patient from Pancurx resectable cohort ##################################
### load in data
data4 <- readRDS("survival/pancurx_resectable.rds")

### calculate the patient scores for gene signatures
dat <- data4$data.norm
genes <- intersect(colnames(dat), rownames(fit.gbcd$F$lfc))
dat <- dat[, genes]
pancurx <- CreateSeuratObject(counts = t(dat), project = "pancurx")
for(k in 1:length(gene.lists)){
  pancurx <- AddModuleScore(object = pancurx, features = list(gene.lists[[k]]), name = names(gene.lists)[k])
}
score4 <- pancurx@meta.data[, -c(1:3)]



######################################## preprocess and clean the data from all cohorts for survival analysis #########################################
### tcga cohort
info1 <- data1$metadata
futime1 <- as.numeric(info1$`Follow up days`)
event1 <- ifelse(info1$`Follow up vital status`=="Dead", 1, 0)
age1 <- info1$`Age at initial pathologic diagnosis`
sex1 <- info1$Gender
sex1[info1$Gender=="female"] <- "Female"
sex1[info1$Gender=="male"] <- "Male"
stage1 <- info1$`AJCC pathologic tumor stage`
stage1[grep("iia", info1$`AJCC pathologic tumor stage`)] <- "I-IIa"
stage1[grep("iib", info1$`AJCC pathologic tumor stage`)] <- "IIb"
stage1[grep("iii", info1$`AJCC pathologic tumor stage`)] <- "III-higher"
histo1 <- rep(0, nrow(info1))
histo1[grep("Ductal", info1$`Histological type by RHH`)] <- 1

### cptac cohort
info2 <- data2$metadata
idx.rm <- info2$pathologic_staging_distant_metastasis_pm=="pM1"
info2 <- info2[!idx.rm,]
score2 <- score2[!idx.rm,]
futime2 <- as.numeric(info2$follow_up_days)
event2 <- rep(NA, nrow(info2))
event2[grep("Living", info2$vital_status)] <- 0
event2[grep("Deceased", info2$vital_status)] <- 1
age2 <- info2$age
sex2 <- info2$sex
stage2 <- info2$tumor_stage_pathological
stage2[stage2=="Stage IA" | stage2=="Stage IB" | stage2=="Stage IIA"] <- "I-IIa"
stage2[stage2=="Stage IIB"] <- "IIb"
stage2[stage2=="Stage III" | stage2=="Stage IV"] <- "III-higher"
histo2 <- rep(0, nrow(info2))
histo2[info2$histology_diagnosis=="PDAC"] <- 1

### qcmg cohort
info3 <- data3$metadata
futime3 <- round(as.numeric(info3$`Length of Follow Up (months)`)*(365/12))
event3 <- rep(NA, nrow(info3))
event3[grep("Alive", info3$Status)] <- 0
event3[grep("Deceased", info3$Status)] <- 1
age3 <- as.numeric(info3$`Age at Diagnosis in Years`)
sex3 <- info3$Gender
stage3 <- info3$`AJCC Pathology Stage`
stage3[stage3=="IA" | stage3=="IB" | stage3=="IIA"] <- "I-IIa"
stage3[stage3=="IIB"] <- "IIb"
stage3[stage3=="III" | stage3=="IV"] <- "III-higher"
histo3 <- rep(0, nrow(info3))
histo3[grep("Ductal", info3$HistoSubtype)] <- 1

### pancurx resectable cohort
info4 <- data4$metadata
info4 <- info4[!is.na(info4$OSFromSurgery), ]
info4 <- info4[info4$OSFromSurgery!="Unknown" & info4$Age.at.Diagnosis!="NA" & info4$PathologicalStagingMetsM!="M1 Distant metastasis", ]
score4 <- score4[info4$PCSI, ]
futime4 <- as.numeric(info4$OSFromSurgery)
event4 <- ifelse(info4$Deceased=="TRUE", 1, 0)
age4 <- as.numeric(info4$Age.at.Diagnosis)
sex4 <- ifelse(info4$Sex=="M", "Male", "Female")
stage4 <- info4$PathologicalStaging
stage4[stage4=="IA" | stage4=="IB" | stage4=="IIA"] <- "I-IIa"
stage4[stage4=="IIB"] <- "IIb"
stage4[stage4=="III" | stage4=="IV" | stage4=="I - III"] <- "III-higher"
histo4 <- rep(0, nrow(info4))
histo4[info4$HistologicalType=="Ductal adenocarcinoma"] <- 1


### normalize the GEP scores separately for each cohort
colnames(score1) <- names(gene.lists)
colnames(score2) <- names(gene.lists)
colnames(score3) <- names(gene.lists)
colnames(score4) <- names(gene.lists)
score1 <- scale(score1)
score2 <- scale(score2)
score3 <- scale(score3)
score4 <- scale(score4)

### combine the data from cohorts
score <- rbind(score1, score2, score3, score4)
surdata <- data.frame(futime=c(futime1, futime2, futime3, futime4), event=c(event1, event2, event3, event4),
                      pdac=c(histo1, histo2, histo3, histo4), stage=c(stage1, stage2, stage3, stage4),
                      age=c(age1, age2, age3, age4), sex=c(sex1, sex2, sex3, sex4))
data.combined <- cbind(surdata, score)

### remove the observations which are not PDAC by histology or have missing values
idx.keep <- rowSums(is.na(data.combined))==0
data.combined <- data.combined[idx.keep,]
idx.keep <- rowSums(data.combined=="NA")==0
data.combined <- data.combined[idx.keep,]
idx.keep <- data.combined$pdac==1
data.combined <- data.combined[idx.keep,]
data.combined$sex <- factor(data.combined$sex)
data.combined$stage <- factor(data.combined$stage)
saveRDS(data.combined, "survival/combined_bulk.rds")



################# perform survival analysis on bulk RNA-seq data using shared GEPs and literature-derived subtype signatures (Fig. 4) ##################
### perform stepwise variable selection in a proportional hazard model while always including age, sex and tumor stage in the model
data.combined$male <- ifelse(data.combined$sex=="Male", 1, 0)
data.combined$stageIIb <- ifelse(data.combined$stage=="IIb", 1, 0)
data.combined$stageIII <- ifelse(data.combined$stage=="III-higher", 1, 0)
fit.surv <- My.stepwise.coxph(Time = "futime", Status = "event", variable.list = colnames(data.combined)[7:32],
                              in.variable = c("age", "male", "stageIIb", "stageIII"), data = data.combined)

### plot Kaplan Meier curves for dichotomized GEP14
data.combined$hypoxia_ISR <- as.factor(ifelse(data.combined$GEP14 > median(data.combined$GEP14), "high", "low"))
p <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined),
                legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", 
                pval = TRUE, pval.size = 8, pval.coord = c(1800, 0.8), conf.int = FALSE, surv.median.line = "hv", break.x.by = 250,
                ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival")
p <- ggpar(p, font.legend = c(16))
p

### plot Kaplan Meier curves for dichotomoized GEP14 within patient subgroups stratified by tumor stage and dichotomized GEP4 that is basal-related
data.combined$subtype <- as.factor(ifelse(data.combined$GEP4 > median(data.combined$GEP4), "basal", "classical"))
p1 <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined[data.combined$subtype=="basal" & data.combined$stage=="I-IIa",]),
                 legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", conf.int = FALSE, 
                 pval = TRUE, pval.size = 6, pval.coord = c(1800, 0.9), surv.median.line = "hv", xlim=c(0, max(data.combined$futime)), break.x.by = 250,
                 ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival", 
                 title = "GEP4-high (basal), tumor stage IIA or lower") 
p1$plot <- p1$plot + 
  ggplot2::annotate("text", x = 500, y = 0.2, 
  label = paste0("n = ", sum(data.combined$subtype=="basal" & data.combined$stage=="I-IIa" & data.combined$hypoxia_ISR=="high")), 
  size = 6, color = "#F8766D") + 
  ggplot2::annotate("text", x = 800, y = 0.6, 
  label = paste0("n = ", sum(data.combined$subtype=="basal" & data.combined$stage=="I-IIa" & data.combined$hypoxia_ISR=="low")), 
  size = 6, color = "#00BFC4")
p1 <- ggpar(p1, font.legend = c(16), font.main = c(12), font.xtickslab = c(9)) 

p2 <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined[data.combined$subtype=="basal" & data.combined$stage=="IIb",]),
                 legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", conf.int = FALSE, 
                 pval = TRUE, pval.size = 6, pval.coord = c(1800, 0.9), surv.median.line = "hv", xlim=c(0, max(data.combined$futime)), break.x.by = 250,
                 ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival", 
                 title = "GEP4-high (basal), tumor stage IIB") 
p2$plot <- p2$plot + 
  ggplot2::annotate("text", x = 200, y = 0.3, 
  label = paste0("n = ", sum(data.combined$subtype=="basal" & data.combined$stage=="IIb" & data.combined$hypoxia_ISR=="high")), 
  size = 6, color = "#F8766D") + 
  ggplot2::annotate("text", x = 600, y = 0.6, 
  label = paste0("n = ", sum(data.combined$subtype=="basal" & data.combined$stage=="IIb" & data.combined$hypoxia_ISR=="low")),
  size = 6, color = "#00BFC4")
p2 <- ggpar(p2, font.legend = c(16), font.main = c(12), font.xtickslab = c(9)) 

p3 <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined[data.combined$subtype=="basal" & data.combined$stage=="III-higher",]),
                 legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", conf.int = FALSE, 
                 pval = TRUE, pval.size = 6, pval.coord = c(1800, 0.9), surv.median.line = "hv", xlim=c(0, max(data.combined$futime)), break.x.by = 250,
                 ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival", 
                 title = "GEP4-high (basal), tumor stage III or higher") 
p3$plot <- p3$plot + 
  ggplot2::annotate("text", x = 200, y = 0.1, 
  label = paste0("n = ", sum(data.combined$subtype=="basal" & data.combined$stage=="III-higher" & data.combined$hypoxia_ISR=="high")), 
  size = 6, color = "#F8766D") + 
  ggplot2::annotate("text", x = 600, y = 0.6, 
  label = paste0("n = ", sum(data.combined$subtype=="basal" & data.combined$stage=="III-higher" & data.combined$hypoxia_ISR=="low")), 
  size = 6, color = "#00BFC4")
p3 <- ggpar(p3, font.legend = c(16), font.main = c(12), font.xtickslab = c(9)) 

p4 <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined[data.combined$subtype=="classical" & data.combined$stage=="I-IIa",]),
                 legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", conf.int = FALSE, 
                 pval = TRUE, pval.size = 6, pval.coord = c(1800, 0.9), surv.median.line = "hv", xlim=c(0, max(data.combined$futime)), break.x.by = 250,
                 ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival", 
                 title = "GEP4-low (classical), tumor stage IIA or lower") 
p4$plot <- p4$plot + 
  ggplot2::annotate("text", x = 300, y = 0.3,
  label = paste0("n = ", sum(data.combined$subtype=="classical" & data.combined$stage=="I-IIa" & data.combined$hypoxia_ISR=="high")),
  size = 6, color = "#F8766D") + 
  ggplot2::annotate("text", x = 1200, y = 0.68, 
  label = paste0("n = ", sum(data.combined$subtype=="classical" & data.combined$stage=="I-IIa" & data.combined$hypoxia_ISR=="low")), 
  size = 6, color = "#00BFC4")
p4 <- ggpar(p4, font.legend = c(16), font.main = c(12), font.xtickslab = c(9)) 

p5 <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined[data.combined$subtype=="classical" & data.combined$stage=="IIb",]),
                 legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", conf.int = FALSE, 
                 pval = TRUE, pval.size = 6, pval.coord = c(1800, 0.9), surv.median.line = "hv", xlim=c(0, max(data.combined$futime)), break.x.by = 250,
                 ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival", 
                 title = "GEP4-low (classical), tumor stage IIB") 
p5$plot <- p5$plot + 
  ggplot2::annotate("text", x = 400, y = 0.3, 
  label = paste0("n = ", sum(data.combined$subtype=="classical" & data.combined$stage=="IIb" & data.combined$hypoxia_ISR=="high")), 
  size = 6, color = "#F8766D") + 
  ggplot2::annotate("text", x = 1000, y = 0.6, 
  label = paste0("n = ", sum(data.combined$subtype=="classical" & data.combined$stage=="IIb" & data.combined$hypoxia_ISR=="low")), 
  size = 6, color = "#00BFC4")
p5 <- ggpar(p5, font.legend = c(16), font.main = c(12), font.xtickslab = c(9)) 

p6 <- ggsurvplot(survfit(Surv(futime, event) ~ hypoxia_ISR, data = data.combined[data.combined$subtype=="classical" & data.combined$stage=="III-higher",]),
                 legend.title = "GEP14", legend.labs = c("high", "low"), linetype = "strata", conf.int = FALSE, 
                 pval = TRUE, pval.size = 6, pval.coord = c(1800, 0.9), surv.median.line = "hv", xlim=c(0, max(data.combined$futime)), break.x.by = 250,
                 ggtheme = theme_bw(base_size = 16), xlab = "Time to Death (Days)", ylab = "Percentage Survival", 
                 title = "GEP4-low (classical), tumor stage III or higher") 
p6$plot <- p6$plot + 
  ggplot2::annotate("text", x = 200, y = 0.3, 
  label = paste0("n = ", sum(data.combined$subtype=="classical" & data.combined$stage=="III-higher" & data.combined$hypoxia_ISR=="high")), 
  size = 6, color = "#F8766D") + 
  ggplot2::annotate("text", x = 800, y = 0.9, 
  label = paste0("n = ", sum(data.combined$subtype=="classical" & data.combined$stage=="III-higher" & data.combined$hypoxia_ISR=="low")),
  size = 6, color = "#00BFC4")
p6 <- ggpar(p6, font.legend = c(16), font.main = c(12), font.xtickslab = c(9)) 

### plot all the panels together
plot_grid(p1$plot, p2$plot, p3$plot, p4$plot, p5$plot, p6$plot, nrow = 2)

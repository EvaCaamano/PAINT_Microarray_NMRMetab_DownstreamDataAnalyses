# Packages (install if missing)
# -------------------------

# load key libraries
require(ggplot2)
require(limma)
require(reshape)
library(variancePartition)
library(ggpubr)
require(statmod)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)
library(tidyr)
library(fgsea)
library(data.table)
library(ReactomePA)
library(data.table)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer)

##### MICROARRAY ANALYSES----
#read in data file for microarray and metadata. Note this study had samples from both non supplemented and Arg-supplemented 
#patients. We will only be focusing on the non-supplemented patients for this study
LTU<-read.csv("MicroarrayData_noMar.csv",header = T,stringsAsFactors = F,row.names = 1)
lauMeta<-read.csv("meta1.csv",header=T,stringsAsFactors = F)

#Discard columns pertaining to resuLTUs from other statistical analyses/hypotheses
LTU_data<-LTU[,11:ncol(LTU)]

#format the columns to extract relevant information
samplenames_LTU<-colnames(LTU_data)
samplenames_LTUF<-strsplit2(samplenames_LTU,"D")
samplenames_LTUF[,2]<-strsplit2(samplenames_LTUF[,2],"A..")[,1]
samplenames_LTUF[,1]<-gsub("X","",samplenames_LTUF[,1])#remove x
samplenames_LTUF[,1]<-gsub("\\.","",samplenames_LTUF[,1])#remove dot
samplenames_LTUF[,2]<-gsub("\\.","",samplenames_LTUF[,2])#remove dot
samplenames_LTUF<-as.data.frame(samplenames_LTUF)
colnames(samplenames_LTUF)<-c("Patient","Time")
samplenames_LTUF$PT<-paste("P",samplenames_LTUF[,1],"_t",samplenames_LTUF[,2],sep="")

#a bit more info on the grouping
samplenames_LTUF$Patient<-as.integer(samplenames_LTUF$Patient)
samplenames_LTUF$Intervention<-lauMeta[match(samplenames_LTUF[,1],lauMeta$Study.No),"Arginine.Supplemented"]
samplenames_LTUF$Sex<-lauMeta[match(samplenames_LTUF[,1],lauMeta$Study.No),"Sex"]

samplenames_LTUF$All<-paste(samplenames_LTUF$PT,samplenames_LTUF$Intervention,samplenames_LTUF$Sex,sep="_")

colnames(LTU_data)<-samplenames_LTUF$All

#then we select only the non supplemented samples and create appropriate factors
LTU_NS<-LTU_data[,strsplit2(colnames(LTU_data),"_")[,3]%in%c("N")]
TimeFactor_NS<-strsplit2(colnames(LTU_NS),"_")[,2]
Sex_NS<-strsplit2(colnames(LTU_NS),"_")[,4]
Patient_NS<-strsplit2(colnames(LTU_NS),"_")[,1]

#Differential Gene analsyis 
set.seed(800)
data_for_limma<-data.frame("Time"=TimeFactor_NS,
                           "Sex"=Sex_NS,
                           "PatientID"=Patient_NS,
                           t(LTU_NS))

data_for_limma<-data_for_limma[complete.cases(data_for_limma),]
data_for_limma$Time<-factor(data_for_limma$Time,levels=c("t3","t10"))
data_for_limma$PatientID<-factor(data_for_limma$PatientID)
data_for_limma$Sex<-factor(data_for_limma$Sex,levels=c("M","F"))
numeric_data<-data_for_limma[,4:ncol(data_for_limma)]

#### Analysis of variance to determine covariates----


form <- ~(1|Time)+(1|Sex)+(1|PatientID)

#note as per package recommendation, all categorical variables are modeled as random effects with continuous as fixed effects

varPart <- fitExtractVarPartModel(t(numeric_data), form, data_for_limma[,1:3])

vp <- sortCols(varPart,last="Residuals")
vp_bars <- plotPercentBars(vp[vp$Sex>0.2 &vp$Time>0.2 &vp$PatientID>0.25, ])+theme(legend.position="bottom")#top most affected by time and sex
vp_barsT <- plotPercentBars(vp[vp$Sex<0.01 &vp$Time>0.55 &vp$PatientID<0.01, ])+theme(legend.position="bottom")#most affected by time
vp_barsS <- plotPercentBars(vp[vp$Sex>0.55 &vp$Time<0.05&vp$PatientID<0.01, ])+theme(legend.position="bottom")#most affected by sex
vp_barsP <- plotPercentBars(vp[vp$Sex<0.01 &vp$Time<0.01&vp$PatientID>0.9992, ])+theme(legend.position="bottom")#most affected by patientID
vp_violin <- plotVarPart(vp) 

ggarrange(vp_bars, vp_barsT,vp_barsS,vp_barsP, vp_violin, 
          labels = c("A", "B","C","D","E"),
          ncol = 2, nrow = 3)

micro_varPar_plot<-ggarrange(vp_bars, vp_barsT,vp_barsS,vp_barsP, vp_violin, 
                             labels = c("A", "B","C","D","E"),
                             ncol = 2, nrow = 3)
ggsave(plot = micro_varPar_plot,filename = "micro_varPar_plot_wPID.tiff", units="in", width=9, height=11, dpi=300, compression = 'lzw',bg="white")


# DIFFERENTIAL GENE EXPRESSION ####
##### limma with duplicated correlations, group only
# DGE using the correlation structure to account for repeated measurements, no taking into account sex

# Group first, intercept removed
designNoPatNoSex<-model.matrix(~0+Time,data=data_for_limma)

colnames(designNoPatNoSex)<-c("t3","t10")#for sex F is 1

colnames(designNoPatNoSex)


# Correlation structure (accounting for repeated measures)
corfit_noPatnoSex <- duplicateCorrelation(t(numeric_data), designNoPatNoSex, block = data_for_limma$PatientID)

fit_corfit_noSex <- lmFit(t(numeric_data), designNoPatNoSex, block = data_for_limma$PatientID, correlation = corfit_noPatnoSex$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrix2 <- makeContrasts(Treatment_vs_Control = t10 - t3,
  levels = designNoPatNoSex
)

fit2_corfit_noSex <- contrasts.fit(fit_corfit_noSex, contrast.matrix2)
fit2_corfit_noSex <- eBayes(fit2_corfit_noSex)
fit2_corfit_noSex_Res<-topTable(fit2_corfit_noSex, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_data),)

dim(fit2_corfit_noSex_Res[fit2_corfit_noSex_Res$adj.P.Val<0.05,])#785 genes

#comparison with SAM results
length(rownames(fit2_corfit_noSex_Res[fit2_corfit_noSex_Res$adj.P.Val<0.05,])[rownames(fit2_corfit_noSex_Res[fit2_corfit_noSex_Res$adj.P.Val<0.05,])%in%strsplit2(samDaysOutAll[samDaysOutAll$q.value...<5,"Gene.Name"],"_")[,1]])#555 overlap



##### limma with duplicated correlations, group and sex
# Group and SEX- with group first, intercept removed
designNoPat<-model.matrix(~0+Time+Sex,data=data_for_limma)

colnames(designNoPat)<-c("t3","t10","Sex")#for sex F is 1

colnames(designNoPat)


# Correlation structure (accounting for repeated measures- No SEX)
corfit <- duplicateCorrelation(t(numeric_data), designNoPat, block = data_for_limma$PatientID)

fit_corfit <- lmFit(t(numeric_data), designNoPat, block = data_for_limma$PatientID, correlation = corfit$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrix <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                 levels = designNoPat
)

fit2_corfit <- contrasts.fit(fit_corfit, contrast.matrix)
fit2_corfit <- eBayes(fit2_corfit)
fit2_corfit_Res<-topTable(fit2_corfit, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_data),)

dim(fit2_corfit_Res[fit2_corfit_Res$adj.P.Val<0.05,])#only 154

#comparison with SAM results
length(rownames(fit2_corfit_Res[fit2_corfit_Res$adj.P.Val<0.05,])[rownames(fit2_corfit_Res[fit2_corfit_Res$adj.P.Val<0.05,])%in%strsplit2(samDaysOutAll[samDaysOutAll$q.value...<5,"Gene.Name"],"_")[,1]])#119 overlap

#comparison sex and no sex
length(intersect(rownames(fit2_corfit_noSex_Res[fit2_corfit_noSex_Res$adj.P.Val<0.05,]),rownames(fit2_corfit_Res[fit2_corfit_Res$adj.P.Val<0.05,])))#149 overlap

write.csv(fit2_corfit_noSex_Res,"UnivariateResults_MicroarrayData_TimePatient.csv")
write.csv(fit2_corfit_Res,"UnivariateResults_MicroarrayData_TimeSexPatient.csv")

####### FUNCTIONAL ENRICHMENT#####


##### on results from model with time only #####
IDs<-mapIds(org.Hs.eg.db,keys=rownames(fit2_corfit_noSex_Res),
             keytype = "SYMBOL",#what we have
             column="ENTREZID",
             multiVals=first)#what to do if multiple mappings-we get the first in thsi case

fit2_corfit_noSex_Res$EntrezID<-IDs

ranked_gene_list_Time <- fit2_corfit_noSex_Res %>%
  # Calculate rank metric
  mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>% 
  #Sort in descending order
  arrange(desc(rank_metric)) %>% 
  # Extract ranked gene list
  pull(rank_metric, EntrezID) 


reactome_pathways_Time <- reactomePathways(names(ranked_gene_list_Time))

set.seed(800)
fgseaRes_Time <- fgsea(pathways = reactome_pathways_Time,
                          stats = ranked_gene_list_Time,
                          minSize = 10,
                          maxSize = 300)

fgseaRes_Time <- fgseaRes_Time[order(fgseaRes_Time$padj), ]
head(fgseaRes_Time)

#add -log10 transformed FDR values
fgseaRes_Time$log10_p.adjust<--log10(fgseaRes_Time$padj)

fgseaRes_Time_filtered <- fgseaRes_Time[fgseaRes_Time$padj<0.05 & (fgseaRes_Time$NES>2.5 | fgseaRes_Time$NES<= -1) ,] %>%
  arrange(abs(NES)) %>%
  slice_head(n = 50) 

# fgseaRes_Time_filtered <- fgseaRes_Time[fgseaRes_Time$padj<0.02 & abs(fgseaRes_Time$NES)>1.2,] %>%
#   arrange(desc(log10_p.adjust)) %>%
#   slice_head(n = 50)

# Create dot plot
gt <- ggplot(fgseaRes_Time_filtered,
              aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = 1, lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised enrichment score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal(base_size = 14)+theme(
    text = element_text(color = "black", face = "bold"),        # all text black
    axis.text = element_text(color = "black"),   # axis tick labels
    axis.title = element_text(color = "black"),  # axis titles
  ) 

gt
ggsave(plot = gt,filename = "GSEA_Reactome_Time.tiff", units="in", width=10, height=11, dpi=300, compression = 'lzw',bg="white")

dim(fgseaRes_Time[fgseaRes_Time$padj<0.05,])
fgseaRes_Time[grep("Metallothionein",fgseaRes_Time$pathway),]#just under the limit of significance (0.09 adjusted)

#format the leading edge so it is not a list
fgseaRes_Time_tidy <- fgseaRes_Time %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

write.csv(fgseaRes_Time_tidy,"Reactome_GSEARes_TimeOnly.csv",row.names = F)
##### on results from model with sex and time #####
test<-mapIds(org.Hs.eg.db,keys=rownames(fit2_corfit_Res),
             keytype = "SYMBOL",#what we have
             column="ENTREZID",
             multiVals=first)#what to do if multiple mappings-we get the first in thsi case

fit2_corfit_Res$EntrezID<-test

ranked_gene_list_TimeSex <- fit2_corfit_Res %>%
  # Calculate rank metric
  mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>% 
  #Sort in descending order
  arrange(desc(rank_metric)) %>% 
  # Extract ranked gene list
  pull(rank_metric, EntrezID) 


reactome_pathways_TimeSex <- reactomePathways(names(ranked_gene_list_TimeSex))

set.seed(800)
fgseaRes_TimeSex <- fgsea(pathways = reactome_pathways_TimeSex,
                          stats = ranked_gene_list_TimeSex,
                          minSize = 10,
                          maxSize = 300)

fgseaRes_TimeSex <- fgseaRes_TimeSex[order(fgseaRes_TimeSex$padj), ]
head(fgseaRes_TimeSex)

#add -log10 transformed FDR values
fgseaRes_TimeSex$log10_p.adjust<--log10(fgseaRes_TimeSex$padj)

fgseaRes_TimeSex_filtered <- fgseaRes_TimeSex[fgseaRes_TimeSex$padj<0.05 & (fgseaRes_TimeSex$NES>2.4 | fgseaRes_TimeSex$NES<= -1) ,] %>%
  arrange(abs(NES)) %>%
  slice_head(n = 50) 

# Create dot plot
gst <- ggplot(fgseaRes_TimeSex_filtered,
              aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = 1, lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised enrichment score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal(base_size = 14)+
  theme(
    text = element_text(color = "black", face = "bold"),        # all text black
    axis.text = element_text(color = "black"),   # axis tick labels
    axis.title = element_text(color = "black"),  # axis titles
    ) 

gst
ggsave(plot = gst,filename = "GSEA_Reactome_TimeSEX.tiff", units="in", width=10, height=11, dpi=300, compression = 'lzw',bg="white")



#format the leading edge so it is not a list
fgseaRes_TimeSex_tidy <- fgseaRes_TimeSex %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

write.csv(fgseaRes_TimeSex_tidy,"Reactome_GSEARes_TimeSex.csv",row.names = F)

gts<-ggarrange(gt, gst, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
ggsave(plot = gts,filename = "GSEA_Reactome_TIMEA_TimeSEXB.tiff", units="in", width=20, height=11, dpi=300, compression = 'lzw',bg="white")


# Other figures ----

#Heatmap

expr_signif_Time<-t(data_for_limma[,colnames(data_for_limma)%in%rownames(fit2_corfit_noSex_Res[fit2_corfit_noSex_Res$adj.P.Val<0.05,])])
expr_signif_TimeSex<-t(data_for_limma[,colnames(data_for_limma)%in%rownames(fit2_corfit_Res[fit2_corfit_Res$adj.P.Val<0.05,])])


ordered_samples_Time <- order(data_for_limma$Time)
expr_signif_Time_ordered <- expr_signif_Time[, ordered_samples_Time]
expr_signif_TimeSex_ordered <- expr_signif_TimeSex[, ordered_samples_Time]
group_factor_ordered <- data_for_limma$Time[ordered_samples_Time]

# Create annotation data frame
annotation_col <- data.frame(Group = group_factor_ordered)
rownames(annotation_col) <- colnames(expr_signif_Time_ordered)



heat_noSex<-pheatmap(expr_signif_Time_ordered,
         scale = "row",                 # optionally scale genes
         clustering_distance_rows = "euclidean",  # distance metric for genes
         clustering_method = "complete",          # clustering method
         cluster_cols=F,
         show_rownames = F,
         show_colnames = F,
         annotation_col = annotation_col, # adds color bar for sample groups
         color = colorRampPalette(c("blue", "white", "red"))(50)
)

ggsave(plot = heat_noSex,filename = "HeatmapSignifGenes_TimeOnly.tiff", units="in", width=8, height=8, dpi=300, compression = 'lzw',bg="white")


heat_sex<-pheatmap(expr_signif_TimeSex_ordered,
         scale = "row",                 # optionally scale genes
         clustering_distance_rows = "euclidean",  # distance metric for genes
         clustering_method = "complete",          # clustering method
         cluster_cols=F,
         show_rownames = F,
         show_colnames = F,
         annotation_col = annotation_col, # adds color bar for sample groups
         color = colorRampPalette(c("blue", "white", "red"))(50)
)

ggsave(plot = heat_sex,filename = "HeatmapSignifGenes_TimeSex.tiff", units="in", width=8, height=8, dpi=300, compression = 'lzw',bg="white")


#PCA plot

plot_pca <- function(data, color_vec = NULL, shape_vec = NULL,
                     color_legend = "Color", shape_legend = "Shape",
                     color_values=rainbow(length(unique(color_vec)))) {
  
  # Check vector lengths
  n <- nrow(data)
  if (!is.null(color_vec) && length(color_vec) != n) {
    stop("color_vec must have the same length as the number of rows in data")
  }
  if (!is.null(shape_vec) && length(shape_vec) != n) {
    stop("shape_vec must have the same length as the number of rows in data")
  }
  
  
  # Perform PCA
  pca_res <- prcomp(data, scale. = T)
  
  # % variance explained
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
  
  # Prepare data frame for plotting
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )
  
  if (!is.null(color_vec)) pca_df$Color <- color_vec
  if (!is.null(shape_vec)) pca_df$Shape <- shape_vec
  
  # Build ggplot mapping
  mapping <- aes(x = PC1, y = PC2)
  if (!is.null(color_vec)) mapping$colour <- pca_df$Color
  if (!is.null(shape_vec)) {
    # Only use shape if categorical
    if (is.factor(shape_vec) || is.character(shape_vec)) {
      mapping$shape <- pca_df$Shape
    } else {
      warning("shape_vec is not categorical. Ignoring shape mapping.")
    }
  }
  
  # Generate plot
  p <- ggplot(pca_df, mapping = mapping) +
    geom_point(size = 3) +
    xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(size = 14),
      legend.title = element_text(size = 12)
    )
  
  # Set legend titles
  if (!is.null(color_vec)) p <- p + labs(color = color_legend)
  if (!is.null(shape_vec) && (is.factor(shape_vec) || is.character(shape_vec))) {
    p <- p + labs(shape = shape_legend)
  }
  # Set color scale if provided
  if (!is.null(color_vec) && !is.null(color_values)) {
    if (is.character(color_values) && length(color_values) == 1) {
      # Treat as palette name
      p <- p + scale_color_brewer(palette = color_values)
    } else if (is.character(color_values) && length(color_values) > 1) {
      # Treat as manual color vector
      p <- p + scale_color_manual(values = color_values)
    }
  }
  
  return(p)
}

  

PCA_signifGenesTimeOnly<-plot_pca(data = t(expr_signif_Time),color_vec = data_for_limma$Time,shape_vec = data_for_limma$Sex,color_legend = "Time",shape_legend = "Sex",color_values = c("#D81B60","#1E88E5"))

PCA_signifGenesTimeSex<-plot_pca(data = t(expr_signif_TimeSex),color_vec = data_for_limma$Time,shape_vec = data_for_limma$Sex,color_legend = "Time",shape_legend = "Sex",color_values = c("#D81B60","#1E88E5"))

ggsave(plot = PCA_signifGenesTimeOnly,filename = "PCA_signifGenesTimeOnly.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")
ggsave(plot = PCA_signifGenesTimeSex,filename = "PCA_signifGenesTimeSex.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")


# METABOLOMICS DATA ANALYSIS -----


MetabNoNas_extracols<-read.csv("MetabNoNas_updated.csv",header=T)
MetabNoNas_matchTranscr_notSupp<-MetabNoNas[MetabNoNas$Sample.ID%in%c(strsplit2(colnames(LT_NS),"_")[,1]),]


#metabolomics PCA using metabNoNAs ####
MetabNoNas_matchTranscr_notSupp<-MetabNoNas[MetabNoNas$Sample.ID%in%c(strsplit2(colnames(LT_NS),"_")[,1]),]

plot_pca(data = log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)]),color_vec = factor(TimeFactor_NS,ordered = T,levels=c("t3","t10")),shape_vec = NULL,color_legend = "Time",shape_legend = "Sex",color_values = c("#D81B60","#1E88E5"))

#metabnonas_matcht_normagain<-PQN(data = MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)])
#plot_pca(data = metabnonas_matcht_normagain,color_vec = factor(TimeFactor_NS,ordered = T,levels=c("t3","t10")),shape_vec = NULL,color_legend = "Time",shape_legend = "Sex",color_values = c("#D81B60","#1E88E5"))

#analysis of variance metabolomics ----
formMETAB <- ~(1|Time)+(1|Sex)+(1|PatientID)

#note as per package recommendation, all categorical variables are modeled as random effects with continuous as fixed effects

varPartMETAB <- fitExtractVarPartModel(t(numeric_dataMETAB), formMETAB, data_for_limmaMETAB[,1:3])

vpMETAB <- sortCols(varPartMETAB,last="Residuals")
vpMETAB_bars <-plotPercentBars(vpMETAB)
vpMETAB_violin <- plotVarPart(vpMETAB) 

ggarrange(vpMETAB_bars, vpMETAB_violin, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

metab_varPar_plot<-ggarrange(vpMETAB_bars, vpMETAB_violin, 
                             labels = c("A", "B"),
                             ncol = 2, nrow = 1)
ggsave(plot = metab_varPar_plot,filename = "metab_varPar_plot_wPID.tiff", units="in", width=11, height=9, dpi=300, compression = 'lzw',bg="white")



# differential abundance metabolomics ----

set.seed(800)
data_for_limmaMETAB<-data.frame("Time"=MetabNoNas_matchTranscr_notSupp$Day,
                                "Sex"=MetabNoNas_matchTranscr_notSupp$Sex,
                                "PatientID"=MetabNoNas_matchTranscr_notSupp$Sample.ID,
                                log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)]))


data_for_limmaMETAB<-data_for_limmaMETAB[complete.cases(data_for_limmaMETAB),]
data_for_limmaMETAB$Time<-factor(data_for_limmaMETAB$Time,levels=c("3","10"))
data_for_limmaMETAB$PatientID<-factor(data_for_limmaMETAB$PatientID)
data_for_limmaMETAB$Sex<-factor(data_for_limmaMETAB$Sex,levels=c("Male","Female"))
numeric_dataMETAB<-log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)])


# Group first, intercept removed
designNoPatNoSexMETAB<-model.matrix(~0+Time,data=data_for_limmaMETAB)

colnames(designNoPatNoSexMETAB)<-c("t3","t10")#for sex F is 1

colnames(designNoPatNoSexMETAB)


# Correlation structure (accounting for repeated measures)
corfit_noPatnoSexMETAB <- duplicateCorrelation(t(numeric_dataMETAB), designNoPatNoSexMETAB, block = data_for_limmaMETAB$PatientID)

fit_corfit_noSexMETAB <- lmFit(t(numeric_dataMETAB), designNoPatNoSexMETAB, block = data_for_limmaMETAB$PatientID, correlation = corfit_noPatnoSexMETAB$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrix2METAB <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                       levels = designNoPatNoSexMETAB
)

fit2_corfit_noSexMETAB <- contrasts.fit(fit_corfit_noSexMETAB, contrast.matrix2METAB)
fit2_corfit_noSexMETAB <- eBayes(fit2_corfit_noSexMETAB)
fit2_corfit_noSexMETAB_Res<-topTable(fit2_corfit_noSexMETAB, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_data),)

dim(fit2_corfit_noSexMETAB_Res[fit2_corfit_noSexMETAB_Res$adj.P.Val<0.05,])#nothing significant, two above 0.2



pca_metab_signif_noSex<-plot_pca(data = numeric_dataMETAB[,colnames(numeric_dataMETAB)%in%rownames(fit2_corfit_noSexMETAB_Res[fit2_corfit_noSexMETAB_Res$P.Value<0.05 ,])],
         color_vec = factor(data_for_limmaMETAB$Time,ordered = T,levels=c("3","10")),
         shape_vec = factor(data_for_limmaMETAB$Sex),
         color_legend = "Time",shape_legend = "sex",color_values = c("#D81B60","#1E88E5"))
ggsave(plot = pca_metab_signif_noSex,filename = "PCA_signifMetabTimeOnly.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")

write.csv(fit2_corfit_noSexMETAB_Res,"DifferentialAbundanceMetabolites_ModelTimePatient.csv")

## Now with Sex in the models

set.seed(800)
data_for_limmaMETAB<-data.frame("Time"=MetabNoNas_matchTranscr_notSupp$Day,
                                "Sex"=MetabNoNas_matchTranscr_notSupp$Sex,
                                "PatientID"=MetabNoNas_matchTranscr_notSupp$Sample.ID,
                                log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)]))


data_for_limmaMETAB<-data_for_limmaMETAB[complete.cases(data_for_limmaMETAB),]
data_for_limmaMETAB$Time<-factor(data_for_limmaMETAB$Time,levels=c("3","10"))
data_for_limmaMETAB$PatientID<-factor(data_for_limmaMETAB$PatientID)
data_for_limmaMETAB$Sex<-factor(data_for_limmaMETAB$Sex,levels=c("Male","Female"))
numeric_dataMETAB<-log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)])

# Group first, intercept removed
designNoPatSexMETAB<-model.matrix(~0+Time+Sex,data=data_for_limmaMETAB)

colnames(designNoPatSexMETAB)<-c("t3","t10","Sex")#for sex F is 1

colnames(designNoPatSexMETAB)


# Correlation structure (accounting for repeated measures)
corfit_noPatSexMETAB <- duplicateCorrelation(t((numeric_dataMETAB)), designNoPatSexMETAB, block = data_for_limmaMETAB$PatientID)

fit_corfit_SexMETAB <- lmFit(t(numeric_dataMETAB), designNoPatSexMETAB, block = data_for_limmaMETAB$PatientID, correlation = corfit_noPatSexMETAB$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrix2METAB_Sex <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                           levels = designNoPatSexMETAB
)

fit2_corfit_SexMETAB <- contrasts.fit(fit_corfit_SexMETAB, contrast.matrix2METAB_Sex)
fit2_corfit_SexMETAB <- eBayes(fit2_corfit_SexMETAB)
fit2_corfit_SexMETAB_Res<-topTable(fit2_corfit_SexMETAB, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_dataMETAB),)

dim(fit2_corfit_SexMETAB_Res[fit2_corfit_SexMETAB_Res$adj.P.Val<0.05,])#o-acetylcarnitine is significant.


#pca plot with those handful of metabolites significant prior correction
pca_metab_signif_Sex<-plot_pca(data = numeric_dataMETAB[,colnames(numeric_dataMETAB)%in%rownames(fit2_corfit_SexMETAB_Res[fit2_corfit_SexMETAB_Res$P.Value<0.05 ,])],
         color_vec = factor(data_for_limmaMETAB$Time,ordered = T,levels=c("3","10")),
         shape_vec = factor(data_for_limmaMETAB$Sex),
         color_legend = "Time",shape_legend = "sex",color_values = c("#D81B60","#1E88E5"))
ggsave(plot = pca_metab_signif_Sex,filename = "PCA_signifMetabTimeSex.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")

write.csv(fit2_corfit_SexMETAB_Res,"DifferentialAbundanceMetabolites_ModelTimeSexPatient.csv")

#------------------------------------------------------
#---- EXTRA ANALYSES TAKING INTO ACCOUNG GESTATION ----
#------------------------------------------------------


set.seed(800)
LTU_NS<-LTU_data[,strsplit2(colnames(LTU_data),"_")[,3]%in%c("N")]
TimeFactor_NS<-strsplit2(colnames(LTU_NS),"_")[,2]
Sex_NS<-strsplit2(colnames(LTU_NS),"_")[,4]
Patient_NS<-strsplit2(colnames(LTU_NS),"_")[,1]
lauMeta2<-lauMeta
lauMeta2$PatientID<-as.character(lauMeta$Study.No)
lauMeta2$PatientID[lauMeta2$Study.No<10]<-paste0("P0",lauMeta$Study.No[lauMeta2$Study.No<10])
lauMeta2$PatientID[lauMeta2$Study.No>10]<-paste0("P",lauMeta$Study.No[lauMeta2$Study.No>10])

Gestational_Age<-lauMeta[match(Patient_NS,lauMeta2$PatientID),"Gest"]


data_for_limmaGEST<-data.frame("Time"=TimeFactor_NS,
                               "Sex"=Sex_NS,
                               "PatientID"=Patient_NS,
                               "Gest"=Gestational_Age,
                               t(LTU_NS))

data_for_limmaGEST<-data_for_limmaGEST[complete.cases(data_for_limmaGEST),]
data_for_limmaGEST$Time<-factor(data_for_limmaGEST$Time,levels=c("t3","t10"))
data_for_limmaGEST$PatientID<-factor(data_for_limmaGEST$PatientID)
data_for_limmaGEST$Sex<-factor(data_for_limmaGEST$Sex,levels=c("M","F"))
numeric_data<-data_for_limmaGEST[,5:ncol(data_for_limmaGEST)]

#### Analysis of variance to determine covariates----


formGEST <- ~(1|Time)+Gest+(1|PatientID)+(1|Sex)



#note as per package recommendation, all categorical variables are modeled as random effects with continuous as fixed effects

varPartGEST <- fitExtractVarPartModel(t(numeric_data), formGEST, data_for_limmaGEST[,1:4])

vpGEST <- sortCols(varPartGEST,last="Residuals")
vpGEST_bars <- plotPercentBars(vpGEST[vpGEST$Sex>0.15 &vpGEST$Time>0.15 &vpGEST$PatientID>0.15 &vpGEST$Gest>0.15, ])+theme(legend.position="bottom")#top most affected by time and sex
vpGEST_barsT <- plotPercentBars(vpGEST[vpGEST$Sex<0.01 &vpGEST$Time>0.48 &vpGEST$PatientID<0.01 &vpGEST$Gest<0.01, ])+theme(legend.position="bottom")#most affected by time
vpGEST_barsS <- plotPercentBars(vpGEST[vpGEST$Sex>0.35 &vpGEST$Time<0.05&vpGEST$PatientID<0.01 &vpGEST$Gest<0.01, ])+theme(legend.position="bottom")#most affected by sex
vpGEST_barsP <- plotPercentBars(vpGEST[vpGEST$Sex<0.01 &vpGEST$Time<0.01&vpGEST$PatientID>0.995 &vpGEST$Gest<0.01, ])+theme(legend.position="bottom")#most affected by patientID
vpGEST_barsG <- plotPercentBars(vpGEST[vpGEST$Sex<0.01 &vpGEST$Time<0.01&vpGEST$PatientID<0.01 &vpGEST$Gest>0.4, ])+theme(legend.position="bottom")
vpGEST_violin <- plotVarPart(vpGEST) 

ggarrange(vpGEST_bars, vpGEST_barsT,vpGEST_barsS,vpGEST_barsP,vpGEST_barsG, vpGEST_violin, 
          labels = c("A", "B","C","D","E", "F"),
          ncol = 2, nrow = 3)

micro_varPar_plotGEST<-ggarrange(vpGEST_bars, vpGEST_barsT,vpGEST_barsS,vpGEST_barsP,vpGEST_barsG, vpGEST_violin, 
                                 labels = c("A", "B","C","D","E", "F"),
                                 ncol = 2, nrow = 3)
ggsave(plot = micro_varPar_plotGEST,filename = "micro_varPar_plot_wPIDGEST.tiff", units="in", width=10, height=11.5, dpi=300, compression = 'lzw',bg="white")


# DIFFERENTIAL GENE EXPRESSION ####
##### limma with duplicated correlations, group only
# DGE using the correlation structure to account for repeated measurements, no taking into account sex

# Group first, intercept removed
designNoPatYesGestNoSexYesGest<-model.matrix(~0+Time+Gest,data=data_for_limmaGEST)

colnames(designNoPatYesGestNoSexYesGest)<-c("t3","t10","Gest")#for sex F is 1

colnames(designNoPatYesGestNoSexYesGest)


# Correlation structure (accounting for repeated measures)
corfitYesGest_noPatnoSexYesGest <- duplicateCorrelation(t(numeric_data), designNoPatYesGestNoSexYesGest, block = data_for_limmaGEST$PatientID)

fit_corfitYesGest_noSexYesGest <- lmFit(t(numeric_data), designNoPatYesGestNoSexYesGest, block = data_for_limmaGEST$PatientID, correlation = corfitYesGest_noPatnoSexYesGest$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrixYesGest2YesGest <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                                levels = designNoPatYesGestNoSexYesGest
)

fit2_corfitYesGest_noSexYesGest <- contrasts.fit(fit_corfitYesGest_noSexYesGest, contrast.matrixYesGest2YesGest)
fit2_corfitYesGest_noSexYesGest <- eBayes(fit2_corfitYesGest_noSexYesGest)
fit2_corfitYesGest_noSexYesGest_Res<-topTable(fit2_corfitYesGest_noSexYesGest, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_data),)

dim(fit2_corfitYesGest_noSexYesGest_Res[fit2_corfitYesGest_noSexYesGest_Res$adj.P.Val<0.05,])#836 genes

#comparison with SAM results and previous limma results
length(rownames(fit2_corfitYesGest_noSexYesGest_Res[fit2_corfitYesGest_noSexYesGest_Res$adj.P.Val<0.05,])[rownames(fit2_corfitYesGest_noSexYesGest_Res[fit2_corfitYesGest_noSexYesGest_Res$adj.P.Val<0.05,])%in%strsplit2(samDaysOutAll[samDaysOutAll$q.value...<5,"Gene.Name"],"_")[,1]])#582 overlap
length(rownames(fit2_corfitYesGest_noSexYesGest_Res[fit2_corfitYesGest_noSexYesGest_Res$adj.P.Val<0.05,])[rownames(fit2_corfit_noSexYesGest_Res[fit2_corfitYesGest_noSexYesGest_Res$adj.P.Val<0.05,])%in%rownames(fit2_corfit_noSex_Res[fit2_corfit_noSex_Res$adj.P.Val<0.05,])])#733/836 overlapping



##### limma with duplicated correlations, group and sex
# Group and SEX- with group first, intercept removed
designNoPatYesGest<-model.matrix(~0+Time+Gest+Sex,data=data_for_limmaGEST)

colnames(designNoPatYesGest)<-c("t3","t10","Gest","Sex")#for sex F is 1

colnames(designNoPatYesGest)


# Correlation structure (accounting for repeated measures- No SEX)
corfitYesGest <- duplicateCorrelation(t(numeric_data), designNoPatYesGest, block = data_for_limmaGEST$PatientID)

fit_corfitYesGest <- lmFit(t(numeric_data), designNoPatYesGest, block = data_for_limmaGEST$PatientID, correlation = corfitYesGest$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrixYesGest <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                        levels = designNoPatYesGest
)

fit2_corfitYesGest <- contrasts.fit(fit_corfitYesGest, contrast.matrixYesGest)
fit2_corfitYesGest <- eBayes(fit2_corfitYesGest)
fit2_corfitYesGest_Res<-topTable(fit2_corfitYesGest, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_data),)

dim(fit2_corfitYesGest_Res[fit2_corfitYesGest_Res$adj.P.Val<0.05,])#only 98

#comparison with SAM results
length(rownames(fit2_corfitYesGest_Res[fit2_corfitYesGest_Res$adj.P.Val<0.05,])[rownames(fit2_corfitYesGest_Res[fit2_corfitYesGest_Res$adj.P.Val<0.05,])%in%strsplit2(samDaysOutAll[samDaysOutAll$q.value...<5,"Gene.Name"],"_")[,1]])#74 overlap

#comparison sex and no sex
length(intersect(rownames(fit2_corfit_noSexYesGest_Res[fit2_corfit_noSexYesGest_Res$adj.P.Val<0.05,]),rownames(fit2_corfitYesGest_Res[fit2_corfitYesGest_Res$adj.P.Val<0.05,])))#96/98 overlap

write.csv(fit2_corfit_noSexYesGest_Res,"UnivariateResults_MicroarrayData_TimePatient_Gestation.csv")
write.csv(fit2_corfitYesGest_Res,"UnivariateResults_MicroarrayData_TimeSexPatient_Gestation.csv")

####### FUNCTIONAL ENRICHMENT#####


##### on results from model with time only #####
IDs_YesGest<-mapIds(org.Hs.eg.db,keys=rownames(fit2_corfit_noSexYesGest_Res),
                    keytype = "SYMBOL",#what we have
                    column="ENTREZID",
                    multiVals=first)#what to do if multiple mappings-we get the first in thsi case

fit2_corfit_noSexYesGest_Res$EntrezID<-IDs_YesGest

ranked_gene_list_Time_YesGest <- fit2_corfit_noSexYesGest_Res %>%
  # Calculate rank metric
  mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>% 
  #Sort in descending order
  arrange(desc(rank_metric)) %>% 
  # Extract ranked gene list
  pull(rank_metric, EntrezID) 


reactome_pathways_Time_YesGest <- reactomePathways(names(ranked_gene_list_Time_YesGest))

set.seed(800)
fgseaRes_Time_YesGest <- fgsea(pathways = reactome_pathways_Time_YesGest,
                               stats = ranked_gene_list_Time_YesGest,
                               minSize = 10,
                               maxSize = 300)

fgseaRes_Time_YesGest <- fgseaRes_Time_YesGest[order(fgseaRes_Time_YesGest$padj), ]
head(fgseaRes_Time_YesGest)

#add -log10 transformed FDR values
fgseaRes_Time_YesGest$log10_p.adjust<--log10(fgseaRes_Time_YesGest$padj)

fgseaRes_Time_YesGest_filtered <- fgseaRes_Time_YesGest[fgseaRes_Time_YesGest$padj<0.05 & (fgseaRes_Time_YesGest$NES>2.5 | fgseaRes_Time_YesGest$NES<= -1) ,] %>%
  arrange(abs(NES)) %>%
  slice_head(n = 50) 

# fgseaRes_Time_YesGest_filtered <- fgseaRes_Time_YesGest[fgseaRes_Time_YesGest$padj<0.02 & abs(fgseaRes_Time_YesGest$NES)>1.2,] %>%
#   arrange(desc(log10_p.adjust)) %>%
#   slice_head(n = 50)

# Create dot plot
gt_YesGest <- ggplot(fgseaRes_Time_YesGest_filtered,
                     aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = 1, lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised enrichment score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal(base_size = 14)+theme(
    text = element_text(color = "black", face = "bold"),        # all text black
    axis.text = element_text(color = "black"),   # axis tick labels
    axis.title = element_text(color = "black"),  # axis titles
  ) 

gt_YesGest
ggsave(plot = gt_YesGest,filename = "GSEA_Reactome_Time_YesGest.tiff", units="in", width=11, height=11.5, dpi=300, compression = 'lzw',bg="white")

dim(fgseaRes_Time_YesGest[fgseaRes_Time_YesGest$padj<0.05,])
fgseaRes_Time_YesGest[grep("Metallothionein",fgseaRes_Time_YesGest$pathway),]#just under the limit of significance (0.09 adjusted)

#format the leading edge so it is not a list
fgseaRes_Time_YesGest_tidy <- fgseaRes_Time_YesGest %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

write.csv(fgseaRes_Time_YesGest_tidy,"Reactome_GSEARes_TimeANDGESTATION.csv",row.names = F)
##### on results from model with sex and time #####
test_YesGest<-mapIds(org.Hs.eg.db,keys=rownames(fit2_corfitYesGest_Res),
                     keytype = "SYMBOL",#what we have
                     column="ENTREZID",
                     multiVals=first)#what to do if multiple mappings-we get the first in thsi case

fit2_corfitYesGest_Res$EntrezID<-test_YesGest

ranked_gene_list_Time_YesGestSex <- fit2_corfitYesGest_Res %>%
  # Calculate rank metric
  mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>% 
  #Sort in descending order
  arrange(desc(rank_metric)) %>% 
  # Extract ranked gene list
  pull(rank_metric, EntrezID) 


reactome_pathways_Time_YesGestSex <- reactomePathways(names(ranked_gene_list_Time_YesGestSex))

set.seed(800)
fgseaRes_Time_YesGestSex <- fgsea(pathways = reactome_pathways_Time_YesGestSex,
                                  stats = ranked_gene_list_Time_YesGestSex,
                                  minSize = 10,
                                  maxSize = 300)

fgseaRes_Time_YesGestSex <- fgseaRes_Time_YesGestSex[order(fgseaRes_Time_YesGestSex$padj), ]
head(fgseaRes_Time_YesGestSex)

#add -log10 transformed FDR values
fgseaRes_Time_YesGestSex$log10_p.adjust<--log10(fgseaRes_Time_YesGestSex$padj)

fgseaRes_Time_YesGestSex_filtered <- fgseaRes_Time_YesGestSex[fgseaRes_Time_YesGestSex$padj<0.05 & (fgseaRes_Time_YesGestSex$NES>2.4 | fgseaRes_Time_YesGestSex$NES<= -1) ,] %>%
  arrange(abs(NES)) %>%
  slice_head(n = 50) 

# Create dot plot
gst_YesGest <- ggplot(fgseaRes_Time_YesGestSex_filtered,
                      aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = 1, lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised enrichment score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal(base_size = 14)+
  theme(
    text = element_text(color = "black", face = "bold"),        # all text black
    axis.text = element_text(color = "black"),   # axis tick labels
    axis.title = element_text(color = "black"),  # axis titles
  ) 

gst_YesGest
ggsave(plot = gst_YesGest,filename = "GSEA_Reactome_TimeSEXandGESTATION.tiff", units="in", width=11, height=11.5, dpi=300, compression = 'lzw',bg="white")



#format the leading edge so it is not a list
fgseaRes_Time_YesGestSex_tidy <- fgseaRes_Time_YesGestSex %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

write.csv(fgseaRes_Time_YesGestSex_tidy,"Reactome_GSEARes_TimeSexAndGestation.csv",row.names = F)

gt_YesGests<-ggarrange(gt_YesGest, gst_YesGest, 
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1)
ggsave(plot = gt_YesGests,filename = "GSEA_Reactome_TIMEA_TimeSEXB_andGESTATION.tiff", units="in", width=21, height=11, dpi=300, compression = 'lzw',bg="white")


# Other figures ----

#Heatmap

expr_signif_Time_YesGestSex<-t(data_for_limmaGEST[,colnames(data_for_limmaGEST)%in%rownames(fit2_corfit_noSexYesGest_Res[fit2_corfit_noSexYesGest_Res$adj.P.Val<0.05,])])
expr_signif_Time_YesGestSexSex<-t(data_for_limmaGEST[,colnames(data_for_limmaGEST)%in%rownames(fit2_corfitYesGest_Res[fit2_corfitYesGest_Res$adj.P.Val<0.05,])])


ordered_samples_TimeGEST <- order(data_for_limmaGEST$Time)
expr_signif_Time_YesGestSex_ordered <- expr_signif_Time_YesGestSex[, ordered_samples_TimeGEST]
expr_signif_Time_YesGestSexSex_ordered <- expr_signif_Time_YesGestSexSex[, ordered_samples_TimeGEST]
group_factor_orderedGEST <- data_for_limmaGEST$Time[ordered_samples_TimeGEST]

# Create annotation data frame
annotation_colGEST <- data.frame(Group = group_factor_orderedGEST)
rownames(annotation_colGEST) <- colnames(expr_signif_Time_YesGestSex_ordered)



heat_noSex_YesGEST<-pheatmap(expr_signif_Time_YesGestSex_ordered,
                             scale = "row",                 # optionally scale genes
                             clustering_distance_rows = "euclidean",  # distance metric for genes
                             clustering_method = "complete",          # clustering method
                             cluster_cols=F,
                             show_rownames = F,
                             show_colnames = F,
                             annotation_colGEST = annotation_colGEST, # adds color bar for sample groups
                             color = colorRampPalette(c("blue", "white", "red"))(50)
)

ggsave(plot = heat_noSex_YesGEST,filename = "HeatmapSignifGenes_TimeAndGestation.tiff", units="in", width=8, height=8, dpi=300, compression = 'lzw',bg="white")


heat_sex_YesGEST<-pheatmap(expr_signif_Time_YesGestSexSex_ordered,
                           scale = "row",                 # optionally scale genes
                           clustering_distance_rows = "euclidean",  # distance metric for genes
                           clustering_method = "complete",          # clustering method
                           cluster_cols=F,
                           show_rownames = F,
                           show_colnames = F,
                           annotation_colGEST = annotation_colGEST, # adds color bar for sample groups
                           color = colorRampPalette(c("blue", "white", "red"))(50)
)

ggsave(plot = heat_sex_YesGEST,filename = "HeatmapSignifGenes_TimeSexGest.tiff", units="in", width=8, height=8, dpi=300, compression = 'lzw',bg="white")


#PCA plot


PCA_signifGenesTimeOnly_YesGEST<-plot_pca(data = t(expr_signif_Time_YesGestSex),color_vec = data_for_limmaGEST$Time,shape_vec = data_for_limmaGEST$Sex,color_legend = "Time",shape_legend = "Sex",color_values = c("#D81B60","#1E88E5"))

PCA_signifGenesTimeSexGEST<-plot_pca(data = t(expr_signif_Time_YesGestSexSex),color_vec = data_for_limmaGEST$Time,shape_vec = data_for_limmaGEST$Sex,color_legend = "Time",shape_legend = "Sex",color_values = c("#D81B60","#1E88E5"))

ggsave(plot = PCA_signifGenesTimeOnly_YesGEST,filename = "PCA_signifGenesTimeOnly_YesGEST.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")
ggsave(plot = PCA_signifGenesTimeSexGEST,filename = "PCA_signifGenesTimeSexGEST.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")

# METAB_YesGestOLOMICS DATA ANALYSIS -----


MetabNoNas_extracols<-read.csv("MetabNoNas_updated.csv",header=T)
MetabNoNas_matchTranscr_notSupp<-MetabNoNas[MetabNoNas$Sample.ID%in%c(strsplit2(colnames(LT_NS),"_")[,1]),]


#metabolomics PCA using metabNoNAs ####
MetabNoNas_matchTranscr_notSupp<-MetabNoNas[MetabNoNas$Sample.ID%in%c(strsplit2(colnames(LT_NS),"_")[,1]),]

#analysis of variance metabolomics ----
formMETAB_YesGest <- ~(1|Time)+Gestation+(1|Sex)+(1|PatientID)

#note as per package recommendation, all categorical variables are modeled as random effects with continuous as fixed effects

varPartMETAB_YesGest <- fitExtractVarPartModel(t(numeric_dataMETAB_YesGest), formMETAB_YesGest, data_for_limmaMETAB_YesGest[,1:4])

vpMETAB_YesGest <- sortCols(varPartMETAB_YesGest,last="Residuals")
vpMETAB_YesGest_bars <-plotPercentBars(vpMETAB_YesGest)
vpMETAB_YesGest_violin <- plotVarPart(vpMETAB_YesGest) 

ggarrange(vpMETAB_YesGest_bars, vpMETAB_YesGest_violin, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

metab_varPar_plot_gest<-ggarrange(vpMETAB_YesGest_bars, vpMETAB_YesGest_violin, 
                                  labels = c("A", "B"),
                                  ncol = 2, nrow = 1)
ggsave(plot = metab_varPar_plot_gest,filename = "metab_varPar_plot_wPID_Gest.tiff", units="in", width=11, height=9, dpi=300, compression = 'lzw',bg="white")



# differential abundance metabolomics ----

set.seed(800)
data_for_limmaMETAB_YesGest<-data.frame("Time"=MetabNoNas_matchTranscr_notSupp$Day,
                                        "Sex"=MetabNoNas_matchTranscr_notSupp$Sex,
                                        "PatientID"=MetabNoNas_matchTranscr_notSupp$Sample.ID,
                                        "Gestation"=MetabNoNas_matchTranscr_notSupp$Gestation,
                                        log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)]))


data_for_limmaMETAB_YesGest<-data_for_limmaMETAB_YesGest[complete.cases(data_for_limmaMETAB_YesGest),]
data_for_limmaMETAB_YesGest$Time<-factor(data_for_limmaMETAB_YesGest$Time,levels=c("3","10"))
data_for_limmaMETAB_YesGest$PatientID<-factor(data_for_limmaMETAB_YesGest$PatientID)
data_for_limmaMETAB_YesGest$Sex<-factor(data_for_limmaMETAB_YesGest$Sex,levels=c("Male","Female"))
numeric_dataMETAB_YesGest<-log2(MetabNoNas_matchTranscr_notSupp[,9:ncol(MetabNoNas_matchTranscr_notSupp)])


# Group first, intercept removed
designNoPatNoSexMETAB_YesGest<-model.matrix(~0+Time+Gestation,data=data_for_limmaMETAB_YesGest)

colnames(designNoPatNoSexMETAB_YesGest)<-c("t3","t10","Gest")#for sex F is 1

colnames(designNoPatNoSexMETAB_YesGest)


# Correlation structure (accounting for repeated measures)
corfit_noPatnoSexMETAB_YesGest <- duplicateCorrelation(t(numeric_dataMETAB_YesGest), designNoPatNoSexMETAB_YesGest, block = data_for_limmaMETAB_YesGest$PatientID)

fit_corfit_noSexMETAB_YesGest <- lmFit(t(numeric_dataMETAB_YesGest), designNoPatNoSexMETAB_YesGest, block = data_for_limmaMETAB_YesGest$PatientID, correlation = corfit_noPatnoSexMETAB_YesGest$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrix2METAB_YesGest <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                               levels = designNoPatNoSexMETAB_YesGest
)

fit2_corfit_noSexMETAB_YesGest <- contrasts.fit(fit_corfit_noSexMETAB_YesGest, contrast.matrix2METAB_YesGest)
fit2_corfit_noSexMETAB_YesGest <- eBayes(fit2_corfit_noSexMETAB_YesGest)
fit2_corfit_noSexMETAB_YesGest_Res<-topTable(fit2_corfit_noSexMETAB_YesGest, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_dataMETAB_YesGest),)

dim(fit2_corfit_noSexMETAB_YesGest_Res[fit2_corfit_noSexMETAB_YesGest_Res$adj.P.Val<0.05,])#nothing significant, two above 0.2



pca_metab_signif_noSex_YesGest<-plot_pca(data = numeric_dataMETAB_YesGest[,colnames(numeric_dataMETAB_YesGest)%in%rownames(fit2_corfit_noSexMETAB_YesGest_Res[fit2_corfit_noSexMETAB_YesGest_Res$P.Value<0.05 ,])],
                                         color_vec = factor(data_for_limmaMETAB_YesGest$Time,ordered = T,levels=c("3","10")),
                                         shape_vec = factor(data_for_limmaMETAB_YesGest$Sex),
                                         color_legend = "Time",shape_legend = "sex",color_values = c("#D81B60","#1E88E5"))
ggsave(plot = pca_metab_signif_noSex_YesGest,filename = "PCA_signifMetabTimeOnlyGEST.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")

write.csv(fit2_corfit_noSexMETAB_YesGest_Res,"DifferentialAbundanceMetabolites_ModelTimePatient_YesGestation.csv")

## Now with Sex in the models


designNoPatSexMETAB_YesGest<-model.matrix(~0+Time+Gestation+Sex,data=data_for_limmaMETAB_YesGest)

colnames(designNoPatSexMETAB_YesGest)<-c("t3","t10","Gest","Sex")#for sex F is 1

colnames(designNoPatSexMETAB_YesGest)


# Correlation structure (accounting for repeated measures)
corfit_noPatSexMETAB_YesGest <- duplicateCorrelation(t((numeric_dataMETAB_YesGest)), designNoPatSexMETAB_YesGest, block = data_for_limmaMETAB_YesGest$PatientID)

fit_corfit_SexMETAB_YesGest <- lmFit(t(numeric_dataMETAB_YesGest), designNoPatSexMETAB_YesGest, block = data_for_limmaMETAB_YesGest$PatientID, correlation = corfit_noPatSexMETAB_YesGest$consensus)

# Contrast of interest: Treatment vs Control
contrast.matrix2METAB_YesGest_Sex <- makeContrasts(Treatment_vs_Control = t10 - t3,
                                                   levels = designNoPatSexMETAB_YesGest
)

fit2_corfit_SexMETAB_YesGest <- contrasts.fit(fit_corfit_SexMETAB_YesGest, contrast.matrix2METAB_YesGest_Sex)
fit2_corfit_SexMETAB_YesGest <- eBayes(fit2_corfit_SexMETAB_YesGest)
fit2_corfit_SexMETAB_YesGest_Res<-topTable(fit2_corfit_SexMETAB_YesGest, coef="Treatment_vs_Control",adjust.method = "BH",number = ncol(numeric_dataMETAB_YesGest),)

dim(fit2_corfit_SexMETAB_YesGest_Res[fit2_corfit_SexMETAB_YesGest_Res$adj.P.Val<0.05,])#o-acetylcarnitine is significant.


#pca plot with those handful of metabolites significant prior correction
pca_metab_signif_SexGEST<-plot_pca(data = numeric_dataMETAB_YesGest[,colnames(numeric_dataMETAB_YesGest)%in%rownames(fit2_corfit_SexMETAB_YesGest_Res[fit2_corfit_SexMETAB_YesGest_Res$P.Value<0.05 ,])],
                                   color_vec = factor(data_for_limmaMETAB_YesGest$Time,ordered = T,levels=c("3","10")),
                                   shape_vec = factor(data_for_limmaMETAB_YesGest$Sex),
                                   color_legend = "Time",shape_legend = "sex",color_values = c("#D81B60","#1E88E5"))
ggsave(plot = pca_metab_signif_SexGEST,filename = "PCA_signifMetabTimeSexGEST.tiff", units="in", width=7, height=6, dpi=300, compression = 'lzw',bg="white")

write.csv(fit2_corfit_SexMETAB_YesGest_Res,"DifferentialAbundanceMetabolites_ModelTimeSexPatientandGEST.csv")




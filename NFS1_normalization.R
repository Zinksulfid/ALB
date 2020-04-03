#neue normalisierungsfunktion

normalize_vsn_anna <- function(se, spikeins) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))  # Variance stabilization transformation on assay data
  se_vsn <- se
  spfit = vsn::vsnMatrix(2 ^ assay(se_vsn)[spikeins,], minDataPointsPerStratum = 1)
  assay(se_vsn) <- vsn::predict(spfit, 2 ^ assay(se_vsn))
  return(se_vsn)
}

#DEP script zum normalisieren

library("DEP")
library("vsn")
library("ggplot2")
library("dplyr")
Interactome <- read.csv("Fly_Interactome_curated.csv", header=TRUE, sep=",")
expdesign <- read.table("NFS1_expdesign.txt", header = TRUE, sep="\t")
proteinGroups <- read.table("proteinGroups.txt", header = TRUE, sep="\t")

data <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+")
dim(data)
colnames(data)
data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any()

LFQ_columns <- grep("LFQ.", colnames(data_unique))# get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
#fehlerbehebung
LFQ_columns
dim(data_unique)
data_unique[102:117]

plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_filt)
plot_coverage(data_filt)


#data_norm <- normalize_vsn_anna(se=data_filt, spikeins=c("CG30022", "Jafrac1", "RpS19a", "Hsp60", "alpha-Spec" )) #ohne wt
#data_norm <- normalize_vsn_anna(data_filt, c( "skap", "Tina-1", "EG:87B1.3", "Q8MT18", "Q7PLE6")) #mit wt
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

plot_missval(data_filt)
plot_detect(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

data_diff <- test_diff(data_imp, type = "control", control = "wt")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

plot_pca(dep, x = 1, y = 2, n = 374, point_size = 4)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
plot_volcano(dep, contrast = "F_vs_wt", label_size = 2, add_names = TRUE)

#regression curves

data <- assay(data_filt)
df <- data.frame(data)
NFS1 <- c("CG12264")
name.to.keep <- c("CG30022", "Jafrac1", "RpS19a", "Hsp60", "alpha-Spec")
df_final <- subset(df, row.names(df) %in% name.to.keep) 
df_NFS1 <- subset(df, row.names(df) %in% NFS1)
df_NFS1 <- gather (df_NFS1, "A", "B")
df_plot <- gather(df_final, "A", "B" )
df_plot_f <- full_join(df_plot, df_NFS1, by = "A")# by= "A")
df_plot_f <- cbind(df_plot_f, rep(name.to.keep,16))
colnames (df_plot_f) <- c("sample", "value", "NFS1", "protein")
df_plot_f %>%
  mutate(sample = str_replace(sample, "_.", "")) -> df_plot_N
#cairo.pdf("regression.pdf")
ggplot(df_plot_N, aes(x =value , y = NFS1)) +
  geom_point()+
  geom_smooth(method='lm', formula= y ~ x)
  facet_grid( sample ~ protein, scale="free") 
#dev.off()

# neue normalisierung unter Berücksichtigung der Steigung 
library("dplyr")
library("ggplot2")
library("DEP")
library("tidyr")
        
proteinGroups <- read.table("proteinGroups.txt", header = TRUE, sep="\t")
data <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+")

data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any()
LFQ <- grep("LFQ.", colnames(data_unique))
data_LFQ <- cbind(data_unique["Majority.protein.IDs"], data_unique[LFQ])
data_LFQ[data_LFQ == 0] <- NA
d_LFQ <- drop_na(data_LFQ)
R_values <- data.frame(1,2)
S_values <- data.frame(1,2)

# funktion 1 -> sort for R value
housekeeping <- function (zeile) {
  d_LFQ %>% 
    select(matches(".wt.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("protein", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_wt
  colnames(data_wt) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".FLAG.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_FLAG
  colnames(data_FLAG) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".K.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_K
  colnames(data_K) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches("F.$") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4)%>%
    spread("Majority.protein.IDs", "CG12264") -> data_F
  colnames(data_F) <- c("LFQ", "NFS1", "Protein")
  
  
  model_wt <- lm(NFS1 ~ Protein, data= data_wt)
  R_wt<-summary(model_wt)$r.squared
  model_FLAG <- lm(NFS1 ~ Protein, data_FLAG)
  R_FLAG<-summary(model_FLAG)$r.squared
  model_K <- lm(NFS1 ~ Protein, data_K)
  R_K<-summary(model_K)$r.squared
  model_F <- lm(NFS1 ~ Protein, data_F)
  R_F<-summary(model_F)$r.squared
  
  vec <- c( R_wt, R_FLAG, R_K, R_F)
  mean_R <- mean(vec)
  
  
  return(mean_R)
}

for(i in c(1:dim(d_LFQ)[1])){
  #i <- 1
  d<-d_LFQ[i,]
  zeile <- d[,"Majority.protein.IDs"]#Spalte undefiniert; hier noch data.frame, würde funktionieren, wenn named vector
  if(zeile == "Q9VKD3"){
    R_values <- rbind(R_values, c(NA,NA))
  }
  else {
    #print(zeile)
    mean_R <- housekeeping(zeile) # Muss in Variable gespeichert werden
    #print(c(as.character(zeile), mean_R))
    R_values <- rbind(R_values, c(as.character(zeile), mean_R)) 
  }
}

colnames(R_values) <- c("zeilen", "R")

#funktion 2 -> kontrolle der Steigung
steigung <- function (zeile) {
  d_LFQ %>% 
    select(matches(".wt.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("protein", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_wt
  colnames(data_wt) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".FLAG.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_FLAG
  colnames(data_FLAG) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".K.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_K
  colnames(data_K) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches("F.$") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4)%>%
    spread("Majority.protein.IDs", "CG12264") -> data_F
  colnames(data_F) <- c("LFQ", "NFS1", "Protein")
  
  
  model_wt <- lm(NFS1 ~ Protein, data= data_wt)
  S_wt<-summary(model_wt)$coefficients[2] 
  model_FLAG <- lm(NFS1 ~ Protein, data_FLAG)
  S_FLAG<-summary(model_FLAG)$coefficients[2] 
  model_K <- lm(NFS1 ~ Protein, data_K)
  S_K<-summary(model_K)$coefficients[2] 
  model_F <- lm(NFS1 ~ Protein, data_F)
  S_F<-summary(model_F)$coefficients[2] 
  
  if(S_F > 0 & S_wt >0 & S_K > 0 & S_FLAG >0) {
    S=TRUE
  }
  else{
    S=FALSE
  }
  return(S)
}

for(i in c(1:dim(d_LFQ)[1])){
  #i <- 1
  d<-d_LFQ[i,]
  zeile <- d[,"Majority.protein.IDs"]#Spalte undefiniert; hier noch data.frame, würde funktionieren, wenn named vector
  if(zeile == "Q9VKD3"){
    S_values <- rbind(S_values, c(NA,NA))
  }
  else {
    #print(zeile)
    S <- steigung(zeile)
    # Muss in Variable gespeichert werden
    #print(c(as.character(zeile), mean_R))
    S_values <- rbind(S_values, c(as.character(zeile), S)) 
  }
}
colnames(S_values) <- c("zeilen", "S")
data_S_R <- full_join (R_values, S_values, by = "zeilen")
data_sub <- subset (data_S_R, S_values == TRUE)
data_final <- arrange(data_sub, desc(R))
data_final[1:5,]

#extension of housekeeping

library("dplyr")
library("ggplot2")
library("DEP")
library("tidyr")
        
proteinGroups <- read.table("proteinGroups.txt", header = TRUE, sep="\t")
data <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+")

data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any()
LFQ <- grep("LFQ.", colnames(data_unique))
data_LFQ <- cbind(data_unique["Majority.protein.IDs"], data_unique[LFQ])
data_LFQ[data_LFQ == 0] <- NA
d_LFQ <- drop_na(data_LFQ)
R_values <- data.frame(1,2)
S_values <- data.frame(1,2)
# funktion 1 -> sort for R value
housekeeping <- function (zeile) {
  d_LFQ %>% 
    select(matches(".wt.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("protein", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_wt
  colnames(data_wt) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".FLAG.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_FLAG
  colnames(data_FLAG) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".K.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_K
  colnames(data_K) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches("F.$") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4)%>%
    spread("Majority.protein.IDs", "CG12264") -> data_F
  colnames(data_F) <- c("LFQ", "NFS1", "Protein")
  
  
  model_wt <- lm(NFS1 ~ Protein, data= data_wt)
  R_wt<-summary(model_wt)$r.squared
  model_FLAG <- lm(NFS1 ~ Protein, data_FLAG)
  R_FLAG<-summary(model_FLAG)$r.squared
  model_K <- lm(NFS1 ~ Protein, data_K)
  R_K<-summary(model_K)$r.squared
  model_F <- lm(NFS1 ~ Protein, data_F)
  R_F<-summary(model_F)$r.squared
  
  vec <- c( R_wt, R_FLAG, R_K, R_F)
  mean_R <- mean(vec)
  
  
  return(mean_R)
}

for(i in c(1:dim(d_LFQ)[1])){
  #i <- 1
  d<-d_LFQ[i,]
  zeile <- d[,"Majority.protein.IDs"]#Spalte undefiniert; hier noch data.frame, würde funktionieren, wenn named vector
  if(zeile == "Q9VKD3"){
    R_values <- rbind(R_values, c(NA,NA))
  }
  else {
    #print(zeile)
    mean_R <- housekeeping(zeile) # Muss in Variable gespeichert werden
    #print(c(as.character(zeile), mean_R))
    R_values <- rbind(R_values, c(as.character(zeile), mean_R)) 
  }
}

colnames(R_values) <- c("zeilen", "R")

#funktion 2 -> kontrolle der Steigung
steigung <- function (zeile) {
  d_LFQ %>% 
    select(matches(".wt.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("protein", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_wt
  colnames(data_wt) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".FLAG.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_FLAG
  colnames(data_FLAG) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".K.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_K
  colnames(data_K) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches("F.$") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==zeile | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Muc11A", "CG12264", 1:4)%>%
    spread("Majority.protein.IDs", "CG12264") -> data_F
  colnames(data_F) <- c("LFQ", "NFS1", "Protein")
  
  
  model_wt <- lm(NFS1 ~ Protein, data= data_wt)
  S_wt<-summary(model_wt)$coefficients[2] 
  model_FLAG <- lm(NFS1 ~ Protein, data_FLAG)
  S_FLAG<-summary(model_FLAG)$coefficients[2] 
  model_K <- lm(NFS1 ~ Protein, data_K)
  S_K<-summary(model_K)$coefficients[2] 
  model_F <- lm(NFS1 ~ Protein, data_F)
  S_F<-summary(model_F)$coefficients[2] 
  
  if(S_F > 0 & S_wt >0 & S_K > 0 & S_FLAG >0) {
    S=TRUE
  }
  else{
    S=FALSE
  }
  return(S)
}

for(i in c(1:dim(d_LFQ)[1])){
  #i <- 1
  d<-d_LFQ[i,]
  zeile <- d[,"Majority.protein.IDs"]#Spalte undefiniert; hier noch data.frame, würde funktionieren, wenn named vector
  if(zeile == "Q9VKD3"){
    S_values <- rbind(S_values, c(NA,NA))
  }
  else {
    #print(zeile)
    S <- steigung(zeile)
    # Muss in Variable gespeichert werden
    #print(c(as.character(zeile), mean_R))
    S_values <- rbind(S_values, c(as.character(zeile), S)) 
  }
}
colnames(S_values) <- c("zeilen", "S")
data_S_R <- full_join (R_values, S_values, by = "zeilen")
data_S_R %>%
  subset(S == TRUE) -> data_sub
data_final <- arrange(data_sub, desc(R))
data_final$name <- data_unique[match(data_final$zeilen, data_unique$Majority.protein.IDs),]$name
d_F <- subset(data_final, R >0.75)
spikeins <- c(d_F [,4])


LFQ_columns <- grep("LFQ.", colnames(data_unique))# get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
#fehlerbehebung
LFQ_columns
dim(data_unique)
data_unique[102:117]

plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_filt)
plot_coverage(data_filt)


#data_norm <- normalize_vsn_anna(se=data_filt, spikeins=c("CG30022", "Jafrac1", "RpS19a", "Hsp60", "alpha-Spec" )) #ohne wt
data_norm <- normalize_vsn_anna(data_filt, spikeins) #mit wt
#data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

plot_missval(data_filt)
plot_detect(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

data_diff <- test_diff(data_imp, type = "control", control = "wt")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

plot_pca(dep, x = 1, y = 2, n = 374, point_size = 4)
#plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
            k = 6, col_limit = 4, show_row_names = FALSE,
            indicate = c("condition", "replicate"))
plot_volcano(dep, contrast = "F_vs_wt", label_size = 2, add_names = TRUE)


#eigener volcanoplot
Interactome <- read.csv("Fly_Interactome_curated.csv", header=TRUE, sep=",")
data_int <- Interactome [,1] 

data_results <- get_results(dep)
data_results$Gene.ID <- mapIds(org.Dm.eg.db, keys=str_replace_all(data_results$ID, ";.*", ""),
                              column="ENSEMBL", keytype="UNIPROT",  multiVals="first")
data_results %>%
  dplyr::select(matches(".ratio")) -> data_ratio
#data_r <- log2(data_ratio)
data_rf <- tidyr::gather(data_ratio, "sample1", "ratio")
data_rf <- cbind(data_rf, data_results$Gene.ID)
colnames(data_rf) <- c("sampler", "ratio", "protein")
data_results %>%
  dplyr::select(matches(".p.val")) -> data_p_value
data_p <- -log10(data_p_value)
data_pf <- tidyr::gather(data_p, "sample", "pval")
data_pf <- cbind(data_pf, data_results$Gene.ID)
colnames(data_pf) <- c("samplep", "pval", "protein")
data_results %>%
  dplyr::select(matches(".significant")) -> data_significant
data_s <- tidyr::gather(data_significant, "sample", "pval")
data_sf <- cbind(data_s, data_results$Gene.ID)
colnames(data_sf) <- c("samples", "significant", "protein")
data_final <- full_join (data_pf, data_rf, by = "protein")
data_final <- full_join (data_final, data_sf, by = "protein")
data_final %>%
  mutate(samplep = str_replace(samplep, "._pval", " ")) -> data_plot
mit <- data.frame()
#colnames(mit) <- c( "Mitochondrial protein")
for(i in c(1:dim(data_final)[1])){
  if(is.element(data_final[i,3], data_int)){
    mit <- rbind(mit, TRUE)
  }
  else {
    mit <- rbind(mit, FALSE) 
  }
}
colnames(mit) <- c( "Mitochondrial_protein")
data_plot <- cbind(data_final, mit)
#cairo_pdf("volcano_plot.pdf", 15,5)
ggplot(data_plot, aes(x = ratio, y = pval, color=significant,  shape=Mitochondrial_protein ))+# shape= mitochondrial protein)) +
  geom_point()+#color = ifelse(data_plot$significant == TRUE, "orange", "grey50")), shape= ifelse(data_plot$Mitochondrial_protein == TRUE, 15, 16))+
  facet_wrap(~ samplep, ncol=3)+
  xlab("log2 Fold change")+
  ylab("-log10 p-value")+
  geom_line(aes(x=0), colour="#990000")+
  theme_bw()+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  geom_label_repel(data = subset(data_plot, significant== TRUE), aes(label=protein),   segment.size  = 0.2,
                   segment.color = "grey50", size=1)#data=significant, aes(x = ratio, y =pval, label=protein)) 
#dev.off()
  
  







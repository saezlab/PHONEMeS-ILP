library(readr)
library(dplyr)
library(viper)

Phospho_STY_Sites <- as.data.frame(read_delim("~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/Phospho (STY)Sites.txt", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE))

Phospho_STY_Sites$`Gene names` <- gsub("[;].*","",Phospho_STY_Sites$`Gene names`)
Phospho_STY_Sites$`Positions within proteins` <- gsub("[;].*","",Phospho_STY_Sites$`Positions within proteins`)

phospho <- as.data.frame(Phospho_STY_Sites[,c(8,150,161,167,168,169,170,171,172,173,151,152,153,154,155,156,157,158,159,160,162,163,164,165,166)])
plot(hist(phospho$`Localization prob`, breaks = 100))

phospho <- as.data.frame(Phospho_STY_Sites[,c(150,161,167,168,169,170,171,172,173,151,152,153,154,155,156,157,158,159,160,162,163,164,165,166)])

Phospho_STY_Sites$ID <- apply(Phospho_STY_Sites[ , c(6,112) ] , 1 , paste , collapse = "_" )
Phospho_STY_Sites$ID <- apply(Phospho_STY_Sites[ , c(364,2) ] , 1 , paste , collapse = "" )
phospho$ID <- Phospho_STY_Sites$ID

length(unique(phospho$ID)) - length(phospho$ID)


batches <- phospho

batches <- batches %>% group_by(ID) %>% summarise_each(funs(sum(., na.rm = TRUE)))
batches <- as.data.frame(batches)

row.names(batches) <- batches$ID
batches <- batches[,-1]

library(omicToolsTest)

batches[batches == 0] <- NA
batches <- batches[rowSums(is.na(batches)) < 22,]
names(batches) <- c("HK2_1","HK2_2","HK2_3","HK2_4","HK2_5",
                    "786-O+EV_1","786-O+EV_2","786-O+EV_3","786-O+EV_4","786-O+EV_5",
                    "786-O+VHL_1","786-O+VHL_2","786-O+VHL_3","786-O+VHL_4","786-O+VHL_5",
                    "786-M1A+EV_1","786-M1A+EV_2","786-M1A+EV_3","786-M1A+EV_4","786-M1A+EV_5",
                    "786-M1A+VHL_2","786-M1A+VHL_3","786-M1A+VHL_4","786-M1A+VHL_5")

targets <- as.data.frame(matrix(NA,ncol = 2, nrow = 24))
names(targets) <- c("sample","condition")
targets$sample <- names(batches)
targets$condition<- c(rep("HK2", 5), rep("786-O+EV",5), rep("786-O+VHL",5), rep("786-M1A+EV", 5), rep("786-M1A+VHL",4))

library(vsn)
fit <- vsnMatrix(as.matrix(batches))
meanSdPlot(fit)
batches <- as.data.frame(vsn::predict(fit,as.matrix(batches)))

library(limma)
unique(targets$condition)
limma_res <- runLimma(batches, targets, comparisons = list(c(2,-1),c(4,-1),c(4,-2),c(5,-4,-4,2)))

ttop_list <- list()
for(i in 1:4)
{
  ttop_list[[i]] <- ttopFormatter(topTable(limma_res[[1]],coef = i, number = 13670-5732, adjust.method = "fdr"))
}
View(ttop_list[[1]])
View(ttop_list[[2]])
View(ttop_list[[3]])
View(ttop_list[[4]])


write_csv(ttop_list[[1]],"~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/ttop_cancer_vs_healthy.csv")
write_csv(ttop_list[[2]],"~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/ttop_metastasis_vs_healthy.csv")
write_csv(ttop_list[[3]],"~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/ttop_metastasis_vs_cancer.csv")
write_csv(ttop_list[[4]],"~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/ttop_metastasis_resolution.csv")

url <- paste0(
  'http://omnipathdb.org/ptms?',
  'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omnipath_ptm <- download_omnipath()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

KSN_viper <- df_to_viper_regulon(KSN)

eset <- merge(ttop_list[[1]][,c(1,4)], ttop_list[[2]][,c(1,4)], by = "ID")
eset <- merge(eset,ttop_list[[3]][,c(1,4)], by = "ID")
eset <- merge(eset, ttop_list[[4]][,c(1,4)], by = "ID")

row.names(eset) <- eset$ID
eset <- eset[,-1]
names(eset) <- c("cancer_vs_healthy","metastasis_vs_healthy","metastasis_vs_cancer","metastasis_resolution")

viperRes <- as.data.frame(viper(eset = eset, regulon = KSN_viper, minsize = 5, adaptive.size = F, eset.filter = F, cores = 3))

kinase_for_phonemes_list <- list()

for(i in 1:4)
{
  viperRes <- viperRes[order(abs(viperRes[,i]), decreasing = T),]
  # kinase_for_phonemes_list[[i]] <- as.data.frame(t(viperRes[1:(length(viperRes[,1])/5),i]))
  # names(kinase_for_phonemes_list[[i]]) <- row.names(viperRes[1:(length(viperRes[,1])/5),])
  kinase_for_phonemes_list[[i]] <- row.names(viperRes[1:(length(viperRes[,1])/5),])
}

for_enio <- list()

for(i in 1:4)
{
  for_enio[[i]] <- list()
  for_enio[[i]][[1]] <- KSN
  for_enio[[i]][[2]] <- ttop_list[[i]]
  for_enio[[i]][[3]] <- kinase_for_phonemes_list[[i]]
  names(for_enio[[i]]) <- c("KSN","ttop","kinases")
}
names(for_enio) <- c("cancer_vs_healthy","metastasis_vs_healthy","metastasis_vs_cancer","metastasis_resolution")

save(file = "~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/input_PHONEMES_UD.Rdata", for_enio)

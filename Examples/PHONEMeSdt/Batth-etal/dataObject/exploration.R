## Limma analysis of the data to estimate their regulation level and significance

nicettop <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1:6)]
  return(ttop)
}

library(readr)

mmc2 <- as.data.frame(mmc2 <- read_delim("mmc2.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

mmc2$`Gene names` <- gsub("[;].*","",mmc2$`Gene names`)
mmc2$ID <- paste(mmc2$`Gene names`, mmc2$`Amino acid`, sep = "_")
mmc2$ID <- paste(mmc2$ID, mmc2$Position, sep = "")
mmc2$ID <- toupper(mmc2$ID)

# Taking only the measurements with their intensities
batches <- mmc2[,c(63, 1:44)]

library(dplyr)

batches <- batches %>% group_by(ID) %>% summarise_each(funs(sum(., na.rm = TRUE)))
batches <- as.data.frame(batches)
row.names(batches) <- batches$ID
batches <- batches[,-1]
batches[batches == 0] <- NA
batches <- batches[rowSums(is.na(batches)) != 44,]

targets <- as.data.frame(matrix(NA,nrow = 44, ncol = 2))
names(targets) <- c("sample","condition")
targets[c(1:7),2] <- "control"
targets[c(8:14),2] <- "3min_PDGF"
targets[c(15:20),2] <- "3min_IGF1"
targets[c(21:26),2] <- "3min_FGF2"
targets[c(27:31),2] <- "15min_PDGF"
targets[c(32:37),2] <- "15min_IGF1"
targets[c(38:44),2] <- "15min_FGF2"

j <- 0
for (i in 1:7)
{
  targets[i,1] <- paste(targets[i,2],i, sep = "_")
  j <- j+1
}
for (i in 1:7)
{
  targets[i+7,1] <- paste(targets[i+7,2],i, sep = "_")
  j <- j+1
}
for (i in 1:6)
{
  targets[i+14,1] <- paste(targets[i+14,2],i, sep = "_")
  j <- j+1
}
for (i in 1:6)
{
  targets[i+20,1] <- paste(targets[i+20,2],i, sep = "_")
  j <- j+1
}
for (i in 1:5)
{
  targets[i+26,1] <- paste(targets[i+26,2],i, sep = "_")
  j <- j+1
}
for (i in 1:6)
{
  targets[i+31,1] <- paste(targets[i+31,2],i, sep = "_")
  j <- j+1
}
for (i in 1:7)
{
  targets[i+37,1] <- paste(targets[i+37,2],i, sep = "_")
  j <- j+1
}

names(batches) <- targets$sample

library(vsn)

fit <- vsnMatrix(as.matrix(batches))
meanSdPlot(fit)
batches <- as.data.frame(predict(fit,as.matrix(batches)))

print(unique(targets$condition))

source("../../../../Public/limmaWrapper.R")
limmaRes <- runLimma(measurements = batches, targets = targets, comparisons = list(c(2,-1),c(3,-1),c(4,-1),c(5,-1),c(6,-1),c(7,-1)))

ttop_list <- list()

ttop_PDGF_3min <- topTable(limmaRes[[1]], coef = 1, number = length(which(complete.cases(limmaRes[[1]]$coefficients))), adjust.method = "fdr") # which is all NA's in the coefficients matrix
ttop_list[[1]] <- nicettop(ttop_PDGF_3min)
ttop_IGF1_3min <- topTable(limmaRes[[1]], coef = 2, number = length(which(complete.cases(limmaRes[[1]]$coefficients))), adjust.method = "fdr")
ttop_list[[2]] <- nicettop(ttop_IGF1_3min)
ttop_FGF2_3min <- topTable(limmaRes[[1]], coef = 3, number = length(which(complete.cases(limmaRes[[1]]$coefficients))), adjust.method = "fdr")
ttop_list[[3]] <- nicettop(ttop_FGF2_3min)
ttop_PDGF_15min <- topTable(limmaRes[[1]], coef = 4, number = length(which(complete.cases(limmaRes[[1]]$coefficients))), adjust.method = "fdr")
ttop_list[[4]] <- nicettop(ttop_PDGF_15min)
ttop_IGF1_15min <- topTable(limmaRes[[1]], coef = 5, number = length(which(complete.cases(limmaRes[[1]]$coefficients))), adjust.method = "fdr")
ttop_list[[5]] <- nicettop(ttop_IGF1_15min)
ttop_FGF2_15min <- topTable(limmaRes[[1]], coef = 6, number = length(which(complete.cases(limmaRes[[1]]$coefficients))), adjust.method = "fdr")
ttop_list[[6]] <- nicettop(ttop_FGF2_15min)

names(ttop_list) <- c("PDGF_3min", "IGF1_3min", "FGF2_3min", "PDGF_15min", "IGF1_15min", "FGF2_15min")

save(ttop_list, file = "ttop_list.RData")

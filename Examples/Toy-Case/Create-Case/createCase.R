load(file = "allD.RData")
load(file = "dataGMM.RData")

allD_Toy <- allD[1:7, ]

allD_Toy$S.AC <- c("M1", "M2", "M3", "M4", "I", "M4", "M5")
allD_Toy$S.ID <- c("M1", "M2", "M3", "M4", "I", "M4", "M5")
allD_Toy$K.AC <- c("S", "M1", "M1", "M3", "M1", "I", "I")
allD_Toy$K.ID <- c("S", "M1", "M1", "M3", "M1", "I", "I")
allD_Toy$res <- "S"
allD_Toy$pos <- "123"
allD_Toy$S.cc <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "I_S123", "M4_S123", "M5_S123")

save(allD_Toy, file = "allD_Toy.RData")

GMM.Toy <- list()
GMM.Toy.wFC <- list()

# Measurement - 1
GMM.Toy[[1]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[1]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[1]]) <- colnames(GMM[[1]])

GMM.Toy[[1]][, 1] <- c("-10", "10", "-10")
GMM.Toy[[1]][, 2] <- c("P", "C", "P")
GMM.Toy[[1]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy[[1]][, 4] <- c("OK", "FP", "OK")

GMM.Toy.wFC[[1]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[1]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[1]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[1]][, 1] <- c("-10", "10", "-10")
GMM.Toy.wFC[[1]][, 2] <- c("P", "C", "P")
GMM.Toy.wFC[[1]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy.wFC[[1]][, 4] <- c("OK", "FP", "OK")
GMM.Toy.wFC[[1]][, 5] <- c("2", "0", "2")

# Measurement - 2
GMM.Toy[[2]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[2]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[2]]) <- colnames(GMM[[1]])

GMM.Toy[[2]][, 1] <- c("-10", "10", "-10")
GMM.Toy[[2]][, 2] <- c("P", "C", "P")
GMM.Toy[[2]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy[[2]][, 4] <- c("OK", "FP", "OK")

GMM.Toy.wFC[[2]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[2]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[2]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[2]][, 1] <- c("-10", "10", "-10")
GMM.Toy.wFC[[2]][, 2] <- c("P", "C", "P")
GMM.Toy.wFC[[2]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy.wFC[[2]][, 4] <- c("OK", "FP", "OK")
GMM.Toy.wFC[[2]][, 5] <- c("2", "0", "2")

# Measurement - 3
GMM.Toy[[3]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[3]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[3]]) <- colnames(GMM[[1]])

GMM.Toy[[3]][, 1] <- c("-10", "10", "-10")
GMM.Toy[[3]][, 2] <- c("P", "C", "P")
GMM.Toy[[3]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy[[3]][, 4] <- c("OK", "FP", "OK")

GMM.Toy.wFC[[3]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[3]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[3]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[3]][, 1] <- c("-10", "10", "-10")
GMM.Toy.wFC[[3]][, 2] <- c("P", "C", "P")
GMM.Toy.wFC[[3]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy.wFC[[3]][, 4] <- c("OK", "FP", "OK")
GMM.Toy.wFC[[3]][, 5] <- c("2", "0", "2")

# Measurement - 4
GMM.Toy[[4]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[4]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[4]]) <- colnames(GMM[[1]])

GMM.Toy[[4]][, 1] <- c("10", "-10", "-10")
GMM.Toy[[4]][, 2] <- c("C", "P", "P")
GMM.Toy[[4]][, 3] <- c("1", "0.01", "0.01")
GMM.Toy[[4]][, 4] <- c("FP", "OK", "OK")

GMM.Toy.wFC[[4]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[4]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[4]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[4]][, 1] <- c("10", "-10", "-10")
GMM.Toy.wFC[[4]][, 2] <- c("C", "P", "P")
GMM.Toy.wFC[[4]][, 3] <- c("1", "0.01", "0.01")
GMM.Toy.wFC[[4]][, 4] <- c("FP", "OK", "OK")
GMM.Toy.wFC[[4]][, 5] <- c("0", "2", "2")

# Measurement - 5
GMM.Toy[[5]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[5]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[5]]) <- colnames(GMM[[1]])

GMM.Toy[[5]][, 1] <- c("10", "-10", "-10")
GMM.Toy[[5]][, 2] <- c("C", "P", "P")
GMM.Toy[[5]][, 3] <- c("1", "0.01", "0.01")
GMM.Toy[[5]][, 4] <- c("FP", "OK", "OK")

GMM.Toy.wFC[[5]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[5]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[5]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[5]][, 1] <- c("10", "-10", "-10")
GMM.Toy.wFC[[5]][, 2] <- c("C", "P", "P")
GMM.Toy.wFC[[5]][, 3] <- c("1", "0.01", "0.01")
GMM.Toy.wFC[[5]][, 4] <- c("FP", "OK", "OK")
GMM.Toy.wFC[[5]][, 5] <- c("0", "2", "2")

names(GMM.Toy) <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")
names(GMM.Toy.wFC) <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")

GMM.Toy.ID <- GMM.ID[1:5, ]
GMM.Toy.ID$dataID <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")
GMM.Toy.ID$UPID <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")
GMM.Toy.ID$S.cc <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")

save(file="dataObjects_Toy.RData", 
     list=c("GMM.Toy.ID","GMM.Toy", "GMM.Toy.wFC"))


#######################################################################################
load(file = "allD.RData")
load(file = "dataGMM.RData")

allD_Toy <- allD[1:7, ]

allD_Toy$S.AC <- c("M1", "M2", "M3", "M4", "I", "M4", "M5")
allD_Toy$S.ID <- c("M1", "M2", "M3", "M4", "I", "M4", "M5")
allD_Toy$K.AC <- c("S", "M1", "M1", "M3", "M1", "I", "I")
allD_Toy$K.ID <- c("S", "M1", "M1", "M3", "M1", "I", "I")
allD_Toy$res <- "S"
allD_Toy$pos <- "123"
allD_Toy$S.cc <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "I_S123", "M4_S123", "M5_S123")

save(allD_Toy, file = "allD_Toy.RData")

GMM.Toy <- list()
GMM.Toy.wFC <- list()

# Measurement - 1
GMM.Toy[[1]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[1]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[1]]) <- colnames(GMM[[1]])

GMM.Toy[[1]][, 1] <- c("-10", "10", "-10")
GMM.Toy[[1]][, 2] <- c("P", "C", "P")
GMM.Toy[[1]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy[[1]][, 4] <- c("OK", "FP", "OK")

GMM.Toy.wFC[[1]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[1]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[1]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[1]][, 1] <- c("-10", "10", "-10")
GMM.Toy.wFC[[1]][, 2] <- c("P", "C", "P")
GMM.Toy.wFC[[1]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy.wFC[[1]][, 4] <- c("OK", "FP", "OK")
GMM.Toy.wFC[[1]][, 5] <- c("2", "0", "2")

# Measurement - 2
GMM.Toy[[2]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[2]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[2]]) <- colnames(GMM[[1]])

GMM.Toy[[2]][, 1] <- c("-10", "10", "-10")
GMM.Toy[[2]][, 2] <- c("P", "C", "P")
GMM.Toy[[2]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy[[2]][, 4] <- c("OK", "FP", "OK")

GMM.Toy.wFC[[2]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[2]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[2]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[2]][, 1] <- c("-10", "10", "-10")
GMM.Toy.wFC[[2]][, 2] <- c("P", "C", "P")
GMM.Toy.wFC[[2]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy.wFC[[2]][, 4] <- c("OK", "FP", "OK")
GMM.Toy.wFC[[2]][, 5] <- c("2", "0", "2")

# Measurement - 3
GMM.Toy[[3]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[3]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[3]]) <- colnames(GMM[[1]])

GMM.Toy[[3]][, 1] <- c("-10", "10", "-10")
GMM.Toy[[3]][, 2] <- c("P", "C", "P")
GMM.Toy[[3]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy[[3]][, 4] <- c("OK", "FP", "OK")

GMM.Toy.wFC[[3]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[3]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[3]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[3]][, 1] <- c("-10", "10", "-10")
GMM.Toy.wFC[[3]][, 2] <- c("P", "C", "P")
GMM.Toy.wFC[[3]][, 3] <- c("0.01", "1", "0.01")
GMM.Toy.wFC[[3]][, 4] <- c("OK", "FP", "OK")
GMM.Toy.wFC[[3]][, 5] <- c("2", "0", "2")

# Measurement - 4
GMM.Toy[[4]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[4]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[4]]) <- colnames(GMM[[1]])

GMM.Toy[[4]][, 1] <- c("10", "-10", "10")
GMM.Toy[[4]][, 2] <- c("C", "P", "C")
GMM.Toy[[4]][, 3] <- c("1", "0.01", "1")
GMM.Toy[[4]][, 4] <- c("FP", "OK", "FP")

GMM.Toy.wFC[[4]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[4]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[4]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[4]][, 1] <- c("10", "-10", "10")
GMM.Toy.wFC[[4]][, 2] <- c("C", "P", "C")
GMM.Toy.wFC[[4]][, 3] <- c("1", "0.01", "1")
GMM.Toy.wFC[[4]][, 4] <- c("FP", "OK", "FP")
GMM.Toy.wFC[[4]][, 5] <- c("0", "2", "0")

# Measurement - 5
GMM.Toy[[5]] <- matrix(data = , nrow = 3, ncol = 4)
rownames(GMM.Toy[[5]]) <- paste0("cond", 1:3)
colnames(GMM.Toy[[5]]) <- colnames(GMM[[1]])

GMM.Toy[[5]][, 1] <- c("10", "-10", "10")
GMM.Toy[[5]][, 2] <- c("C", "P", "C")
GMM.Toy[[5]][, 3] <- c("1", "0.01", "1")
GMM.Toy[[5]][, 4] <- c("FP", "OK", "FP")

GMM.Toy.wFC[[5]] <- matrix(data = , nrow = 3, ncol = 5)
rownames(GMM.Toy.wFC[[5]]) <- paste0("cond", 1:3)
colnames(GMM.Toy.wFC[[5]]) <- colnames(GMM.wFC[[1]])

GMM.Toy.wFC[[5]][, 1] <- c("10", "-10", "10")
GMM.Toy.wFC[[5]][, 2] <- c("C", "P", "C")
GMM.Toy.wFC[[5]][, 3] <- c("1", "0.01", "1")
GMM.Toy.wFC[[5]][, 4] <- c("FP", "OK", "FP")
GMM.Toy.wFC[[5]][, 5] <- c("0", "2", "0")

names(GMM.Toy) <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")
names(GMM.Toy.wFC) <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")

GMM.Toy.ID <- GMM.ID[1:5, ]
GMM.Toy.ID$dataID <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")
GMM.Toy.ID$UPID <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")
GMM.Toy.ID$S.cc <- c("M1_S123", "M2_S123", "M3_S123", "M4_S123", "M5_S123")

save(file="dataObjects_Toy_2.RData", 
     list=c("GMM.Toy.ID","GMM.Toy", "GMM.Toy.wFC"))
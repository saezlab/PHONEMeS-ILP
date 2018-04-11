GMM.ID <- dataGMM@IDmap
sites <- intersect(GMM.ID$S.cc, unique(c(resultsSIF[, 1], resultsSIF[, 3])))
nodesAttributes <- matrix(data = , nrow = length(as.character(unlist(targets.P)))+length(sites), ncol = 2)
colnames(nodesAttributes) <- c("Species", "nodesP")
nodesAttributes[1:length(as.character(unlist(targets.P))), 1] <- as.character(unlist(targets.P))
nodesAttributes[1:length(as.character(unlist(targets.P))), 2] <- "D"
nodesAttributes[(length(as.character(unlist(targets.P)))+1):nrow(nodesAttributes), 1] <- sites
nodesAttributes[(length(as.character(unlist(targets.P)))+1):nrow(nodesAttributes), 2] <- "P"
write.table(x = nodesAttributes, file = "nodesAttributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
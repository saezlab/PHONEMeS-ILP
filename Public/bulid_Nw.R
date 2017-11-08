build_Nw<-function(data.On,targets.On, bg,
                  nK=c("all","no", "drugs2data", "data")){
  pkn<-build_PKN(data.On, targets.On, bg,nK=nK) 
  pknList<-PKNlist(pkn, targets.On, data.On)
  idx <- which(pknList@interactions$S.ID==pknList@interactions$K.ID)
  rem <- c()
  if(length(idx) > 0){
    for(i in 1:length(idx)){
      pknList@interactions <- pknList@interactions[-intersect(which(pknList@interactions$K.ID==pknList@interactions[idx[i], ]$S.cc), which(pknList@interactions$S.cc==pknList@interactions[idx[i], ]$K.ID)), ]
    }
  }
  return(pknList)   
}

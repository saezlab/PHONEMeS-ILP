build_Nw_Inv<-function(data.On,targets.On, bg,
                  nK=c("all","no", "drugs2data", "data")){
  pkn<-build_PKN_Inv(data.On, targets.On, bg,nK=nK)
  if(is.null(pkn)){
    
    return(NULL)
    
  }
  else{
    
    View(pkn@interactions)
    pknList<-PKN_list(pkn, targets.On, data.On)
    
    pkn <- pknList@interactions[complete.cases(pknList@interactions), ]
    
    sites <- unique(pkn$K.ID)
    
    for(i in 1:length(sites)){
      
      if(!(sites[i]%in%pknList@interactions$S.cc)){
        
        idx1 <- which(pknList@interactions$K.ID==sites[i])
        idx2 <- which(pknList@species==sites[i])
        
        pknList@interactions <- pknList@interactions[-idx1, ]
        pknList@species <- pknList@species[-idx2]
        
      }
      
    }
    
    
    return(pknList)
    
  }
  
}

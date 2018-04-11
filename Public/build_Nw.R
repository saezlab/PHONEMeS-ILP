build_Nw<-function(data.On,targets.On, bg,
                  nK=c("all","no", "drugs2data", "data")){
  pkn<-build_PKN(data.On, targets.On, bg,nK=nK) 
  pknList<-PKNlist(pkn, targets.On, data.On)           
  return(pknList)   
}
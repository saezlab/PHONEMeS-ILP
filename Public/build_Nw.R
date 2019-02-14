#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E. Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: https://saezlab.github.io/PHONEMeS/
#
##############################################################################
# $Id$

# This function makes a starting network object, it just calls the build_PKN and PKNlist functions

build_Nw<-function(data.On,targets.On, bg,
                   nK=c("all","no", "drugs2data", "data")){
  pkn<-build_PKN(data.On, targets.On, bg,nK=nK) 
  if(!is.null(pkn)){
    # pknList<-PKNlist(pkn, targets.On, data.On)
    
    interactions <- pkn@interactions
    
    idxNA <- which(is.na(interactions$S.AC))
    
    kinases <- unique(interactions$K.ID[-idxNA])
    
    kinase2remove <- c()
    for(ii in 1:length(kinases)){
      
      if((kinases[ii]%in%interactions$S.cc)==FALSE){
        
        kinase2remove <- c(kinase2remove, kinases[ii])
        
      }
      
    }
    
    targets <- unique(unlist(targets.On))
    
    kinase2remove <- setdiff(kinase2remove, targets)
    
    if(length(kinase2remove)>0){
      
      idx2remove <- which(pkn@interactions$K.ID%in%kinase2remove)
      pkn@interactions <- pkn@interactions[-idx2remove, ]
      
    }
    
    return(pkn)
  } else {
    
    return(NULL)
  }
}
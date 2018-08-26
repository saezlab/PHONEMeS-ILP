#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmai.com
#
##############################################################################
# 12:56 23/03/2018
# This function makes a starting network object, it just calls the build_PKN and PKNlist functions
build_Nw<-function(data.On,targets.On, bg,
                  nK=c("all","no", "drugs2data", "data")){
  pkn<-build_PKN(data.On, targets.On, bg,nK=nK) 
  # pknList<-PKNlist(pkn, targets.On, data.On)
  
  return(pkn)   
}
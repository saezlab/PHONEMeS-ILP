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
  # pknList<-PKNlist(pkn, targets.On, data.On)
  
  return(pkn)   
}
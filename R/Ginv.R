# 
# Copyright 2001 Mayo Foundation for Medical Education and Research. 
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
# 02111-1307, USA.
# 
# 
Ginv<-function(x){

savesvd<-svd(x)                                  
U.svd<-savesvd$u                                   
V.svd<-savesvd$v
d.svd<-savesvd$d

eps<-.000001                                       
maxd<-max(d.svd)                                  
w<-ifelse((d.svd/maxd) < eps, rep(0,length(d.svd)), 1/d.svd)
df<-sum(d.svd/maxd >= eps)

Ginv<-V.svd %*% diag(w) %*% t(U.svd)

list(Ginv=Ginv,rank=df)
}

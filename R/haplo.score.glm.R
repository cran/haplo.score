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
haplo.score.glm <- function(y,mu,a,v,x.adj,nreps,x.post, post, x){

   u.mtx  <- (y-mu)*x.post / a
   u.score <- apply(u.mtx,2,sum)

   # Var matrix for x.adj covariates

   v.11 <- t(x.adj * v) %*% x.adj

   # Var matrix for covar(x.adj, x.post)
   v.21 <- t(x.post) %*% (x.adj * v)

   # Var matrix for haplo scores
   res <- ( (y - mu)/a ) ^2
   t1 <- rep( (v-res) ,nreps) * post
   v.22 <- t(x*t1) %*% x + t(u.mtx) %*% u.mtx
           
   # Var matrix for haplo scores, adjusted for x.adj

   v.score <- v.22 - v.21 %*% solve(v.11) %*% t(v.21) 
    
   return(list(u.score=u.score, v.score=v.score))

}


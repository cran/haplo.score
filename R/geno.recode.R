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
geno.recode <- function(geno, miss.val=0){

# loop over all loci to recode pairs of allele vectors

n.loci <- ncol(geno)/2
alist <- vector("list",n.loci)
grec <- NULL

for(i in 1:n.loci){
   t <- (i-1)*2 + 1
   tmp <- allele.recode(geno[,t],geno[,t+1], miss.val=miss.val)
   grec <- cbind(grec,tmp$a1,tmp$a2)
   alist[[i]] <- list(allele=tmp$allele.label)
}

return(list(grec=grec, alist=alist))

}




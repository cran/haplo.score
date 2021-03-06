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
allele.recode <- function(a1,a2,miss.val=NA){

n <- length(a1)

# convert factor to character:
if(is.factor(a1)) a1 <- as.character(a1)
if(is.factor(a2)) a2 <- as.character(a2)

# test if either a1 or a2 are character. If either, then 
# treat allele levels as ordered charcter. If neither, then
# treat allele levels as ordered numeric.

is.ch <- is.character(a1) | is.character(a2)

if (is.ch){
   t <- factor(c(a1,a2),exclude=miss.val)
  }


if (!is.ch) {
   lev <- sort(unique(c(a1,a2)))
   t <- factor(c(a1,a2),levels=lev,exclude=miss.val)
  }

allele.label <- levels(t)
t <- as.numeric(t)
a1 <- t[1:n]
a2 <- t[(n+1):(2*n)]

return(list(a1=a1,a2=a2,allele.label=allele.label))

}


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
ungcode <- function(x) {

# uncodes unique index for a one-locus genotype 
# based on formula from:
#
#   Zykin D., Zhivotovsky L. and Weir B, "Exact tests for
#   association between alleles at arbitrary numbers of
#   loci"  Human Identification: The Use of DNA Markers,
#   169-178, 1995.

  i <- 1 + floor((sqrt(8*x+1)-1)/2)
  j <- x - i*(i-1)/2
  
  # the following 2 lines were added as a patch to make this ungcode work:
  i <- ifelse(j==0,i-1,i)
  j <- ifelse(j==0,i,j)

  return(cbind(j,i))
}

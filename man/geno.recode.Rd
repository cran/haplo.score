% 
% Copyright 2001 Mayo Foundation for Medical Education and Research. 
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
% 02111-1307, USA.
% 
% 
\name{geno.recode}
\alias{geno.recode}
\title{
  Recode Genotypes
}
\description{
For all loci as pairs of columns in a matrix, recode
alleles
}
\usage{
geno.recode(geno, miss.val=0)
}
\arguments{


\item{geno}{
Matrix of alleles, such that each locus has a  pair  of
adjacent  columns  of  alleles.  If there are K loci, then
ncol(geno) = 2*K. Rows represent alleles for each subject.
}
\item{miss.val }{
Vector of codes for missing values of alleles.
}
}
\value{
List with components:

\item{grec}{
Matrix of recoded alleles - see allele.recode
}
\item{alist}{
List of allele labels.  For K loci, there are K
components in the list, and the kth  component is a vector
of sorted unique allele labels for the kth locus.
}
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
allele.recode
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

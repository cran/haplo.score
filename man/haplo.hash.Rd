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
\name{haplo.hash}
\alias{haplo.hash}
\title{
  Integer Rank Codes for Haplotypes
}
\description{
Create a vector of integer codes for the input  matrix  of
haplotypes.   The haplotypes in the input matrix are
converted to character strings, and if there are C unique
strings, the integer codes for the haplotypes will be 1,
2, ..., C.
}
\usage{
haplo.hash(hap)
}
\arguments{



\item{hap}{
A matrix of haplotypes. If there are N haplotypes for K
loci, hap have dimensions N x K.
}
}
\value{
List with elements:

\item{hash}{
Vector of integer codes for the input data  (hap).  The
value of hash is the row number of the unique haplotypes
given in the returned matrix hap.mtx.
}
\item{hap.mtx}{
Matrix of unique haplotypes.
}
}
\section{Side Effects}{

}
\details{
The alleles that make  up  each  row  in  hap  are  pasted
together  as character strings, and the unique strings are
sorted so that the rank order of  the  sorted  strings  is
used as the integer code for the unique haplotypes.
}
\section{References}{

}
\seealso{
haplo.em
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

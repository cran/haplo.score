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
\name{ungcode}
\alias{ungcode}
\title{
  Uncode Numeric Genotype Codes to Alleles
}
\description{
Convert a vector with numeric codes that represents the
geneotypes to two vectors of alleles for genotypes
(alleles coded as rank order).
}
\usage{
ungcode(x)
}
\arguments{


\item{x}{
Vector of integer codes for genotypes.
}
None.
}
\value{
A matrix, with the first column for the "first" allele,
and the second column for the "second" allele.
}
\section{Side Effects}{

}
\details{
The input vector of genotypes must have been coded by  the
function gcode.
}
\section{References}{

Zykin D., Zhivotovsky L. and Weir B,  "Exact tests for
association  between alleles at arbitrary numbers of loci"
Human Identification: The Use of DNA Markers, 169-178, 1995.


The algorithm by Zykin  was  not  quite  correct,  and  we
needed  to  add  a  patch  in  the  code  (see comments in
ungcode).

}
\seealso{
gcode, allele.recode, geno.recode
}
\examples{
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

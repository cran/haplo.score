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
\name{gcode}
\alias{gcode}
\title{
  Numeric Code for Genotypes
}
\description{
For a one-locus genotype, composed of two vectors of
alleles, create a single vector of numeric codes that
represent the genotypes.
}
\usage{
gcode(a1, a2)
}
\arguments{



\item{a1}{
Vector of "first" alleles.
}
\item{a2}{
Vector of "second" alleles.
}
}
\value{
Vector of integer codes for genotypes.
}
\section{Side Effects}{

}
\details{
Alleles are assumed to be coded as rank order (see
allele.recode).
}
\section{References}{
Zykin D., Zhivotovsky L. and  Weir  B,  "Exact tests for
association  between alleles at arbitrary numbers of loci"
Human Identification: The Use of DNA Markers, 169-178, 1995.
}
\seealso{
ungcode, allele.recode, geno.recode
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

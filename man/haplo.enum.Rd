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
\name{haplo.enum}
\alias{haplo.enum}
\title{
  Enumerate all possible pairs of haplotypes  that  are  consistent
with a set of un-phased multilocus markers
}
\description{
Given input vectors h1 and h2  (one  possible  pair  of
haplotypes), enumerate all possible pairs of haplotypes,
and return enumerated pairs in matrices h1 and h2.
}
\usage{
haplo.enum(h1, h2)
}
\arguments{



\item{h1}{
A haplotype vector of length K (K = number of loci), with
values of alleles for each of the loci.
}
\item{h2}{
A second haplotype vector.
}
}
\value{
List with components:

\item{h1}{
A matrix of enumerated haplotypes.  If there are N
enumerations, h1 will have dimension N x K.
}
\item{h2}{
Similar to h1, a matrix of enumerated haplotypes for the
second members of the pairs of haplotypes.
}
}
\section{Side Effects}{

}
\details{
  For a pair of haplotypes, if there are H sites that are heterozygous,
  then there are \eqn{2^{(H-1)}}{2^(H-1)} possible pairs to
  enumerate. To achieve this, the algorithm moves across the loci that
  are heterozygous (after the 1st heterozygous locus), flipping alleles
  at heterozygous locations to enumerate all possible pairs of
  haplotpes, and appending results as rows of the output matrices h1,
  and h2.
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

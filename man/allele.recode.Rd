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
\name{allele.recode}
\alias{allele.recode}
\title{
  Recode allele values to integer ranks
}
\description{
Genotypes for subjects represented by a pair  of  vectors,
with the vectors containing allele values (either numeric,
factor, or character), are recoded to the  rank  order  of
allele values.
}
\usage{
allele.recode(a1, a2, miss.val=NA)
}
\arguments{


\item{a1}{
Vector of "first" alleles.
}
\item{a2}{
Vector of "second" alleles.
}
\item{miss.val }{
Vector of missing value codes for alleles.
}
}
\value{
List with components:

\item{a1}{
Vector of recoded "first" alleles.
}
\item{a2}{
Recode of recoded "second" alleles.
}
\item{allele.label}{
Vector of labels for unique alleles.
}
}
\section{Side Effects}{

}
\details{
If alleles are numeric, they are recoded to the rank order
of the alleles.  If the alleles are factor or character,
they are recoded to interger values that correspond to the
indices of the sorted values of the unique alleles, but
sorted as character values.
}
\section{References}{

}
\seealso{
geno.recode
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

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
\name{print.haplo.score}
\alias{print.haplo.score}
\title{
  Print a haplo.score object

}
\description{
Method function to print a class of type haplo.score
}
\usage{
print.haplo.score(x, ...)
}
\arguments{
\item{x}{ The object returned from haplo.score (which has class haplo.score).}
\item{...}{Optional argunents.}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
This is a print method function used to print information from
haplo.score class, with haplotype-specific information given in a
table. Because haplo.score is a class, the generic print function 
can be used, which in turn calls this print.haplo.score function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland  GA.
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Submitted to Amer J Hum Genet.
}
\seealso{
haplo.score
}
\examples{
\dontrun{
save <- haplo.score(y, geno, trait.type = "gaussian")

# Example illustrating generic print function:
print(save)

# Example illustrating specific method print function:
print.haplo.score(save)
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

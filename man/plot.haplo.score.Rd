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
\name{plot.haplo.score}
\alias{plot.haplo.score}
\title{
  Plot Haplotype Frequencies versus Haplotype Score Statistics
}
\description{
Method function to plot a class of type haplo.score
}
\usage{
plot.haplo.score(x, ...)
}
\arguments{
\item{x}{The object returned from haplo.score (which has class haplo.score).}
\item{...}{Optional arguments}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
This is a plot method function used to plot haplotype frequencies on
the x-axis and haplotype-specific scores on the y-axis. Because
haplo.score is a class, the generic plot function 
can be used, which in turn calls this plot.haplo.score function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Submitted to Amer J Hum Genet.
}
\seealso{
haplo.score
}
\examples{
\dontrun{
save <- haplo.score(y, geno, trait.type = "gaussian")

# Example illustrating generic plot function:
plot(save)

# Example illustrating specific method plot function:
plot.haplo.score(save)
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

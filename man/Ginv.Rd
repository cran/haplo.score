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
\name{Ginv}
\alias{Ginv}
\title{
  Compute Generalized Inverse of Input Matrix
}
\description{
Singular value decomposition (svd) is used to compute a
generalized inverse of input matrix.
}
\usage{
Ginv(x)
}
\arguments{



\item{x}{
A matrix.
}
}
\value{
List with components:

\item{Ginv}{
Generalized inverse of x.
}
\item{rank}{
Rank of matrix x.
}
}
\section{Side Effects}{

}
\details{
The function svd is used to compute the singular values of
the input matrix, and the rank of the matrix is determined
by the number of singular values  that  are  at  least  as
large   as  max(svd)*eps,  where  eps  is  a  small  value
(currently eps = .000001).
}
\section{References}{
Press  WH,  Teukolsky  SA,  Vetterling  WT,  Flannery  BP.
Numerical  recipes  in C. The art of scientific computing.
2nd ed. Cambridge University Press,  Cambridge.1992.  page
61.
}
\seealso{
svd
}
\examples{
\dontrun{
save <- Ginv(x)
ginv.x <- save$Ginv
rank.x <- save$rank
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

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
\name{haplo.score.glm}
\alias{haplo.score.glm}
\title{
  Compute haplotype score statistics for GLM
}
\description{
This function is used by haplo.score when analyzing traits
by a GLM score.
}
\usage{
haplo.score.glm(y, mu, a, v, x.adj, nreps, x.post, post, x)
}
\arguments{


\item{y}{
Vector of trait values.
}
\item{mu}{
Expected value of y.
}
\item{a}{
scale parameter
}
\item{v}{
v= b''/a for a GLM.
}
\item{x.adj }{
Matrix of non-genetic covariates used  to  adjust  the
score  statistics. Note that intercept should  be included
in this matrix.
}
\item{nreps}{
Vector for the count of haplotype  pairs  that  map  to
each subject's marker genotypes (see haplo.em).
}
\item{x.post }{
Matrix for posterior mean of x per subject.
}
\item{post}{
Vector of posterior probabilities of pairs of
haplotypes for a person, given thier marker phenotypes
(see haplo.em).
}
\item{x}{
Matrix  of  scores  for  enumerated  haplotypes  for  each
subject,  with  elements  0,  1,  2  (counts  of  specific
haplotypes).
}
None.
}
\value{
List with components:

\item{u.score}{
Vector of scores for the chosen haplotypes
}
\item{v.score}{
Covariance matrix for u.score
}
}
\section{Side Effects}{

}
\details{
Using posterior probabilities of pairs of haplotypes,  the
"design" matrix for the haplotype effects, and the GLM
residuals, compute the score vector and its variance
matrix, adjusted for the non-genetic covariates.
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

}
\keyword{}
% docclass is functionhaplo
% Converted by Sd2Rd version 1.21.

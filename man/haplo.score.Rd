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
\name{haplo.score}
\alias{haplo.score}
\title{
  Score Statistics for Association of Traits with Haplotypes
}
\description{
Compute score statistics to evaluate the association of  a
trait  with  haplotypes, when linkage phase is unknown and
diploid marker phenotypes  are  observed  among  unrelated
subjects. For now, only autosomal loci are considered.
}
\usage{
haplo.score(y, geno, trait.type="gaussian", offset=NA, x.adj=NA, skip.haplo=0.005, locus.label=NA, miss.val=0, n.sim=0)
}
\arguments{


\item{y}{
Vector of trait values. For  trait.type  =  "binomial",  y
must have values of 1 for event, 0 for no event.
}
\item{geno}{
Matrix of alleles, such that each locus has a  pair  of
adjacent  columns  of  alleles,  and  the order of columns
corresponds to the order of  loci  on  a  chromosome.   If
there  are  K  loci, then ncol(geno) = 2*K. Rows represent
alleles for each subject.
}
\item{trait.type }{
Character string  defining  type  of  trait,  with
values of "gaussian", "binomial", "poisson", "ordinal".
}
\item{offset}{
Vector of offset when trait.type = "poisson"
}
\item{x.adj }{
Matrix of non-genetic covariates  used  to  adjust  the
score  statistics.  Note  that  intercept  should  not  be
included, as it will be added in this function.
}
\item{skip.haplo }{
Skip score statistics for haplotypes with frequencies < skip.haplo
}
\item{locus.label }{
Vector of labels for loci, of length K (see definition of geno matrix).
}
\item{miss.val }{
Vector of codes for missing values of alleles.
}
\item{n.sim }{
Number of simulations for empirical p-values.  If
n.sim=0, no empirical p-values are computed.
}
None.
}
\value{
List with the following components:

\item{score.global}{
Global statistic to test association of trait
with haplotypes that have frequencies >= skip.haplo.
}
\item{df}{
Degrees of freedom for score.global.
}
\item{score.global.p}{
P-value of score.global based on chi-square
distribution, with degrees of freedom equal to df.
}
\item{score.global.p.sim}{
P-value of score.global based on simulations
(set equal to NA when n.sim=0).
}
\item{score.haplo}{
Vector of score statistics for individual
haplotypes that have frequencies >= skip.haplo.
}
\item{score.haplo.p}{
Vector of p-values for score.haplo, based on a
chi-square distribution with 1 df.
}
\item{score.haplo.p.sim}{
Vector of p-values for score.haplo, based on
simulations (set equal to NA when n.sim=0).
}
\item{score.max.p.sim}{
P-value  of  maximum  score.haplo, based on
simulations  (set equal to NA when n.sim=0).
}
\item{haplotype}{
Matrix of hapoltypes  analyzed.  The ith row of
haplotype corresponds to the ith item of score.haplo,
score.haplo.p, and score.haplo.p.sim.
}
\item{hap.prob}{
Vector of haplotype probabilies, corresponding to the
haplotypes in the matrix haplotype.
}
\item{locus.label}{
Vector of labels for loci, of length  K  (same  as
input argument).
}
\item{n.sim}{
Number of simulations.
}
\item{n.val.global}{
Number of valid simulated global statistics.
}
\item{n.val.haplo}{
Number of valid simulated score statistics
(score.haplo) for individual haplotypes.
}
}
\section{Side Effects}{

}
\details{
Using haplo.enum, all possible haplotypes consistent with
the observed marker phenotypes are enumerated, and
haplo.em is used to compute the maximum likelihood
estimates  of  the  haplotype frequencies, as well as the
posterior probabilities of the pairs of haplotypes for
each subject. These posterior probabilities are used to
compute the score statistics for the association of
(ambiguous) haplotypes  with  traits. The glm function is
used to compute residuals of the regression of the trait
on the non-genetic covariates.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Submitted to Amer J Hum Genet.
}
\seealso{
plot.haplo.score, print.haplo.score
}
\examples{
\dontrun{
# For quantitative, normally distributed trait:

  haplo.score(y, geno, trait.type = "gaussian")

# For a  binary trait

  haplo.score(y, geno, trait.type = "binomial")

# For a binary trait, adjusted for sex and age:

  x <- cbind(sex, age)  
  haplo.score(y, geno, trait.type = "binomial", x.adj=x)

# For poisson trait:

  haplo.score(y, geno, trait.type="poisson", offset = personYears)

# For ordinal trait:

  haplo.score(y, geno, trait.type="ordinal")
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

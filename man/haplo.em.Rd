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
\name{haplo.em}
\alias{haplo.em}
\title{
  EM Computation of Haplotype Probabilities
}
\description{
For genotypes measured on unrelated subjects, with linkage
phase  unknown,  compute  maximum  likelihood estimates of
haplotype probabilities. Because linkage phase is unknown,
there  may be  more than one pair of haplotypes that are
consistent with  the oberved marker phenotypes, so
posterior probabilities of pairs of haplotypes for each
subject are also computed.
}
\usage{
haplo.em(geno, locus.label=NA, converge.eps=1e-06, maxiter=500)
}
\arguments{


\item{geno}{
Matrix of alleles, such that each locus has a  pair  of
adjacent  columns  of  alleles,  and  the order of columns
corresponds to the order of  loci  on  a  chromosome.   If
there  are  K  loci, then ncol(geno) = 2*K. Rows represent
alleles for each subject.
}
\item{locus.label }{
Vector of  labels  for  loci,  of  length  K  (see
definition of geno matrix).
}
\item{converge.eps }{
Convergence criterion, based on absolute  change
in log likelihood (lnlike).
}
\item{maxiter}{
Maximum number of iterations of EM.
}
}
\value{
List with components:

\item{converge}{
Indicator of convergence of the EM algorithm
(1=converge, 0 = failed).
}
\item{niter}{
Number of iterations completed in the EM alogrithm.
}
\item{locus.info}{
A list with  a  component for each locus.  Each
component is also a list, and  the  items of a locus-
specific list are the locus name and a vector for the
unique alleles for the locus.
}
\item{locus.label}{
Vector of  labels  for  loci,  of  length  K  (see
definition of input values).
}
\item{haplotype}{
Matrix of unique haplotypes. Each row represents a
unique  haplotype, and the number of columns is the number
of loci.
}
\item{hap.prob}{
Vector of mle's of haplotype probabilities.  The ith
element of hap.prob corresponds to the ith row of 
haplotype.
}
\item{hap.prob.noLD}{
Similar to hap.prob, but assuming no linkage
disequilibrium.
}
\item{lnlike}{
Value of lnlike at last EM iteration (maximum lnlike if converged).
}
\item{lr}{
Likelihood ratio statistic to test no linkage
disequilibrium among all loci.
}
\item{indx.subj}{
Vector for index of subjects, after  expanding  to
all possible  pairs  of  haplotypes  for  each person. If
indx=i, then i is the ith row of input matrix geno. If the
ith subject has  n possible  pairs  of haplotypes that
correspond to their marker phenotype, then i is repeated n
times.
}
\item{nreps}{
Vector for the count of haplotype pairs that map to
each subject's marker genotypes.
}
\item{hap1code}{
Vector of codes for each subject's first haplotype.
The values in hap1code are the row numbers of the unique
haplotypes in the returned matrix haplotype.
}
\item{hap2code}{
Similar to hap1code, but for  each  subject's  second
haplotype.
}
\item{post}{
Vector of posterior probabilities of pairs of
haplotypes for a person, given thier marker phenotypes.
}
}
\section{Side Effects}{

}
\details{
The input data are arranged  as  a  matrix,  with  N  rows
representing N subjects, and 2K columns representing pairs
of alleles for K loci whose phase is  unknown.  The  input
data  matrix  is  reduced to the distinguishable un-phased
multilocus marker phenotypes, along with their counts. For
each  distinguishable  phenotype,  all  possible  pairs of
haplotypes are enumerated. Maximum likelihood  estimation,
implemented    by    the   expectation-maximization   (EM)
algorithm, proceeds by assuming Hardy-Weinberg proportions
of underlying genotypes, so that the probability of a pair
of haplotypes is the product of their probabilities (times
2  if  haplotypes differ), and then relative probabilities
are assigned to the list of possible underlying  pairs  of
haplotypes for each genotype. The haplotypes are "counted"
from the enumerated list of  all  possibilities,  but  the
relative  probabilities  are  used  as  weights. These new
counts are used to determine  new  haplotype  frequencies,
which   in   turn   are   used   to  update  the  relative
probabilities  to  new  values.  This   cyclic   iteration
continues until the likelihood is maximized (i.e., minimal
change in the lnlike).
}
\section{References}{

Excoffier, L., and Slatkin, M.,  1995,  Maximum-likelihood
estimation of molecular haplotype frequencies in a diploid
population, Mol. Biol. Evol. 12(5):921-927.

Hawley, M. E., and Kidd, K. K.,  1995,  HAPLO:  a  program
using  the  EM  algorithm  to  estimate  the frequenciesof
multi-site haplotypes, J.Heredity. 86:409-411.

Long, J. C., Williams, R. C., and Urbanek,  M.,  1995,  An
E-M  algorithm  and  testing  strategy  for multiple-locus
haplotypes, Am.J.Hum.Genet. 56:799-810.

Terwilliger, J. D., and Ott, J., 1994, Handbook  of  human
gentic linkage, Johns Hopkins University Press, Baltimore.

}
\seealso{
haplo.enum, haplo.hash, haplo.score
}
\examples{
\dontrun{
haplo <- haplo.em(geno)
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

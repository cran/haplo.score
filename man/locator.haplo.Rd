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
\name{locator.haplo}
\alias{locator.haplo}
\title{
  Find Location from Mouse Clicks and Print Haplotypes on Plot
}
\description{
Much like the Splus locator() is used to find x-y coordinates on a plot,
locator.haplo() finds all x-y coordinates that are clicked on by a
user, and then prints haplotypes at the chosen positions.
}
\usage{
locator.haplo(obj)
}
\arguments{
\item{obj}{
An object (of class haplo.score) which contains the analysis results
that are returned from the function haplo.score.
}
}
\value{
List with the following components:

\item{x.coord}{
Vector of x-coordinates.
}
\item{y.coord}{
Vector of y-coordinates.
}
\item{hap.txt}{
Vector of character strings for haplotypes.
}
}
\details{
After plotting the results in obj, as from plot(obj), the function
locator.haplo is used to place on the plot the text strings for
haplotypes of interest. After the function call (e.g.,
locator.haplo(obj)), the user can click, with the left mouse button,
on as many points in the plot as desired. Then, clicking with the
middle mouse button will cause the haplotypes to be printed on the
plot. The format of a haplotype is "a:b:c", where a, b, and c are
alleles, and the separator ":" is used to separate alleles on a
haplotype. The algorithm chooses the closest point that the user
clicks on, and prints the haplotype either above the point (for points
on the lower-half of the plot) or below the point (for points in the
upper-half of the plot).
}
\seealso{
haplo.score
}
\examples{
\dontrun{
save <-  haplo.score(y, geno, trait.type = "binomial")

plot(save)

locator.haplo(save)
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.

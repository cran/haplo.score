\name{hla.demo}
\alias{hla.demo}
\title{HLA Gene Variantion and Serologic Response to Measles Vaccination}
\description{This data set gives markers of HLA gene variantion and
  antibody levels in response to measles vaccination for 220 unrelated
  subjects, from a study by the Mayo Clinic.
}
\usage{data(hla.demo)}
\format{A data frame containing 220 rows and 10 columns:
  \describe{
    \item{resp}{quantitative antibody response to measles vaccination}
    \item{resp.cat}{a factor with levels "low", "normal", "high", for
      categorical antibody response}
    \item{male}{sex code with 1="male" ,  0="female"}
    \item{age}{age (months) at immunization}
  }
  
  The remaining variables are genotypes for 3 HLA loci, with a prefix
  name (e.g., "DQB") and a suffix for each of two alleles (".a1" and
  ".a2").
}
\source{
  Poland G, Ovsyannikova I, Jacobson R, Vierkant R, Jacobsen S,
  Pankratz V, Schaid D (2002) Identification of an association between
  HLA class II alleles and low antibody levels after measles
  immunization. Vaccine 20:430?438
}
\keyword{datasets}

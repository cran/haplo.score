# 
# Copyright 2001 Mayo Foundation for Medical Education and Research. 
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
# 02111-1307, USA.
# 
# 
haplo.score <- function(y, geno, trait.type="gaussian",
                         offset = NA, x.adj = NA, skip.haplo=.005,
                         locus.label=NA, miss.val=0, n.sim=0){

# Choose trait type:
   trait.int <- charmatch(trait.type, c("gaussian", "binomial", "poisson",
                          "ordinal"))

   if(is.na(trait.int)) stop("Invalid trait type")
   if(trait.int == 0)   stop("Ambiguous trait type")


# Check dims of y and geno
  if(length(y)!=nrow(geno)) stop("Dims of y and geno are not compatible")

  n.loci <- ncol(geno)/2

  if(n.loci != (floor(ncol(geno)/2)) )stop("Odd number of cols of geno")

# Check for missing genotypes

   miss <- apply(is.na(geno),1,any)

   if(!all(is.na(miss.val))) {
      for(mval in miss.val){
         miss <- miss | apply(geno==mval, 1, any)
      }
    }


# Check if Adjusted by Reg

  adjusted <- T

  if( all(is.na(x.adj)) ) adjusted <- F

  if(adjusted){
    x.adj <- as.matrix(x.adj)
    if(nrow(x.adj)!=length(y)) stop("Dims of y and x.adj are not compatible")
  }
        

# General checks for  missing data (more later for specific trait.type)

   miss <- miss | is.na(y) 
   if(adjusted) miss <- miss| apply(is.na(x.adj),1,any)


# If Poisson, check for errors
   if(trait.int==3) {
       if(all(is.na(offset))) stop("Missing offset")
       miss <- miss | is.na(offset)
       offset <- offset[!miss]
   }

# Subset to non-missing values:
       y <- as.numeric(y[!miss])
       geno <- geno[!miss,]
       if(adjusted) x.adj <- x.adj[!miss,,drop=F]


# If binomial, check for errors
   if(trait.int==2) {
      if(!all(y==1|y==0)) stop("Invalid y values")
      if(all(y==1) | all(y==0)) stop("No variation in y values")
    }

# If Proportional Odds, recode and check y values
   if(trait.int==4){
      y <- factor(y)
      y.lev <- levels(y)
      y <- as.numeric(y)
      if(max(y) < 3) stop("Less than 3 levels for y values")
    }



   n.subj <- length(y)

# Create a haplo object (using EMhaps)
   haplo <- haplo.em(geno, locus.label, converge.eps=0.000001, maxiter=5000) 

# Check convergence of EM
   if(!haplo$converge) stop("EM for haplo failed to converge")

# Compute haplotype score statistics

   hap1 <- haplo$hap1code
   hap2 <- haplo$hap2code
   indx <- haplo$indx.subj
   post <- haplo$post
   nreps <- as.vector(haplo$nreps)
   uhap  <- sort(unique(c(hap1,hap2)))

# Don't score haplotypes that have probabilities <  skip.haplo

   which.haplo <- haplo$hap.prob >= skip.haplo
   uhap <- uhap[which.haplo]

   x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")

   n.x <- ncol(x)
   x.post<-matrix(rep(NA, n.subj * n.x), ncol=n.x)

   for(j in 1:n.x){
      x.post[,j] <- tapply(x[,j]*post, indx, sum)
   }


# Scores for GLM's

   if(trait.int <= 3){ 

     # If not adjusted, use summaries of y

     if(!adjusted){
        mu <- switch(trait.int, mean(y), mean(y), sum(y)/sum(offset) )
        a  <- switch(trait.int, var(y), 1, 1)

        # if not adjusted, still need to set up x.adj, for intercept,
        # to be used by haplo.score.glm to compute score stat adjusted
        # for intercept

        x.adj <- matrix(rep(1,n.subj),ncol=1)
     }

    # If adjusted, use GLM to get fitted and residuals
   
     if(adjusted){

         reg.out <- glm(y ~ x.adj, family=trait.type)

         # bind col of 1's for intercept:
         x.adj <- cbind(rep(1,n.subj),x.adj)

         mu <- reg.out$fitted.values
         a  <- switch(trait.int, 
                 sum(reg.out$residuals^2)/reg.out$df.residual,
                 1, 1)
      }

     v <- switch(trait.int, 1/a, mu*(1-mu), mu )
  
     # Now compute score statistics
 
     tmp <- haplo.score.glm(y, mu, a, v, x.adj, nreps, x.post, post, x)
     u.score <- tmp$u.score
     v.score <- tmp$v.score
   }


# Scores for Proportional Odds
   if(trait.int ==4) {

      if(adjusted){
         library("Design")
         library("Hmisc")
         reg.out <- lrm(y ~ x.adj)
         K <- max(y)
         n.xadj <- ncol(x.adj)
         alpha <- reg.out$coef[1:(K-1)]
         beta <- reg.out$coeff[K:(K-1 + n.xadj)]

         tmp <- haplo.score.podds(y, alpha, beta, x.adj, nreps, x.post,
                              post, x)
       }

      if(!adjusted){
         tbl <- table(y)
         s <- 1- (cumsum(tbl)-tbl)/n.subj
         alpha <-  - log((1-s[-1])/s[-1])
         tmp <- haplo.score.podds(y, alpha, beta=NA, x.adj=NA, nreps, x.post,
                              post, x)
       }

      u.score <- tmp$u.score
      v.score <- tmp$v.score
    }
   
# Compute Score Statistics:

   tmp<-Ginv(v.score)
   df <- tmp$rank
   g.inv <- tmp$Ginv
   score.global <- u.score%*% g.inv %*%u.score
   score.haplo <- u.score / sqrt(diag(v.score))
   score.max <-  max(score.haplo^2)

# Compute empirical p-values if n.sim > 0

   if(n.sim==0){
      score.global.p.sim <- NA
      score.haplo.p.sim <- rep(NA,length(score.haplo))
      score.max.p.sim <- NA
      n.val.global <- NA
      n.val.haplo <- NA
   }

   if(n.sim > 0){

      # initialize rejection counts 
      score.global.rej <- 0
      score.haplo.rej  <- rep(0,length(score.haplo))
      score.max.rej    <- 0

      # initialize valid simulation counts 
      n.val.global <- 0
      n.val.haplo <- 0

      if(trait.int<=3){

         # Initialize mu.rand and v.rand, in case not adjusted (if adjusted, 
         # then these will be over-ridden in simulation loop

         mu.rand <- mu
         v.rand <- v  
       }


      for(i in 1:n.sim){

         # random order
         rand.ord <- order(runif(n.subj))

         if(trait.int <=3){ 
     
           if(adjusted){
              mu.rand <- mu[rand.ord]
              v.rand <- switch(trait.int, v, v[rand.ord], v[rand.ord])
            }

           tmp <- haplo.score.glm(y[rand.ord], mu.rand, a, v.rand, 
                                x.adj[rand.ord,], nreps, x.post, post, x)
         }

         if(trait.int ==4){

            if(adjusted){             
               tmp <- haplo.score.podds(y[rand.ord], alpha, beta, 
                               x.adj[rand.ord,,drop=F],nreps, x.post, post, x)
            }

            if(!adjusted) {
               tmp <- haplo.score.podds(y[rand.ord], alpha, beta=NA, 
                               x.adj=NA,nreps, x.post, post, x)
             }

          }

         u.score <- tmp$u.score
         v.score <- tmp$v.score
         
         # Now compute score statistics
     
         tmp <- Ginv(v.score)
         g.inv <- tmp$Ginv  
         score.global.sim <- u.score %*% g.inv %*% u.score

         # note that score.haplo.sim is squared, where score.haplo is not
         # because we want to keep the sign of score.haplo for returned 
         # values

         score.haplo.sim  <- (u.score / sqrt(diag(v.score)))^2

         score.max.sim <- max(score.haplo.sim)

         if(!is.na(score.global.sim)) {
            n.val.global <- n.val.global +1
            if(score.global.sim >= score.global) score.global.rej <- score.global.rej +1
          }

         if(!any(is.na(score.haplo.sim))){
            n.val.haplo <- n.val.haplo + 1

            score.haplo.rej <- score.haplo.rej +
                               ifelse(score.haplo.sim >= score.haplo^2, 1, 0)

            if(score.max.sim >= score.max) score.max.rej <- score.max.rej +1
          }

      }

      score.global.p.sim <- score.global.rej /  n.val.global
      score.haplo.p.sim <- score.haplo.rej / n.val.haplo
      score.max.p.sim <- score.max.rej / n.val.haplo

    }


   score.global.p <- 1 - pchisq(score.global,df)
   score.haplo.p <- 1-pchisq(score.haplo^2,1)

# Create locus label if missing:
   if(all(is.na(locus.label))) {
      locus.label<- paste("loc-",1:n.loci,sep="")
    }

   obj <- (list(score.global=score.global, df=df,score.global.p=score.global.p,
       score.global.p.sim=score.global.p.sim,
       score.haplo=score.haplo,score.haplo.p=score.haplo.p,
       score.haplo.p.sim=score.haplo.p.sim,
       score.max.p.sim=score.max.p.sim,
       haplotype=haplo$haplotype[which.haplo,],
       hap.prob=haplo$hap.prob[which.haplo],
       locus.label=locus.label,
       n.sim=n.sim, n.val.global=n.val.global, n.val.haplo=n.val.haplo))

   class(obj) <- "haplo.score"
   
   return(obj)

}

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
haplo.em<-function(geno, locus.label=NA, converge.eps=0.000001, maxiter=500){

# Assume no missing values in geno, and replace geno with recoded alleles of rank order
   tmp <- geno.recode(geno,miss.val=NA)
   geno <- tmp$grec

   n.loci <- ncol(geno)/2

# Create locus label if missing:
   if( all(is.na(locus.label)) ) locus.label<- paste("loc-",1:n.loci,sep="")

# Collect info for locus.info (names of loci and alleles)
   locus.info <- vector("list",n.loci)
   for(i in 1:n.loci){
     locus.info[[i]]$name    <- locus.label[i]
     locus.info[[i]]$alleles <- tmp$alist[[i]]$allele
   }

   a.freq <- vector("list",n.loci)
   genocode <- NULL
   for(i in 1:n.loci) {
      j <- (i-1)*2 + 1
      a1<- geno[,j]
      a2<- geno[,(j+1)]
      genocode <- cbind(genocode, gcode(a1,a2))
      p <- table(c(a1,a2))
      p <- p/sum(p)
      a.freq[[i]] <- list(p=p)
   }


# Create "hash" index for multilocus genotypes. This will create a
# vector "ghash", with elements that are the row number of 
# unique multilocus genotypes. The unique multilocus genotypes are in the
# matrix ugeno, so if r=ghash[i], then r is the rth row of ugeno.

   tmp <- haplo.hash(genocode)
   ghash <- tmp$hash
   ugenocode <- tmp$hap.mtx

# Uncode to get unique geno
   ugeno <- NULL
   for (i in 1:n.loci){
      ugeno <- cbind(ugeno,ungcode(ugenocode[,i]))
    }

  nr<-nrow(ugeno)

# create all possible haplotypes for each of the unique genotypes.
# the vector "indx" is a index for the different pairs of haplotypes
# that map to the same multilocus unique genotype. So, the elements
# of indx with the same value, say r, map to the same genotype, 
# and values if indx are the row number of ugeno.

   indx<-NULL
   hap1<-NULL
   hap2<-NULL

   for(i in 1:nr){
      t<-matrix(ugeno[i,],nrow=2)
      t1<-t[1,]
      t2<-t[2,]
      tmp<-haplo.enum(t1,t2)
      nhap<-nrow(tmp$h1)
      hap1<-rbind(hap1,tmp$h1)
      hap2<-rbind(hap2,tmp$h2)
      indx<-c(indx,rep(i,nhap))
   }

# nreps is the count of haplotype pairs that map to each unique
# genotype

   nreps<-tapply(indx,indx,length)

# hash code haplotypes

   tmp <- haplo.hash(rbind(hap1,hap2))
   nr  <- nrow(hap1)
   hash1 <- tmp$hash[1:nr]
   hash2 <- tmp$hash[(nr+1):(2*nr)]
   hap.mtx <- tmp$hap.mtx 

# ngeno is the counts of each unique genotype, determined from
# counting the unique genotypes in the input matrix geno.

  ngeno <- as.vector(table(ghash))

  niter <- 0
  happrob <- NULL

# lngth.ngeno.nreps is the number of unique genotypes
   lngth.ngeno.nreps<-length(ngeno)

# nhash is the number of possible pairs of haplotypes
   nhash<-length(hash1)

# for each possible pair of haplotypes, consistent with a given
# genotype, wt is the corresponding weight according to the
# genotype count.

   wt <- rep(ngeno,nreps)

# orig.hap is a matrix of unique haplotypes, with the original allele
# labels from the input geno data matrix.

   orig.hap <- locus.info[[1]]$alleles[hap.mtx[,1]]

   for(j in 2:n.loci){
      orig.hap <- cbind(orig.hap,locus.info[[j]]$alleles[hap.mtx[,j]])
   }

# init probs for genotypes, based on allele freqs, assuming no LD:

  hap.prob.noLD <- a.freq[[1]]$p[hap.mtx[,1]]

   for(j in 2:n.loci){
      hap.prob.noLD <- hap.prob.noLD * a.freq[[j]]$p[hap.mtx[,j]]
   }

   num.uhaplo <- length(hap.prob.noLD)
   prior <- hap.prob.noLD[hash1]*hap.prob.noLD[hash2]
   prior <- ifelse(hash1!=hash2, 2*prior,prior)
   ppheno <- tapply(prior,indx,sum)

   lnlike <- sum(ngeno*log(ppheno))
   lnlike.noLD <- lnlike

   den <- rep(ppheno,nreps)
   post <- prior/den

# EM steps:
   converge <- F

# Note that the em.o object code is loaded via the following commands
# at the unix prompt, when first setting up the directory with the
# source code:
#
# Splus60 CHAPTER
# Splus60 make

   em.save <- .C("em",
      eps=as.double(converge.eps),
      maxiter=as.integer(maxiter),
      niter=as.integer(niter),
      converge=as.integer(converge),
      happrob=as.double(hap.prob.noLD),
      nhash=as.integer(nhash),
      wt=as.double(wt),
      post=as.double(post),
      hap1code=as.integer(hash1),
      hap2code=as.integer(hash2),
      indx=as.integer(indx),
      ngeno=as.integer(ngeno),
      nreps=as.integer(nreps),
      lngth.ngeno.nreps=as.integer(lngth.ngeno.nreps),
      lnlike=as.double(lnlike),
      num.uhaplo=as.integer(num.uhaplo)
    )


   lr <- 2*(em.save$lnlike -lnlike.noLD)


# Merge output from em.o (which deals with unique
# genotypes and their counts,) with subject index in order to
# expand results from em.o to subject-specific information, sorted
# according to the original input geno data.

   tmp1 <- data.frame(indx.subj=(1:length(ghash)),indx.geno=ghash)

   tmp2 <- data.frame(indx.geno=em.save$indx,
            hap1code=em.save$hap1code, hap2code=em.save$hap2code,
            post=em.save$post)

   tmp3 <- merge(tmp1,tmp2,by="indx.geno")
   ord  <- order(tmp3$indx.subj)
   tmp3 <- tmp3[ord,]

# create new nreps for expansion to subject-specific info:
   nreps <- tapply(tmp3$indx.subj,tmp3$indx.subj,length)

   return(list(lnlike=em.save$lnlike,lr=lr,hap.prob=em.save$happrob,
             hap.prob.noLD=hap.prob.noLD,converge=em.save$converge,
             locus.info=locus.info, locus.label=locus.label,
             indx.subj=tmp3$indx.subj,post=tmp3$post,
             hap1code=tmp3$hap1code,hap2code=tmp3$hap2code,
             haplotype=orig.hap,nreps=nreps,niter=em.save$niter))
}











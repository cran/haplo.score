# Test Splus code, and steps used to create User Manual


# Attach haplo.score routines and hla.demo, by use of a local libary
# path.
# 
# BE SURE TO CHANGE "INSTALL DIR PATH HERE" TO THE DIRECTORY PATH WHERE
# haplo.score IS INSTALLED.

library(haplo.score)
data(hla.demo)

names(hla.demo)

attach(hla.demo)

# Dataframe matrix of columns of alleles for pairs of loci:
geno <- hla.demo[,-c(1:4)]
locus.label <-c("DQB","DRB","B")

# Haplotype analyses for quantitative trait, resp

hap.gaus <- haplo.score(resp, geno, trait.type="gaussian",
                        offset = NA, x.adj = NA, skip.haplo=.005,
                        locus.label=locus.label, miss.val=0, n.sim=0)

print(hap.gaus)

#motif()
#plot(hap.gaus)
#title("Figure 1. Haplotype Score Statistics\nQuantitative Response",cex=.6)
#loc.gaus <- locator.haplo(hap.gaus)
#dev.off()

# now create a good eps file:

#postscript(file="Fig1.eps",onefile=T, width=5, height=4, print.it=F, horizo=F)
#plot(hap.gaus,cex=.6)
#title("Figure 1. Haplotype Score Statistics\nQuantitative Response",cex=.6)
#text(loc.gaus$x.coord,loc.gaus$y.coord,loc.gaus$hap.txt,cex=.5)
#dev.off()


# analyses with skip.haplo = .01
hap.gaus2 <- haplo.score(resp, geno, trait.type="gaussian",
                         offset = NA, x.adj = NA, skip.haplo=.01,
                         locus.label=locus.label, miss.val=0, n.sim=0)
print(hap.gaus2)


# Now try adjusted scores, adjusted for sex and age.

x.ma <- cbind(male, age) # set-up design matrix, with male=1/0 (male/female) and age in months

hap.gaus.adj <- haplo.score(resp, geno, trait.type="gaussian",
                        offset = NA, x.adj = x.ma, skip.haplo=.005,
                        locus.label=locus.label, miss.val=0, n.sim=0)
print(hap.gaus.adj)


# Haplotype analyses for ordinal trait

y.ord <- as.numeric(resp.cat)

hap.ord <- haplo.score(y.ord, geno, trait.type="ordinal",
                        offset = NA, x.adj = NA, skip.haplo=.005,
                        locus.label=locus.label, miss.val=0, n.sim=0)
print(hap.ord)

#motif()
#plot(hap.ord)
#title("Haplotype Score Statistics\n Ordinal Response")
#loc.ord <- locator.haplo(hap.ord)
#dev.off()

# now create a good eps file:

#postscript(file="Fig2.eps",onefile=T, width=5, height=4, print.it=F, horizo=F)
#plot(hap.ord,cex=.6)
#title("Figure 2. Haplotype Score Statistics\nOrdinal Response",cex=.6)
#text(loc.ord$x.coord,loc.ord$y.coord,loc.ord$hap.txt,cex=.5)
#dev.off()

# Haplotype analyses for binary trait, y.bin

y.bin <-ifelse(y.ord==1,1,0)

hap.bin <- haplo.score(y.bin, geno, trait.type="binomial",
                        offset = NA, x.adj = NA, skip.haplo=.005,
                        locus.label=locus.label, miss.val=0, n.sim=0)
print(hap.bin)

#motif()
#plot(hap.bin)
#title("Haplotype Score Statistics\nBinomial Response")
#loc.bin <- locator.haplo(hap.bin)
#dev.off()

# now create a good eps file:

#postscript(file="Fig3.eps",onefile=T, width=5, height=4, print.it=F, horizo=F)
#plot(hap.bin,cex=.6)
#title("Figure 3. Haplotype Score Statistics\nBinary Response",cex=.6)
#text(loc.bin$x.coord,loc.bin$y.coord,loc.bin$hap.txt,cex=.5)
#dev.off()


# demo of simulation p-values


hap.bin.sim <- haplo.score(y.bin, geno, trait.type="binomial",
                        offset = NA, x.adj = NA, skip.haplo=.005,
                        locus.label=locus.label, miss.val=0, n.sim=1000)
print(hap.bin.sim)

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
print.haplo.score <- function(x, ...){

# print of global score stats:
   banner("Global Score Statistics", banner.width=80, char.perline=60,
                border= "-")
   cat(paste("global-stat = ",round(x$score.global,5),", df = ",x$df,
         ", p-val = ",round(x$score.global.p,5),sep=""))

   if(x$n.sim>0) cat(", sim. p-val = ",x$score.global.p.sim,"\n\n")

   if(x$n.sim>0) cat("max-stat sim. p-val = ",x$score.max.p.sim)
   cat("\n\n")

# create table for haplotype specific stats:

   tbl <- cbind(x$haplotype,round(x$hap.prob,5),round(x$score.haplo,5),
          round(x$score.haplo.p,5))

   if(x$n.sim>0) tbl <- cbind(tbl,x$score.haplo.p.sim)

   ord <- order(x$score.haplo)
   tbl <- tbl[ord,]


   if(x$n.sim == 0) dimnames(tbl) <- list(NULL,c(x$locus.label,"Hap-Freq",
                  "Hap-Score","p-val"))

   if(x$n.sim > 0) dimnames(tbl) <- list(NULL,c(x$locus.label,"Hap-Freq",
                  "Hap-Score","p-val","sim p-val"))

   banner("Haplotype-specific Scores", banner.width=80, char.perline=60,
                border= "-")
   print(tbl,quote=F)

   cat("\n\n")

   invisible()

}

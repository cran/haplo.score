/*

Copyright 2001 Mayo Foundation for Medical Education and Research. 
	.
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
02111-1307, USA.

*/


#include <math.h>
#include <S.h>

void em(double *eps,
        int *maxiter,
        int *niter,
        int *converge,
        double *happrob,
        int *nhash,
        double *wt,
        double *post,
        int *hap1code,
        int *hap2code,
        int *indx,
        int *ngeno,
        int *nreps,
        int *length_of_ngeno_nreps,
        double *lnlike, 
        int *number_of_unique_haplotypes){
 
  double *expecthaps,*hapsum,hapsumtotal,*probpheno,*prior,tmplnlike;

  int i,j,k,index;
 
  expecthaps= (double*) R_alloc(*nhash,sizeof(double));
  hapsum= (double*) R_alloc(*number_of_unique_haplotypes,sizeof(double));
  prior= (double*) R_alloc(*nhash,sizeof(double));
  probpheno= (double*) R_alloc(*length_of_ngeno_nreps,sizeof(double));

  for (k=2;k<=*maxiter;k++){

    /********************************************
     *** expecthaps<-c(post,post)*c(wt,wt)
     *** hap.sum<-tapply(expecthaps,c(hap1code,hap2code),sum
     ********************************************/

    for (i=0;i<*number_of_unique_haplotypes;i++){
      hapsum[i]=0.0;
    }
  
    hapsumtotal=0.0;
    for (i=0;i<*nhash;i++) {
      expecthaps[i]=post[i]*wt[i];
      hapsum[hap1code[i]-1]+=expecthaps[i];
      hapsum[hap2code[i]-1]+=expecthaps[i];
      hapsumtotal+=2.0*expecthaps[i];
    } 

    /********************************************
     *** hap.prob<-hap.sum/sum(hap.sum)
    ********************************************/

    for (i=0;i<*number_of_unique_haplotypes;i++){
      happrob[i]=hapsum[i]/(hapsumtotal);
    }

    /********************************************
     *** prior<-hap.prob[hap1code]*hap.prob[hap2code]
     *** prior<-ifelse(hap1code != hap2code,2 * prior, prior)
    ********************************************/

    for (i=0;i<*nhash;i++) {
      prior[i]=happrob[hap1code[i]-1]*happrob[hap2code[i]-1];
      if (hap1code[i]!=hap2code[i])
        prior[i]=2*prior[i];
    }

    /********************************************
     ***probpheno<-tapply(prior,indx,sum)
     *******************************************/

    for (i=0;i<*length_of_ngeno_nreps;i++)
      probpheno[i]=0.0;

    for (i=0;i<*nhash;i++){
      probpheno[indx[i]-1]+=prior[i];
    }

    /******************************************
     ***lnlike<-sum(ngeno*log(probpheno))
     ******************************************/

    tmplnlike=0.0;
    for (i=0;i<*length_of_ngeno_nreps;i++) {
      tmplnlike+=ngeno[i]*log(probpheno[i]);
    }

    /******************************************
    *** den<-rep(probpheno,nreps)
    *** post<-prior/den
    ******************************************/

    index=0;
    for (i=0;i<*length_of_ngeno_nreps;i++) {
      for (j=0;j<nreps[i];j++) {
        post[index]=prior[index]/probpheno[i];
        index++;
      }
    }

    if (fabs(*lnlike-tmplnlike)<*eps) {
      *converge=1;
       break;
     }

   *lnlike=tmplnlike;
  }

 *niter=k;

 return;

}






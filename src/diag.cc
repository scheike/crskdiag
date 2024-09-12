#include <iostream>
#include <math.h>
#include <Rmath.h>
#include "matrix.h"
#include "util.h"
//#include <vector>

using namespace std;

extern "C" {void diag_lin(double *times, int *cause, double *designX, int *N, int *M, 
					int *Nit, int *NJP, double *TJP, double *betaS, double *beta_var, double *dlamb0,
					double *uniX, int *TC, double *B, double *MAV,  int *nsim, int *nplot, double *pval, int *minor_included);}
					
extern "C" {void diag_prop(double *times, int *cause, double *designX, int *N, int *M, 
					int *Nit, int *NJP, double *TJP, double *betaS, double *beta_var, double *dlamb0,
					double *Bt, double *MAVt,  int *nsim, int *nplot, double *pvalt, int *minor_included);}					

void diag_prop(double *times, int *cause, double *designX, int *N, int *M, 
		int *Nit, int *NJP, double *TJP, double *betaS, double *beta_var, double *dlamb0,
		double *Bt, double *MAVt,  int *nsim, int *nplot, double *pvalt, int *minor_included){
	
	double tmp1;	
	//double gn[*N];
	::vector *gn;
	::vector *sqSI,*E,*tmpv1,*tmpv2,*xi,*rowX; //M; 
	::vector *maxB1t,*maxB2t,*rrt; //M+1
	::vector *z,*Gbeta;   //N
	::vector *S0j,*absB1t,*absB2t,*absB1t_all,*absB2t_all;   //N+1
	matrix *SI;  //M*M*/
	matrix *WX,*Wbeta,*eta;  //N*M
	matrix *Ej,*Uj,*tmpB1t,*tmpB2t,*Ct;     //(N+1)*M
	//::vector *wy[*N],*wdn1[*N],*wdm1[*N];  //(N+1)*N
	matrix *wy,*wdn1,*wdm1;  //(N+1)*N
	//::vector *eta[*N];   //M*N
	//::vector *Ct[*N+1];
	//matrix *W1t[*N+1],*W2t[*N+1];
	matrix3 *W1t,*W2t;   //(N+1)*N*M
	
	malloc_vecs(*M,&sqSI,&E,&tmpv1,&tmpv2,&xi,&rowX,NULL);
	malloc_vecs(*M+1,&maxB1t,&maxB2t,&rrt,NULL);
	malloc_vecs(*N,&z,&Gbeta,&gn,NULL);
	malloc_vecs(*N+1,&S0j,&absB1t,&absB2t,&absB1t_all,&absB2t_all,NULL);
	malloc_mats(*M,*M,&SI,NULL);
	malloc_mats(*N,*M,&WX,&Wbeta,&eta,NULL);
	malloc_mats(*N+1,*M,&Ej,&Uj,&tmpB1t,&tmpB2t,&Ct,NULL);
	//for(int i=0;i<*N;i++){
		//malloc_vecs(*N+1,&wy[i],&wdn1[i],&wdm1[i],NULL);
	//malloc_vecs(*M,&eta[i],NULL);}
	malloc_mats(*N,*N+1,&wy,&wdn1,&wdm1,NULL);
	for(int j=0;j<=*N;j++){
		//malloc_vecs(*M,&Ct[j],NULL);
		//malloc_mats(*N,*M,&W1t[j],&W2t[j],NULL);
	}
	malloc_mat3(*N+1,*N,*M,W1t);
	malloc_mat3(*N+1,*N,*M,W2t);
	
	GetRNGstate();
	
	for(int k=0;k<*M;k++){
		for (int i=0;i<*N;i++){
			ME(WX,i,k)=designX[k*(*N)+i];  //read design matrix into a matrix, read by column
		}
	}	
	est_tmp(times, cause, WX, N, M, Nit, NJP, TJP, betaS, beta_var, dlamb0, Gbeta, eta, wy, wdn1, wdm1, S0j, Ej, Uj, SI, minor_included);
	
	for(int k=0;k<*M;k++) VE(sqSI,k) = pow(ME(SI,k,k),0.5);
	for(int i=0;i<*N;i++){
		extract_row(eta,i,tmpv1);
		vM(SI,tmpv1,tmpv1);
		replace_row(Wbeta,i,tmpv1);
	}
	
	for(int k=0;k<*M;k++){		
		extract_col(WX,k,z); 
		vec_zeros(absB1t);
		for(int j=1;j<=*NJP;j++){
			tmp1 = 0;
			for(int i=0;i<*N;i++){
				tmp1 += VE(z,i)*ME(wdm1,i,j);
			}
			ME(tmpB1t,j,k) = ME(tmpB1t,j-1,k)+tmp1;
			VE(absB1t,j) = fabs(ME(tmpB1t,j,k))*VE(sqSI,k);   //replicate score process 
			Bt[(k)*(*N+1)+j] = ME(tmpB1t,j,k)*VE(sqSI,k);
		}
		VE(maxB1t,k) = vec_max(absB1t);
		MAVt[k]=VE(maxB1t,k);
	}
	
	for(int j=1;j<=*NJP;j++){
		for(int k=0;k<*M;k++){	
			VE(absB1t_all,j) += fabs(ME(tmpB1t,j,k))*VE(sqSI,k);			
		}
	}
	VE(maxB1t,*M) = vec_max(absB1t_all);
	MAVt[*M]=VE(maxB1t,*M);
	
	for(int k=0;k<*M;k++){		
		extract_col(WX,k,z); 
		for(int j=1;j<=*NJP;j++){
			extract_row(Ej,j,E);   		
			for(int i=0;i<*N;i++){
				ME3(W1t,j,i,k) = ME3(W1t,j-1,i,k)+(VE(z,i)-VE(E,k))*ME(wdm1,i,j);
			}
		}		
	
		for(int j=1;j<=*NJP;j++){
			extract_row(Ej,j,E);
			vec_zeros(tmpv2);
			for(int i=0;i<*N;i++){
				extract_row(WX,i,xi);
				vec_subtr(xi,E,rowX);
				tmp1 = VE(z,i)*ME(wy,i,j)*exp(VE(Gbeta,i))*dlamb0[j];
				scl_vec_mult(tmp1,rowX,tmpv1);
				vec_add(tmpv2,tmpv1,tmpv2);   //sum of i
			}
			extract_row(Ct,j-1,tmpv1);
			vec_subtr(tmpv1,tmpv2,tmpv2);			// sum up to t
			replace_row(Ct,j,tmpv2);
			
		}
		
		for(int i=0;i<*N;i++){
			for(int j=1;j<=*NJP;j++){
				extract_row(Wbeta,i,tmpv1);
				extract_row(Ct,j,tmpv2);
				ME3(W2t,j,i,k)=vec_prod(tmpv2,tmpv1);
			}
		}
	}		
	
	for(int r=0; r<*nsim; r++){
		mat_zeros(tmpB2t);	
		for(int k=0;k<*M;k++){	
			vec_zeros(absB2t);
			for(int i=0;i<*N;i++) VE(gn,i) = rnorm(0,1); 	
			for(int j=1;j<=*NJP;j++){
				for(int i=0;i<*N;i++){
					ME(tmpB2t,j,k) += (ME3(W1t,j,i,k)+ME3(W2t,j,i,k))*VE(gn,i);
				}
				VE(absB2t,j) = fabs(ME(tmpB2t,j,k))*VE(sqSI,k);
				if((r>0)&&r<=(*nplot)) Bt[r*(*M)*(*N+1)+(k)*(*N+1)+j] = ME(tmpB2t,j,k)*VE(sqSI,k);  //return to R for plot
			}
			VE(maxB2t,k) = vec_max(absB2t);
			if(VE(maxB2t,k)>=VE(maxB1t,k))  VE(rrt,k) += 1;
		}

		vec_zeros(absB2t_all);
		for(int j=0;j<=*NJP;j++){
			for(int k=0;k<*M;k++){	
				VE(absB2t_all,j) += fabs(ME(tmpB2t,j,k))*VE(sqSI,k);
			}
		}
		VE(maxB2t,*M) = vec_max(absB2t_all);	
		if(VE(maxB2t,*M)>=VE(maxB1t,*M)) VE(rrt,*M) += 1;	
	}
	for(int k=0;k<=*M;k++) pvalt[k] = VE(rrt,k)/(*nsim);
	
	PutRNGstate();
		
	free_vecs(&sqSI,&E,&tmpv1,&tmpv2,&xi,&rowX,NULL);
	free_vecs(&maxB1t,&maxB2t,&rrt,NULL);	
	free_vecs(&z,&Gbeta,&gn,NULL);
	free_vecs(&S0j,&absB1t,&absB2t,&absB1t_all,&absB2t_all,NULL);
	free_mats(&SI,NULL);
	free_mats(&WX,&Wbeta,&eta,NULL);
	free_mats(&Ej,&Uj,&tmpB1t,&tmpB2t,&Ct,NULL);
	//for(int i=0;i<*N;i++){
		//free_vecs(&wy[i],&wdn1[i],&wdm1[i],NULL);
		//free_vecs(&eta[i],NULL);
	//}
	free_mats(&wy,&wdn1,&wdm1,NULL);
	//for(int j=0;j<=*N;j++){
		//free_vecs(&Ct[j],NULL);
		//free_mats(&W1t[j],&W2t[j],NULL);
	//}
	free_mat3(W1t);
	free_mat3(W2t);
	
}			
					
					
void diag_lin(double *times, int *cause, double *designX, int *N, int *M, 
		int *Nit, int *NJP, double *TJP, double *betaS, double *beta_var, double *dlamb0,
		double *uniX, int *TC, double *B, double *MAV,  int *nsim, int *nplot, double *pval, int *minor_included){
	
	double tmp1;
	//double gn[*N];
	//double *gn;
	//::vector<double> gn(*N);	
	::vector *gn;
	::vector *tmpv1,*tmpv2,*xi,*E,*rowX; //M
	::vector *maxB1,*maxB2,*rr;  //M+1
	::vector *z,*Gbeta,*absB1,*absB2; //N
	::vector *S0j; //N+1
	matrix *SI;  //M*M
	matrix *WX,*eta,*C;  //N*M
	matrix *tmpB1a,*W1,*tmpB1,*tmpB2,*tmpB2a,*tmpB2b; //N*(M+1) 
	matrix *g;   //N*(N+1)
	matrix *Ej,*Uj;     //(N+1)*M
	//::vector *ind[*M+1],*zst[*M+1],*uniz[*M+1];
	::vector *ind,*zst,*uniz;   //N
	//::vector *wy[*N],*wdn1[*N],*wdm1[*N];  //(N+1)*N
	matrix *wy,*wdn1,*wdm1;  //(N+1)*N
	//::vector *eta[*N],*C[*N];   //M*N
	//matrix *W[*N],*W2[*N],*W3[*N]; //N*(M+1)*N
	//matrix *W2[*N],*W3[*N]; //N*(M+1)*N
	matrix3 *W2,*W3;  //N*(M+1)*N matrix3
	
	malloc_vecs(*M,&maxB1,&tmpv1,&tmpv2,&xi,&E,&rowX,&maxB2,NULL);
	malloc_vecs(*M+1,&maxB1,&maxB2,&rr,NULL);
	malloc_vecs(*N,&z,&Gbeta,&uniz,&absB1,&absB2,&gn,NULL);
	malloc_vecs(*N+1,&S0j,NULL);
	malloc_mats(*M,*M,&SI,NULL);
	malloc_mats(*N,*M,&WX,&eta,&C,NULL);
	malloc_mats(*N,*M+1,&tmpB1a,&W1,&tmpB1,&tmpB2,&tmpB2a,&tmpB2b,NULL);
	malloc_mats(*N,*N+1,&g,NULL);
	malloc_mats(*N+1,*M,&Ej,&Uj,NULL);
	// for(int k=0;k<=*M;k++){
		// malloc_vecs(*N,&ind[k],&zst[k],&uniz[k],NULL);
	// }
	//malloc_mats(*N,*M,&ind,&zst,&uniz,NULL);
	malloc_vecs(*N,&ind,&zst,&uniz,NULL);
	//for(int i=0;i<*N;i++){
		//malloc_vecs(*N+1,&wy[i],&wdn1[i],&wdm1[i],NULL);
		//malloc_vecs(*M,&eta[i],&C[i],NULL);
		//malloc_mats(*N,*M+1,&W[i],&W2[i],&W3[i],NULL);
		//malloc_mats(*N,*M+1,&W2[i],&W3[i],NULL);
	//}
	malloc_mats(*N,*N+1,&wy,&wdn1,&wdm1,NULL);
	malloc_mat3(*N,*N,*M+1,W2);
	malloc_mat3(*N,*N,*M+1,W3);
	
	GetRNGstate();

	for(int k=0;k<*M;k++){
		for (int i=0;i<*N;i++){
			ME(WX,i,k)=designX[k*(*N)+i];  //read design matrix into a matrix, read by column
		}
	}	
	est_tmp(times, cause, WX, N, M, Nit, NJP, TJP, betaS, beta_var, dlamb0, Gbeta, eta, wy, wdn1, wdm1, S0j, Ej, Uj, SI, minor_included);

	for(int k=0;k<=*M;k++){		
		if(k<*M) extract_col(WX,k,z); 
		if(k==*M) vec_copy(Gbeta,z);
		indsort(z,zst,ind);	
		TC[k] = unisort(z,uniz);
		for(int c=0;c<TC[k];c++) uniX[(k)*(*N)+c]=VE(uniz,c);    //return to R for plot		
				
		vec_zeros(absB1);
		for(int i=0;i<*N;i++){ 
			for(int j=1;j<=*NJP;j++){
				ME(W1,i,k) += ME(wdm1,int(VE(ind,i)),j);
			}
			if(i==0) ME(tmpB1a,i,k) = ME(W1,i,k);
			else ME(tmpB1a,i,k) = ME(tmpB1a,i-1,k)+ME(W1,i,k);
		}
		
		for(int c=0;c<TC[k];c++){
			for(int i=0;i<*N;i++){
				if(VE(zst,i)==VE(uniz,c)) ME(tmpB1,c,k) = ME(tmpB1a,i,k);
				else if(VE(zst,i)>VE(uniz,c)) break;
			}			
		}
						
		for(int c=0;c<TC[k];c++){
			VE(absB1,c) = fabs(ME(tmpB1,c,k));
			B[(k)*(*N)+c] = ME(tmpB1,c,k);  //observed value, return to R for plot	
		}			
		VE(maxB1,k) = vec_max(absB1);
		MAV[k] = VE(maxB1,k);   //Max abs value(observed), return to R
		
		for(int j=1;j<=*NJP;j++){
			for(int c=0;c<*N;c++){
				if(c==0) ME(g,c,j)=-ME(wy,int(VE(ind,c)),j)*exp(VE(Gbeta,int(VE(ind,c))))/VE(S0j,j);
				else ME(g,c,j) = ME(g,c-1,j)-ME(wy,int(VE(ind,c)),j)*exp(VE(Gbeta,int(VE(ind,c))))/VE(S0j,j);	//add c
			}
		}
		
		for(int c=0;c<*N;c++){
			for(int i=0;i<*N;i++){
				for(int j=1;j<=*NJP;j++){	
					//ME(W2[c],i,k) += ME(g,c,j)*VE(wdm1[i],j);   //add j
					ME3(W2,c,i,k) += ME(g,c,j)*ME(wdm1,i,j);   //add j
				}
			}
		}	
		
		for(int c=0;c<*N;c++){  
			extract_row(WX,int(VE(ind,c)),xi);
			vec_zeros(tmpv2);
			for(int j=1;j<=*NJP;j++){
				extract_row(Ej,j,E);				
				tmp1 = ME(wy,int(VE(ind,c)),j)*exp(VE(Gbeta,int(VE(ind,c))))*dlamb0[j];
				vec_subtr(xi,E,rowX);
				scl_vec_mult(tmp1,rowX,tmpv1);
				vec_add(tmpv2,tmpv1,tmpv2);          //add j
			}
			// if(c==0) scl_vec_mult(-1,tmpv2,C[c]);    //add c
			// else vec_subtr(C[c-1],tmpv2,C[c]);	
			
			if(c==0) scl_vec_mult(-1,tmpv2,tmpv1);    //add c
			else{
				extract_row(C,c-1,tmpv1);
				vec_subtr(tmpv1,tmpv2,tmpv1);
			}	
			replace_row(C,c,tmpv1);
		}	
		
		for(int c=0;c<*N;c++){
			for(int i=0;i<*N;i++){
				extract_row(C,c,tmpv1);
				vM(SI,tmpv1,tmpv2);
				//vec_copy(eta[i],tmpv1);
				extract_row(eta,i,tmpv1);
				ME3(W3,c,i,k) = vec_prod(tmpv2,tmpv1);
			}	
		}	
	} //end of k
	
	for(int r=0; r<*nsim; r++){
		mat_zeros(tmpB2);
		mat_zeros(tmpB2a);
		mat_zeros(tmpB2b);			
		for(int k=0;k<=*M;k++){	
			vec_zeros(absB2);
			if(k<*M) extract_col(WX,k,z); 
			else if(k==*M) vec_copy(Gbeta,z);
			TC[k] = unisort(z,uniz);
			indsort(z,zst,ind);	
			for(int i=0;i<*N;i++) VE(gn,i) = rnorm(0.0,1.0);
			for(int c=0;c<*N;c++){
				for(int i=0;i<*N;i++){
					//ME(tmpB2b,c,k) += ME(W2[c],i,k)*gn[i]+ME(W3[c],i,k)*gn[i];
					ME(tmpB2b,c,k) += ME3(W2,c,i,k)*VE(gn,i)+ME3(W3,c,i,k)*VE(gn,i);
				}
			}
				
			for(int i=0;i<*N;i++){
				if(i==0) ME(tmpB2a,i,k) = ME(W1,i,k)*VE(gn,int(VE(ind,i)));
				else ME(tmpB2a,i,k) = ME(tmpB2a,i-1,k)+ME(W1,i,k)*VE(gn,int(VE(ind,i)));
			}
			
			for(int c=0;c<TC[k];c++){
				for(int i=0;i<*N;i++){
					if(VE(zst,i)==VE(uniz,c)) ME(tmpB2,c,k) = ME(tmpB2a,i,k)+ME(tmpB2b,i,k);
					else if(VE(zst,i)>VE(uniz,c)) break;
				}			
			}
			
			for(int c=0;c<TC[k];c++){
				VE(absB2,c) = fabs(ME(tmpB2,c,k));
				if((r>0)&&r<=(*nplot)) B[r*(*M+1)*(*N)+(k)*(*N)+c] = ME(tmpB2,c,k);  //return to R for plot
			}
			VE(maxB2,k) = vec_max(absB2);
			if(VE(maxB2,k)>=VE(maxB1,k))  VE(rr,k) += 1;			
		}  //end of k
	}  //end of r
	
	for(int k=0;k<=*M;k++) pval[k] = VE(rr,k)/(*nsim);
	
	PutRNGstate();
	
	free_vecs(&tmpv1,&tmpv2,&xi,&E,&rowX,NULL);
	free_vecs(&maxB1,&maxB2,&rr,NULL);
	free_vecs(&z,&Gbeta,&absB1,&absB2,&gn,NULL);
	free_vecs(&S0j,NULL);
	free_mats(&SI,NULL);
	free_mats(&WX,&eta,&C,NULL);
	free_mats(&tmpB1a,&W1,&tmpB1,&tmpB2,&tmpB2a,&tmpB2b,NULL);
	free_mats(&g,NULL);
	free_mats(&Ej,&Uj,NULL);
	// for(int k=0;k<=*M;k++){
		// free_vecs(&ind[k],&zst[k],&uniz[k],NULL);
	// }
	free_vecs(&ind,&zst,&uniz,NULL);
	for(int i=0;i<*N;i++){
		//free_vecs(&wy[i],&wdn1[i],&wdm1[i],NULL);
		//free_vecs(&eta[i],&C[i],NULL);
		//free_mats(&W[i],&W2[i],&W3[i],NULL);
		//free_mats(&W2[i],&W3[i],NULL);
	}
	free_mats(&wy,&wdn1,&wdm1,NULL);
	free_mat3(W2);
	free_mat3(W3);
}					
					

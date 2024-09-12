//util.cc
#include <iostream>
#include <math.h>
#include <Rmath.h>
#include "util.h"

void est_tmp(double *times, int *cause, matrix *WX, int *N, int *M, 
		int *Nit, int *NJP, double *TJP, double *betaS, double *beta_var, double *dlamb0, ::vector *Gbeta, 
		matrix *eta, matrix *wy, matrix *wdn1, matrix *wdm1, ::vector *S0j, matrix *Ej, matrix *Uj, matrix *SI, int *minor_included){
					
	int njp;
	double S0,ebz,ss,tmp,dnci,yi;
	::vector *beta,*U,*S1,*dS0,*xi,*E,*tmpv1,*delta,*sqSI,*rowX,*q;  //M
	::vector *Idxt;  //N
	::vector *survc,*dnc,*dlambc,*yl;   //N+1
	matrix *I,*S2,*tmpm1,*tmpm2,*tmpm3,*Sigma; //M*M 
	matrix *qj;
	matrix *minor;  //M*N
	matrix *dmc;	//N*(N+1)
	// ::vector *minor[*N]; //M*N
	// ::vector *dmc[*N]; //N length of N+1 ::vectors

	
	malloc_vecs(*M,&beta,&U,&S1,&dS0,&xi,&E,&tmpv1,&delta,&sqSI,&rowX,&q,NULL);
	malloc_vecs(*N,&Idxt,NULL);
	malloc_vecs(*N+1,&survc,&dnc,&dlambc,&yl,NULL);
	malloc_mats(*M,*M,&I,&S2,&tmpm1,&tmpm2,&tmpm3,&Sigma,NULL);   //
	malloc_mats(*N+1,*M,&qj,NULL);
	malloc_mats(*N,*M,&minor,NULL);
	malloc_mats(*N,*N+1,&dmc,NULL);
	// for(int i=0;i<*N;i++){
		// malloc_vecs(*M,&minor[i],NULL);
		// malloc_vecs(*N+1,&dmc[i],NULL);
	// }

	//i/s=0,N-1
	//j=1,NJP, all related to j defined length as N+1, just in case of N=NJP
	//k,l=0,M-1
	//it=0,Nit-1
	//c=1,NC
	
	// // // /******************
	// // // Difference btw COXmodel and FGmodel
	// // // 1.riskset coding
	// // // 2.variable name:dn1=>wdn1, dm1=>wdm1, 
	// // // 3.minor term R
	// // // ****************************
	for(int k=0;k<*M;k++){
	    VE(beta,k)=betaS[k];
	}
	
	njp = 0;
	TJP[0] = 0.0;
	for(int i=0; i<*N; i++){                  
		if((times[i]>TJP[njp])&&(cause[i]==1)){           //only consider cause 1 jumping time
			njp += 1;
			TJP[njp] = times[i];
		}
		VE(Idxt,i) = njp;
	}
	*NJP = njp;

	//making dmc[i,j], survc[j], yl[j]

	VE(survc,0)=1.0;
	VE(yl,0)=*N;
	for(int j=1; j<=*NJP; j++){
		for(int i=0; i<*N; i++){
			if(times[i]==TJP[j] && cause[i]==0) VE(dnc,j) += 1.0;    
			if(times[i]>=TJP[j]) VE(yl,j) += 1.0; 
		}
		if(VE(yl,j)>0) VE(dlambc,j)=VE(dnc,j)/VE(yl,j);
		else VE(dlambc,j)=0;
		VE(survc,j)=VE(survc,j-1)*(1.0-VE(dlambc,j));
		for(int i=0; i<*N; i++){
			dnci=0.;
			yi=0.;
			if(times[i]==TJP[j] && cause[i]==0) dnci = 1.0;    
			if(times[i]>=TJP[j]) yi = 1.0; 
			ME(dmc,i,j) = dnci - yi*VE(dlambc,j);
		}
	}
	
	
	///////wy[i,j] wdn1[i,j]
	int cnt=0;
	for(int i=0; i<*N; i++){
		for(int j=1; j<=*NJP; j++){
			if(times[i]>=TJP[j]){
				ME(wy,i,j)=1;
				if((times[i]==TJP[j])&&(cause[i]==1)){
					ME(wdn1,i,j)=1;
					cnt += 1;
				}
			}
			else if(times[i]<TJP[j] && cause[i]>1){
				ME(wy,i,j)=VE(survc,j)/VE(survc,(int)VE(Idxt,i));
			}
		}
	}
	
	// //Newton Raphson for beta
	for(int it=0;it<*Nit;it++){
		Mv(WX,beta,Gbeta); 
		mat_zeros(I);
		vec_zeros(U);
		ss = 0;
		for(int j=1; j<=*NJP; j++){
			S0=0.;   
			vec_zeros(dS0);
			vec_zeros(S1);
			mat_zeros(S2);
			for(int i=0; i<*N; i++){
				ebz = exp(VE(Gbeta,i));
				S0 += ME(wy,i,j)*ebz;
				extract_row(WX,i,xi);	
				scl_vec_mult(ME(wy,i,j)*ebz,xi,tmpv1);
				vec_add(S1,tmpv1,S1);
				vec_copy(S1,dS0);
				vtv(xi,xi,tmpm1);
				scl_mat_mult(ME(wy,i,j)*ebz,tmpm1,tmpm1);
				mat_add(S2,tmpm1,S2);	
			}	
			
			for(int i=0; i<*N; i++){
				extract_row(WX,i,xi); 
				for(int k=0; k<*M; k++){
					VE(E,k) = VE(S1,k)/S0;
					VE(U,k) += (VE(xi,k)-VE(E,k))*ME(wdn1,i,j);
					for(int l=0; l<*M; l++){
						ME(tmpm3,k,l) = ME(S2,k,l)/S0;
						ME(tmpm2,k,l) = VE(E,k)*VE(E,l);
						ME(I,k,l) += (ME(tmpm3,k,l)-ME(tmpm2,k,l))*ME(wdn1,i,j);
					}
				}
			}
			
			if(it==*Nit-1){
				VE(S0j,j)=S0;
				replace_row(Ej,j,E);
				replace_row(Uj,j,U);			
			}
		} //end of j

		invertS(I,SI,1); 
		Mv(SI,U,delta);
		vec_add(beta,delta,beta);		
		for (int k=0;k<*M;k++) ss += fabs(VE(U,k)); 
		if((fabs(ss)<0.0001)&&(it<*Nit-2)) it=*Nit-2;
	} //end of it	

	for(int k=0;k<*M;k++){	          //return to R	
		betaS[k] = VE(beta,k);	
		VE(sqSI,k) = pow(ME(SI,k,k),0.5);
	}
	
	for(int j=1;j<=*NJP;j++){
		tmp=0;
		for(int i=0;i<*N;i++){
			tmp += ME(wdn1,i,j);
		}
		dlamb0[j]=tmp/VE(S0j,j);
	}
	
	for(int i=0;i<*N;i++){
		for(int j=1;j<=*NJP;j++){
			ME(wdm1,i,j)=ME(wdn1,i,j)-ME(wy,i,j)*exp(VE(Gbeta,i))*dlamb0[j];
		}	
	}	

	for(int k=0;k<*M;k++){
		for(int i=0;i<*N;i++){
			for(int j=1;j<=*NJP;j++){     
				ME(eta,i,k) += (ME(WX,i,k)-ME(Ej,j,k))*ME(wdm1,i,j);
			}
		}
	}
	
	// ////update eta with the minor term
	if(minor_included[0]==1){
		for(int j=1;j<=*NJP;j++){
			vec_zeros(q);
			for(int t=j;t<=*NJP;t++){
				extract_row(Ej,t,E); 
				for(int i=0;i<*N;i++){
					extract_row(WX,i,xi); 
					if(times[i]<=TJP[j]){                  
						vec_subtr(xi,E,rowX); 
						scl_vec_mult(ME(wdn1,i,j),rowX,tmpv1);
						vec_subtr(q,tmpv1,q);            
					}
				}
			}
			replace_row(qj,j,q);
		}
	
		for(int i=0;i<*N;i++){
			for(int k=0;k<*M;k++){
				for(int j=1;j<=*NJP;j++){
					//VE(minor[i],k) += ME(qj,j,k)*VE(dmc[i],j)/VE(yl,j);
					ME(minor,i,k) += ME(qj,j,k)*ME(dmc,i,j)/VE(yl,j);
				}
			}
		}
			
		for(int k=0;k<*M;k++){
			for(int i=0;i<*N;i++){   
				ME(eta,i,k) = ME(eta,i,k) + ME(minor,i,k);
			}
		}
	} //end if
	
	
	mat_zeros(Sigma);
	for(int i=0;i<*N;i++){
		extract_row(eta,i,tmpv1);
		vtv(tmpv1,tmpv1,tmpm1);
		mat_add(Sigma,tmpm1,Sigma);
	}
	MxA(SI,Sigma,tmpm1);
	MxA(tmpm1,SI,tmpm1);
	for(int k=0;k<*M;k++){	          //return to R	
		for(int l=0;l<*M;l++) beta_var[k*(*M)+l] = ME(tmpm1,l,k);  //by col
	}

	free_vecs(&beta,&U,&S1,&dS0,&xi,&E,&tmpv1,&delta,&sqSI,&rowX,&q,NULL);
	free_vecs(&Idxt,NULL);
	free_vecs(&survc,&dnc,&dlambc,&yl,NULL);
	free_mats(&I,&S2,&tmpm1,&tmpm2,&tmpm3,&Sigma,NULL);   //&SI
	free_mats(&qj,&minor,&dmc,NULL);
	// for(int i=0;i<*N;i++){
		// free_vecs(&minor[i],NULL);
		// free_vecs(&dmc[i],NULL);
	// }
}

void sort(::vector *in, ::vector *out){
	double tmp;
	vec_copy(in,out);
	int n = length_vector(in);
	for(int i=0;i<n-1;i++){
		for(int j=0;j<n-i-1;j++){
			if(VE(out,j)>VE(out,j+1)){  //">" means ties don't change order 
				tmp= VE(out,j);
				VE(out,j) = VE(out,j+1);
				VE(out,j+1) = tmp;
			}
		}
	}
	return;
}

void indsort(::vector *in, ::vector *out, ::vector *ind){	
	int n = length_vector(in);
	::vector *tmpv1;
	malloc_vec(n,tmpv1);
	vec_copy(in,tmpv1);
	sort(in,out);
	for(int j=0;j<n;j++){
		for(int i=0;i<n;i++){
			if(VE(out,j)==VE(in,i)){
				VE(ind,j)=i;
				VE(in,i)=-9999;  
				break;
			}
		}
	}	
	vec_copy(tmpv1,in);
	free_vec(tmpv1);	
	return;
}

int unisort(::vector *in, ::vector *out){
	int cnt=0;
	int m;
	::vector *tmpv1;
	int n = length_vector(in);
	malloc_vec(n,tmpv1);
	sort(in,tmpv1);
	for(int i=0;i<n-1;i++){
		if(VE(tmpv1,i)<VE(tmpv1,i+1)){
			VE(out,cnt)=VE(tmpv1,i);
			cnt += 1;
		}
		if(VE(tmpv1,i)>VE(tmpv1,i+1)) oops("Error in unisort()");
		if(i==n-2){
			VE(out,cnt)=VE(tmpv1,i+1);  //deal with last obs
			cnt += 1;
		}
	}
	m = cnt;
	free_vec(tmpv1);
	return m;
}

double vec_max(::vector *in){
	::vector *tmpv1;
	int n = length_vector(in);
	malloc_vec(n,tmpv1);
	sort(in,tmpv1);
	double max=VE(tmpv1,n-1);
	free_vec(tmpv1);
	return max;
}

double mat_max(matrix *in){
	double maximum = 0;
	int n=nrow_matrix(in);
	int m=ncol_matrix(in);
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			maximum = max(ME(in,i,j), maximum);
		}
	}
			
	return maximum;
}

// m := t(v)*t, JN
matrix *vtv(::vector *v1, ::vector *v2, matrix *m){
	int i,j;
	int n = length_vector(v1);
	if(!(length_vector(v2)==n)){
		oops("Error: dimensions in vtv\n");
	}
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			ME(m,i,j) = VE(v1,i)*VE(v2,j);
		}
	}
	return m;
}

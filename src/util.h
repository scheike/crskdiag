//util.h

#include <iostream>
#include <math.h>
#include <Rmath.h>
#include "matrix.h"
//#include <vector>

using namespace std;

void est_tmp(double *times, int *cause, matrix *WX, int *N, int *M, 
		int *Nit, int *NJP, double *TJP, double *betaS, double *beta_var, double *dlamb0, ::vector *Gbeta, 
		matrix *eta, matrix *wy, matrix *wdn1, matrix *wdm1, ::vector *S0j, matrix *Ej, matrix *Uj, matrix *SI, int *minor_included);
			
void sort(::vector *in, ::vector *out);

void indsort(::vector *in, ::vector *out, ::vector *ind);

int unisort(::vector *in, ::vector *out);

double vec_max(::vector *in);

double mat_max(matrix *in);

matrix *vtv(::vector *v1, ::vector *v2, matrix *m); 

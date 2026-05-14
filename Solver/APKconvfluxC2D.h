#ifndef __APKconvfluxC2D_h
#define __APKconvfluxC2D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKConvfluxIntC2D.h"
#include"APKConvfluxExtCorner2D.h"
#include"APKConvfluxExtCother2D.h"


int convfluxC2D(double *phicellold,double *phicellnew,double alpha,Zone *zone,
		double *Fc_n,double *Fn_d,double *convdia,int vartag,int scheme,int Index[3]);

#endif


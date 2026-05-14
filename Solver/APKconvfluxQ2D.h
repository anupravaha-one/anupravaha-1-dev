#ifndef __APKconvfluxQ2D_h
#define __APKconvfluxQ2D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKConvfluxIntQ2D.h"
#include"APKConvfluxExt2D.h"
#include"APKConvfluxExtQ2D.h"

int convfluxQ2D(double *phicellold,double *phicellnew,double alpha,Zone *zone,
		double *Fc_n,double *Fn_d,double *convdia,int vartag,int scheme,int IndexConv[3]);

#endif


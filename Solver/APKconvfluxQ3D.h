#ifndef __APKconvfluxQ3D_h
#define __APKconvfluxQ3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKConvfluxIntQ3D.h"
#include"APKConvfluxExt3D.h"
#include"APKConvfluxExtQ3D.h"

int convfluxQ3D(double *phicellold,double *phicellnew,double alpha,Zone *zone,double *Fc_n,
				double *Fn_d,double *convdia,int vartag,int scheme,int IndexConv[3]);

#endif


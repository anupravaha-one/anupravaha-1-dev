#ifndef __APKconvfluxC3D_h
#define __APKconvfluxC3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKConvfluxExtCorner3D.h"
#include"APKConvfluxIntC3D.h"
#include"APKConvfluxExtCother3D.h"

int convfluxC3D(double *,double *,double ,Zone *,double *,double *,double *,int ,int ,int IndexConv[3]);

#endif



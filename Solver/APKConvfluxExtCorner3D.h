#ifndef __APKConvfluxExtCorner3D_h
#define __APKConvfluxExtCorner3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double ConvfluxExtCorner3D(Zone *zone,double *,double *,double *,double *,double *,double *,
		int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int ,double );

#endif


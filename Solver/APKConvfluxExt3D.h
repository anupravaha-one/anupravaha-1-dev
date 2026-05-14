#ifndef __APKConvfluxExt3D_h
#define __APKConvfluxExt3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double ConvfluxExt3D(Zone *,double *,double *,double *,double *,double *,double *,
			int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int ,double);

double SOUIntpo(double volu,double voluu,double phiu,double phiuu);

#endif



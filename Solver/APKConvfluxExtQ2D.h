#ifndef __APKConvfluxExtQ2D_h
#define __APKConvfluxExtQ2D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKConvfluxIntQ3D.h"
 
double ConvfluxExtQ2D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,
		double *FUn,int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int face,double Alpha);


#endif


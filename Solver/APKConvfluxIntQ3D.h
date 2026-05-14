#ifndef __APKConvfluxIntQ3D_h
#define __APKConvfluxIntQ3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

int ConvfluxIntQ3D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,double *FUn,
		int index[3],double NCoeff[6],double Diagonal[6],double *slabvolume[6],double *slabPhi[6],double Alpha);

double QVIntpo(double volu,double voluu,double vold,double phiu,double phiuu,double phid);

#endif


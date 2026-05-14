#ifndef __APKDisperisonflux2d_h
#define __APKDisperisonflux2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double Disperisonflux2d(int index[3],double *domain,double *volumes,int *size,
	double *slabVariable[6],double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *SurfArea,double DiffFlux[6],double *GradX,double *GradY,double *slabGradX[6],double *slabGradY[6],double *DeltaInverse,double *slabDeltaInverse[6]);

double Disperisonflux3d(int index[3],double *domain,double *volumes,int *size,double *slabVariable[6],double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *SurfArea,double DiffFlux[6],double *GradX,double *GradY,double *GradZ,double *slabGradX[6],double *slabGradY[6],double *slabGradZ[6],double *DeltaI,double *slabDeltaInverse[6]);
#endif


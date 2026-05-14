#ifndef __APKDiffusion3dWithConstDiffCoeff_optiMPFV_h
#define __APKDiffusion3dWithConstDiffCoeff_optiMPFV_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"
double Diffusion3dWithConstDiffCoeff_optiMPFV(int index[3],Zone *zone,double *domain,double *volumes,int *size,
	double *slabVariable[6],double *slabVolume[6],double *beta,double DiffFlux[6]);

#endif

#ifndef __APKDiffusion2dWithVaryingDiffCoeff_h
#define __APKDiffusion2dWithVaryingDiffCoeff_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double Diffusion2dWithVaryingDiffCoeff(int index[3],double *domain,double *volumes,int *size,
		double *slabVariable[6],double *slabVolume[6],double *diffusionCoeffdomain,
		double *slabDiffusionCoeff[6],double *beta,double DiffFlux[6],double DiagonalN[6]);

#endif


#ifndef __APKDiffusion2d_h
#define __APKDiffusion2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double Diffusion2d(int index[3],double *domain,double *volumes,int *size,double *slabVariable[6],
		double *slabVolume[6],double *beta,double DiffFlux[6],double DiagonalN[6]);

#endif


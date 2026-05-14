#ifndef __APKdatumUnsat2d_h
#define __APKdatumUnsat2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"


int datumUnsat2d(Base *base,Zone *zone,double *volumes,double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *F);

int datumUnsat3d(Base *base,Zone *zone,double *volumes,double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *F);

#endif


#ifndef __APKinitializeSlabDispersionCoeff_h
#define __APKinitializeSlabDispersionCoeff_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

double *initializeSlabDispersionCoeff(int ,Zone *,double *,int ,int ,int *);
double *initializeslab_Delta(int i,int c,int dir,int *Point_Range);
#endif


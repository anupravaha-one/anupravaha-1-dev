#ifndef __APKCalcPressureGrad_Common_h
#define __APKCalcPressureGrad_Common_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

int CalcPressureGrad_Common(Base *base, double *,double *slabVariable[6],double *,double *slabVolume[6],double *,int index[3],int *,double grad[3]);

#endif


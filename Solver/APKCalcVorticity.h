#ifndef __APKCalcVorticity_h
#define __APKCalcVorticity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

#include"APKCalcPressureGrad.h"

int CalcVorticity(Base *base);
int CalcAbsVorticity(Base *base);

#endif


#ifndef __APKCalcVolfracFaceValue_h
#define __APKCalcVolfracFaceValue_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"

double CalcVolfracFaceValue(Base *,Zone *,double *);
double CalcVolfracFaceValue3d(Base *,Zone *,double *);

#endif


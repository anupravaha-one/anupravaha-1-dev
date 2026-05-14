#ifndef __APKUpdatePressure_h
#define __APKUpdatePressure_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMatrixMult.h"
#include"APKMatrixMultCorner.h"

double UpdatePressure(Base *base,Zone *zone);

#endif


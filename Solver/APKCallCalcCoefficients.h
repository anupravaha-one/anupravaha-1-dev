#ifndef __APKCallCalcCoefficients_h
#define __APKCallCalcCoefficients_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKCalcCoefficients.h"
#include"APKCalcCoefficientsCorner.h"
#include"APKCalcCoefficients2d.h"
#include"APKCalcCoefficientsCorner2d.h"

double CallCalcCoefficients(Base *base,Zone *zone,int g);

#endif


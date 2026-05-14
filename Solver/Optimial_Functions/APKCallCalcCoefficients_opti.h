#ifndef __APKCallCalcCoefficients_opti_h
#define __APKCallCalcCoefficients_opti_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKCalcCoefficients_opti.h"
#include"APKCalcCoefficientsCorner_opti.h"
#include"APKCalcCoefficients2d_opti.h"
#include"APKCalcCoefficientsCorner2d_opti.h"

double CallCalcCoefficients_opti(Base *base,Zone *zone,int g);

#endif


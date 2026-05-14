#ifndef __APKCallCalcCoefficientsForVariables_h
#define __APKCallCalcCoefficientsForVariables_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffFlag.h"
#include"APKCalcCoefficientsForVariables.h"
#include"APKCalcCoefficientsCornerForVariables.h"
#include"APKCalcCoefficients2dForVariables.h"
#include"APKCalcCoefficientsCorner2dForVariables.h"
#include"APKCalcCoefficientsForVariablesVarDiffCoeff.h"
#include"APKCalcCoefficientsCornerForVariablesVarDiffCoeff.h"
#include"APKCalcCoefficients2dForVariablesVarDiffCoeff.h"
#include"APKCalcCoefficientsCorner2dForVariablesVarDiffCoeff.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDomainArray.h"
double CallCalcCoefficientsForVariables(Base *base,int g,int vartag);

#endif


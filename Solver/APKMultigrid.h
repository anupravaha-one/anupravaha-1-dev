#ifndef __APKMultigrid_h
#define __APKMultigrid_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKExtractVariables.h"
#include"APKSmoothing.h"
#include"APKProlongation.h"
#include"APKResidFind.h"
#include"APKResidFind0.h"
#include"APKRestriction.h"
#include"APKCorrection.h"
#include"APKMultigridExtractVariables.h"
#include"APKUpdatePressure.h"
#include"APKUpdatePressureVarDensity.h"
#include"APKFillVariable2Dummy.h"
#include"APKDummyVariableProlongate.h"
#include"APKDummyVariableRestrict.h"
#include"APKInterZoneVariableTransfer.h"

double Multigrid(Base *base,int g);

double Initialize_CoarseGrid(double *P,int *Size,double constant,int g);
double Initialize_CoarseGridSlabs(int f,double *P,int *Size,double constant,int gdi);

#endif


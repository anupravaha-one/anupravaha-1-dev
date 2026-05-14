#ifndef __APKMultigrid2d_h
#define __APKMultigrid2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKExtractVariables.h"
#include"APKSmoothing2d.h"
#include"APKProlongation2d.h"
#include"APKResidFind2d.h"
#include"APKResidFind2d0.h"
#include"APKRestriction.h"
#include"APKCorrection.h"
#include"APKMultigridExtractVariables.h"
#include"APKUpdatePressure2d.h"
#include"APKUpdatePressure2dVarDensity.h"
#include"APKFillVariable2Dummy2d.h"
#include"APKDummyVariableProlongate2d.h"
#include"APKDummyVariableRestrict2d.h"
#include"APKInterZoneVariableTransfer2d.h"

double Multigrid2d(Base *base,int g);

double Initialize_CoarseGrid2d(double *P,int *Size,double constant,int g);
double Initialize_CoarseGridSlabs2d(int f,double *P,int *Size,double constant,int gdi);

#endif


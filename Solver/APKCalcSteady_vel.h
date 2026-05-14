#ifndef __APKCalcSteady_vel_h
#define __APKCalcSteady_vel_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKTimeDepend.h"

#include"APKData_Structures.h"
#include"APKCalculateTimeStep.h"
#include"APKMakeDomainArray.h"
#include"APKprintsolver.h"

#include"APKWriteCgnsFile.h"

#include"APKFillNewsolToPrevsol.h"
#include"APKRefreshSurfaceCells.h"
#include"APKSolutions.h"
#include"APKSolutions2d.h"
#include"APKCalcVariable.h"
#include"APKCalcVariableConjugate.h"
#include"APKMassFluxConv.h"
#include"APKVarMassFluxConv.h"
#include"APKShearStress.h"
#include"APKSaveTimeStepData.h"
#include"APKEvaluateVarCoefficients.h"
#include"APKEvaluateEnergysource.h"
#include"APKCDcoeffMemAlloc.h"
#include"APKDiffusionCoeffCalculation.h"
#include"APKDiffusionCoeffCalculationForVariables.h"
#include"APKCDcoeffCalculation.h"

int CalcSteady_vel(Root CgnsRoot);
double CalculateError(Base *base);
#endif


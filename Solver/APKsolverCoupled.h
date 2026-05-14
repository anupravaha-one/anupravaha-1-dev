#ifndef __APKsolverCoupled_h
#define __APKsolverCoupled_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKTimeDepend.h"

#include "APKprintsolver.h"
#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"

#include "APKWriteCgnsFile.h"

#include "APKCalcNusseltNo.h"
#include "APKRefreshSurfaceCells.h"
#include "APKShearStress.h"
#include "APKsavevolumeflux.h"
#include "APKvolumefluxstar.h"
#include "APKPressureIter.h"
#include "APKCalcNewVolumeFlux.h"
#include "APKFindVolumeFluxDifference.h"
#include "APKCalcVariableOneIter.h"
#include "APKCalcVariableOneIterConjugate.h"
#include "APKFillNewsolToPrevsol.h"
#include "APKmassVelocityOneIter.h"
#include "APKMassFluxConv.h"
#include "APKVarMassFluxConv.h"
#include "APKSaveTimeStepData.h"
#include "APKEvaluateVarCoefficients.h"
#include "APKRefreshVariableCoefficients.h"
#include "APKEvaluateEnergysource.h"
#include"APKCDcoeffMemAlloc.h"
#include"APKDiffusionCoeffCalculation.h"
#include"APKDiffusionCoeffCalculationForVariables.h"
#include"APKCDcoeffCalculation.h"
//Size of the arrays are increased to 22 by adhi+shash to incoroporate new variables
double CalcVariableError(int variableTag,Base *base);
int FindTimeStepErrorArray(Base *base,double rmsError[27]);
int FindInsideTimeStepConverToken(double ErrorNew[14],double ErrorBeforeIteration[14]);
int FindtimeStepToken(Base *base,double ErrorTimeStep[27],int CounterTstep);
int FindVelocityIterToken(Base *base,double velocityResidual[3]);
int FindVariableIterToken(double variableResidual[11],Base *base);
int printTimeStepResid(Base *base,double ErrorTimeStep[27]);
int printTimeStepResidOneline(Base *base);

int solverCoupled(Root CgnsRoot);

#endif


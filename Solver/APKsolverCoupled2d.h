#ifndef __APKsolverCoupled2d_h
#define __APKsolverCoupled2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKTimeDepend.h"

#include"APKData_Structures.h"
#include"APKprintsolver.h"

#include"APKWriteCgnsFile.h"

#include"APKCalcNusseltNo.h"
#include"APKCalcVariableOneIterConjugate.h"
#include"APKRefreshSurfaceCells.h"
#include"APKShearStress.h"
#include"APKsolverCoupled.h"
#include"APKPressureIter.h"
#include"APKFindVolumeFluxDifference.h"
#include"APKCalcVariableOneIter.h"
#include"APKFillNewsolToPrevsol.h"
#include"APKmassVelocityOneIter.h"
#include"APKMassFluxConv.h"
#include"APKVarMassFluxConv.h"
#include"APKsavevolumeflux.h"
#include"APKCalcNewVolumeFlux.h"
#include"APKvolumefluxstar.h"
#include"APKSaveTimeStepData.h"
#include"APKEvaluateVarCoefficients.h"
#include"APKRefreshVariableCoefficients.h"
#include"APKEvaluateEnergysource.h"
#include"APKVarMassFluxConv.h"
#include"APKCDcoeffMemAlloc.h"
#include"APKDiffusionCoeffCalculation.h"
#include"APKDiffusionCoeffCalculationForVariables.h"
#include"APKCDcoeffCalculation.h"

int solverCoupled2d(Root CgnsRoot);

#endif



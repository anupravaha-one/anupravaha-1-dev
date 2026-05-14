#ifndef __APKsolverSegregated_h
#define __APKsolverSegregated_h

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
#include"APKCalcConjugate.h"
#include "APKvolumeflux.h"
#include "APKprint_monitors.h"
#include "APKWriteCgnsFile_vertex.h"
#include "APKRadiationSolverWithFlow2D.h"
#include "APKRadiationSolverWithFlow3D.h"

#include "APKCalcSource_opti.h"
#include "APKCallCalcCoefficientsForVariables.h"
#include "APKConv_Coeff_Calculation.h"
#include "APKRefreshVariableCoefficients.h"

#include"APK_File_Open.h"
#include"APKWrite_Solution.h"

#include"APKModiCalcConjugate.h"
#include"APKConj_CalcConjugate.h"

int solverSegregated(Root CgnsRoot);
int FindTimeStepTokenSegr(Base *base,double errorTimeStep[11]);
double CalculateError(Base *base);
double CalculateErrorTemp(Base *base,int vartag);
int make_w_zero(Base* base);
int find_velocities(Base *);

#endif


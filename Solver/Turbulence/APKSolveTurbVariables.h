#ifndef __APKSolveTurbVariables_h
#define __APKSolveTurbVariables_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKCalcNuDiffCoef.h"
#include"APKCalcSource_opti.h"
#include"APKCalcVariable.h"
#include"APKCalcVorticity.h"
#include"APKCallCalcCoefficientsForVariables.h"
#include"APKCDcoeffCalculation.h"
#include"APKConv_Coeff_Calculation.h"
#include"APKData_Structures.h"
#include"APKFillNewsolToPrevsol.h"
#include"APKinitializeTurbVariables.h"
#include"APKLimitVariable.h"
#include"APKRefreshSurfaceCells.h"
#include"APKTimeDepend.h"
#include"APKTURBCalcGamma.h"
#include"APKTURBCalcTurbViscosity.h"
#include"APKTURBCalcUtau.h"
#include"APKTURBCalcC1epsi.h"
#include"APKTURBCalcEpsilonWallTerm.h"
#include"APKTURBCalcKWallTerm.h"
#include"APKTURBCalcMuEff.h"
#include"APKTURBCalcTurbulentProduction.h"
#include"APKTURBCalcVariableOneIterK.h"
#include"APKTURBCalcVariableOneIterEps.h"
#include"APKTURBCalcVariableOneIterOmega.h"
#include"APKTURBCalcVariableOneIterOmegaSST.h"
#include"APKTURBCrossDiffusion.h"
#include"APKTURBDistance2Wall.h"
#include"APKTURBfreeCrossDiff.h"
#include"APKTURBfreeGamma.h"
#include"APKTURBCalcTurbconductivity.h"

int SolveTurbVariables(Root CgnsRoot,double variableResidual[11], int CounterTstep);

#endif


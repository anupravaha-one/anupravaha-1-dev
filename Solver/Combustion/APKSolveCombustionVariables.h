#ifndef __APKSolveCombustionVariables_h
#define __APKSolveCombustionVariables_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>



#include "APKData_Structures.h"
#include "APKCalcVariableOneIter.h"
#include "APKsolverCoupled.h"
#include "APKTURBCalcMuEff.h"
#include "APKCDcoeffCalculation.h"
#include "APKCallCalcCoefficientsForVariables.h"
#include "APKCalcSource_opti.h"
#include "APKCalcVariable.h"

int SolveCombustionVariables(Base *);
int FindGSResidualCombVariables(double *);


#endif


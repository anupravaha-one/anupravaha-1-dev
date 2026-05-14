#ifndef __APKSolverElectricSemiCoupled2d_h
#define __APKSolverElectricSemiCoupled2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"

#include "APKFillNewsolToPrevsol.h" 
#include "APKWriteCgnsFile.h"

#include "APKCalcVariableCharge.h"
#include "APKCalcVariablePotential.h"

#include "APKprintsolver.h"
#include "APKWireCurrentDensity.h"
#include "APKElectricFieldXYZ.h"
#include "APKElectricFluxXYZ.h"
#include "APKSaveTimeStepData.h"
#include "APKCDcoeffCalculation.h"

#include "APKCoeffcalculation_ESP.h"
#include "APKFree_Coefficients_ESP.h"
#include "APKCalcSource_opti.h"
#include "APKESP_Source_Calculation.h"
#include "APKConv_Coeff_Calculation.h"
#include "APKCallCalcCoefficientsForVariables.h"

#include"APK_File_Open.h"

int SolverElectricSemiCoupled2d(Root CgnsRoot);
double CalcErrorTimeStep(Base *base,int vartag, double delta_t);


#endif


#ifndef __APKSolverElectricSemiCoupled_h
#define __APKSolverElectricSemiCoupled_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"

#include "APKFillNewsolToPrevsol.h" 
#include "APKWriteCgnsFile.h"

#include "APKCalcVariablePotential.h"
#include "APKCalcVariableCharge.h"


#include "APKprintsolver.h"
#include "APKElectricFieldXYZ.h"
#include "APKElectricFluxXYZ.h"
#include "APKWireCurrentDensity.h"
#include "APKSaveTimeStepData.h"
#include "APKSolverElectricSemiCoupled2d.h"


#include "APKCoeffcalculation_ESP.h"

#include "APKCalcSource_opti.h"
#include "APKConv_Coeff_Calculation.h"
#include "APKCallCalcCoefficientsForVariables.h"
#include "APKCDcoeffCalculation.h"

#include "APKFree_Coefficients_ESP.h"
#include "APKESP_Source_Calculation.h"

#include"APK_File_Open.h"

int SolverElectricSemiCoupled(Root CgnsRoot);

#endif


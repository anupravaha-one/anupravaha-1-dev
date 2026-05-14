#ifndef __APKsolverLMMHD_h
#define __APKsolverLMMHD_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKCalcUcrossB.h"
//#include"APKCalcUcrossB2d.h"
#include"APKPotentialMHDpoissonRHS.h"
#include"APKPotentialMHDIter.h"
#include"APKCalcCurrentDensity.h"
#include"APKCalcVariableOneIter.h"
//#include"APKPsiIter.h"
#include"APKCalcNewMagneticflux.h"
//#include"APKCalctotalMagneticfield.h"
#include"APKUpdateinducedB.h"
#include"APKUpdateinducedB2d.h"
#include"APKMagneticflux.h"
//#include"APKMagneticflux2d.h"
#include"APKCalcConjugate.h"
#include"APKsavevolumeflux.h"
#include"APKCDcoeffCalculation.h"
#include"APKCalcVariable.h"
#include "APKCalcSource_opti.h"
#include "APKInduced_Source_Calculation.h"
#include "APKCalcPsiSource.h"
#include "APKsavemagneticflux.h"
#include "APKCallCalcCoefficientsForVariables.h"
#include "APKConv_Coeff_Calculation.h"

int solverLMMHD(Base *base);
//int FindBstarIterToken(double BstarResidual[3]);

#endif


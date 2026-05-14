#ifndef __APKSolutions2d_h
#define __APKSolutions2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKmassVelocity.h"
#include"APKPressureIter.h"
#include"APKPressureIterArtificial.h"
#include"APKCalcVariable.h"
#include"APKCalculateErrorPressure.h"
#include"APKfillnewtoprev.h"
#include"APKCalcNewVolumeFlux.h"
#include"APKFindVolumeFluxDifference.h"
#include"APKsavevolumeflux.h"
#include"APKvolumefluxstar.h"
#include"APKCDcoeffCalculation.h"

#include"APKPressure_Grad_Calculation.h"
#include"APKCalcSource_opti.h"
#include"APKConv_Coeff_Calculation.h"

#include"APKConj_PressureIter.h"
#include"APKConj_CalcConjugate.h"
#include"APKModiPressureIter.h"
#include"APKModiCalcConjugate.h"

int Solutions2d(Base *base,int counter_I);

#endif


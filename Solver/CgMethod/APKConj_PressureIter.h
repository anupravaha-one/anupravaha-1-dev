#ifndef __APKConj_PressureIter_h
#define __APKConj_PressureIter_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKRefreshSurfaceCells.h"
#include"APKExtractVariables.h"

#include "APKCalcInitPressureResidue2d.h"
#include "APKCalcLHS.h"
#include "APKUpdatePressure2d.h"

#include"APKUpdatePressure.h"
#include "APKCalcInitPressureResidue.h"
#include "APKCalcLHS3d.h"

#include "APKCalcInitPressureResidue2dVarDensityMassFlux.h"
#include"APKUpdatePressure2dVarDensityMassFlux.h"
#include "APKCalcLHSVarDensity.h"

#include"APKUpdatePressureVarDensity.h"
#include"APKUpdatePressure2dVarDensity.h"
#include"APKUpdatePressureVarDensityMassFlux.h"
#include"APKCalcExtraTermsOfVarDensity2d.h"
#include"APKCalcExtraTermsOfVarDensity.h"
#include"APKCalcExtraTermsOfVarDensityMassflux2d.h"
#include"APKCalcExtraTermsOfVarDensityMassflux.h"


#include"APK_File_Open.h"

double Conj_PressureIter(Base *base);

#endif


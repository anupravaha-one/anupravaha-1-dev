#ifndef __APKinitialize_h
#define __APKinitialize_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKFindZonePtr.h"
#include"APKinitialize_variable.h"
#include"APKinitializeslab.h"
#include"APKinitializeVarCoefficient.h"
#include"APKinitializeSlabVarCoeff.h"
#include"APKinitialize_Effviscosity.h"
#include"APKinitialize_Effviscosityofphase.h"
#include "APKinitialize_Dummyeffviscosityofphase.h"
#include"APKinitialize_T.h"
#include"APKinitialize_enthalpy.h"
//#include"APKinitialize_solidification.h"
#include"APKinitialize_TwoPhase.h"
#include"APKinitializeCombParameters.h"
#include"APKinitializeCombVariables.h"
#include"APKinitializeElectric.h"
#include"APKinitializeLMMHD.h"
#include"APKinitializeMPFV.h"
#include"APKinitializeDisperse.h"
#include"APKData_Structures.h"
#include"APKinitialize_monitors.h"


int initialize(Base *);

#endif


#ifndef __APKinitializeGroundwater_h
#define __APKinitializeGroundwater_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKinitialize_variable.h"
#include"APKinitializeslab.h"


#include"APKinitializeVarCoefficient.h"
#include"APKinitializeSlabVarCoeff.h"
#include"APKinitialize_Effviscosity.h"

#include"APKinitializeDispersionCoeff.h"
#include"APKinitializeSlabDispersionCoeff.h"

#include"APKData_Structures.h"


#include"APKFindInterfaceBCPtrMod.h"
#include"APKupdateBCDeltaInverse.h"
#include"APKFindZonePtr.h"

int initializeGroundwater(Base *);

#endif


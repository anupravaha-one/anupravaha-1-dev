#ifndef __APKTURBCalcTurbViscosity_h
#define __APKTURBCalcTurbViscosity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFillerInterfaceVariables.h"
#include "APKFindZonePtr.h"
#include "APKHomogeneousDirichletIC.h"
#include "APKHomogeneousNeumannIC.h"
#include "APKNonOrthogonalHomogeneousNeumannIC.h"
#include "APKTURBminmax.h"

int CalcTurbViscosity(Base *base);

#endif


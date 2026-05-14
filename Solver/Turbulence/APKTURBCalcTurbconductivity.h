#ifndef __APKTURBCalcTurbconductivity_h
#define __APKTURBCalcTurbconductivity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKHomogeneousDirichletIC.h"
#include "APKHomogeneousNeumannIC.h"
#include "APKNonOrthogonalHomogeneousNeumannIC.h"

int CalcTurbconductivity(Base *base);

#endif

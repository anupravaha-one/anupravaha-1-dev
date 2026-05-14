#ifndef __APKTURBCalcTurbulentProduction_h
#define __APKTURBCalcTurbulentProduction_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

#include"APKTURBgrad3D.h"
#include"APKTURBgrad2D.h"

int CalcTurbulentProduction(Base *base);
int TurbProduction(Base *base,Zone *zone,double *Gk);

#endif


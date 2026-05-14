#ifndef __APKTURBCalcEpsilonWallTerm_h
#define __APKTURBCalcEpsilonWallTerm_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"../Common/APKData_Structures.h"
#include"../Common/APKMakeSlabPtrArray.h"
#include"../Common/APKMakeDomainArray.h"

#include"APKTURBCalcDoubleDerivative2d.h"
#include"APKTURBCalcDoubleDerivative3d.h"

int CalcEpsWallTerm(Base *base);
int EpsilonWall(Base *base,Zone *zone,double *Gk);

#endif

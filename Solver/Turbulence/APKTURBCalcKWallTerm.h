#ifndef __APKTURBCalcKWallTerm_h
#define __APKTURBCalcKWallTerm_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"../Common/APKData_Structures.h"
#include"../Common/APKMakeSlabPtrArray.h"
#include"../Common/APKMakeDomainArray.h"

#include"APKTURBKgradient2D.h"
#include"APKTURBKgradient3D.h"

int CalcKWallTerm(Base *base);
int KWall(Base *base,Zone *zone,double *Gk);

#endif

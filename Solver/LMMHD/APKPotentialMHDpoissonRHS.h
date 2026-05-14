#ifndef __APKPotentialMHDpoissonRHS_h
#define __APKPotentialMHDpoissonRHS_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"
#include"APKCalcUcrossB.h"

int PotentialMHDpoissonRHS(Base *base);
int PotentialMHDpoissonRHS2d(Base *base);

#endif


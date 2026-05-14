#ifndef __APKInterpolateEfield_h
#define __APKInterpolateEfield_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"
#include"APKFindCellCentre.h"

int InterpolateEfield (Base *base,double XP,double YP,double ZP);

#endif 

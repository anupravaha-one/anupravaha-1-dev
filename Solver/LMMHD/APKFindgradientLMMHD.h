#ifndef __APKFindgradientLMMHD_h
#define __APKFindgradientLMMHD_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKCalcGradVar.h"


double FindgradientLMMHD(Base *base,Zone *zone,int index[3],double gradx[3],double grady[3],double gradz[3]);

#endif


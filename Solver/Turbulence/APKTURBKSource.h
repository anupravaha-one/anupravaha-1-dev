#ifndef __APKTURBKSource_h
#define __APKTURBKSource_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKCalcGradVar.h"
#include "APKMakeDomainArray.h"


double KSource(Base *base,Zone *zone,int index[3],int variableTag);

#endif
#ifndef __APKSpalartAllmarasSource_h
#define __APKSpalartAllmarasSource_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKCalcGradVar.h"
#include"APKMakeSlabPtrArray.h"

double SpalartAllmarasSource(Base *base,Zone *zone,double *source,int index[3]);

#endif


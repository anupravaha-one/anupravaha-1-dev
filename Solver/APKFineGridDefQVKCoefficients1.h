#ifndef __APKFineGridDefQVKCoefficients1_h
#define __APKFineGridDefQVKCoefficients1_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKCalcvarCoeffAndPrevCoeff.h"
#include"APKFineGridCoefficients1.h"

double FineGridDefQVKCoefficients1(Base *base,Zone *zone,int vartag);

#endif


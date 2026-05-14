#ifndef __APKCalcSource_h
#define __APKCalcSource_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include "APKCalcPressureGrad.h"
#include "APKFindgradient.h"
#include "APKFindgradientPhiMHD.h"
#include "APKSpalartAllmarasSource.h"
#include "APKTURBCalcGradSquare.h"

double CalcSource(double *source,int variableTag,Zone *zone,Base *base,int index[3]);

#endif

#ifndef __APKCalcSource_opti_h
#define __APKCalcSource_opti_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>


#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKFindgradient.h"
#include "APKFindgradient_Disperse.h"
#include "APKSpalartAllmarasSource.h"
#include "APKTURBCalcGradSquare.h"
#include "APKTURBminmax.h"
#include "APKFind_Mushy_Source.h"
#include "APKTURBKSource.h"

double CalcSource_opti(Base *base,int variableTag);

#endif


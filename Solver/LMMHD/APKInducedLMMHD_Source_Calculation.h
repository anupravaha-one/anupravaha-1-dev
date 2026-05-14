#ifndef __APKInducedLMMHD_Source_Calculation_h
#define __APKInducedLMMHD_Source_Calculation_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKCalcPressureGrad_Common.h"
//#include"APKCalcLMMHDGrad_Common.h"

#include"APKCalcKhiGrad_Common.h"
#include"APKInducedMomentumSource2d.h"
#include"APKInducedMomentumSource.h"
#include"APKFindgradientLMMHD.h"
int InducedLMMHD_Source_Calculation(Base * base);

#endif


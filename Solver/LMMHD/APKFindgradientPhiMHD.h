#ifndef __APKFindgradientPhiMHD_h
#define __APKFindgradientPhiMHD_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKCalcGradVar.h"
#include"APKCalcPressureGrad.h"
#include"APKMakeDomainArray.h"


double FindgradientPhiMHD(Base *base,Zone *zone,int index[3],double Dphi[3]);

#endif


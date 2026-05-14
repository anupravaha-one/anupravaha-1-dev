#ifndef __APKPressureHeadIter_h
#define __APKPressureHeadIter_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"


#include"APKUpdatePressureHead2d.h"
#include"APKUpdatePressureHead3D.h"
#include"APKRefreshSurfaceCells.h"
#include"APKExtractVariables.h"
#include"APKUpdateWaterContent.h"
#include"APKRefreshVariableCoefficients.h"
		
int PressureHeadIter(Base *base);
int PressureHeadIter3D(Base *base);

#endif


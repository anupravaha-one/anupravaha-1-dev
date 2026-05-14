#ifndef __APKPotentialMHDIter_h
#define __APKPotentialMHDIter_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKUpdatePotentialMHD.h"
#include"APKUpdatePotentialMHD2d.h"
#include"APKRefreshSurfaceCells.h"
#include"APKExtractVariables.h"
#include "APKUpdateVarCD2d.h"
#include "APKUpdateVarCD.h"


double PotentialMHDIter(Base *base);

#endif


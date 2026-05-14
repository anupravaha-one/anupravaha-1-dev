#ifndef __APKModiCalcConjugate_h
#define __APKModiCalcConjugate_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKRefreshSurfaceCells.h"
#include"APKExtractVariables.h"

#include "APKCalcCorrection.h"
#include "APKCalcLHS.h"

#include "APKCalcInitConjVarResidue2d.h"

#include "APKEvaluateEnergysource.h"
#include "APKCalcLHS3d.h"
#include "APKCalcCorrection3d.h"
int ModiCalcConjugate(Base *base,int );

#endif


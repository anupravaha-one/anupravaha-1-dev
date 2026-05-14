#ifndef __APKConj_CalcConjugate_h
#define __APKConj_CalcConjugate_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKRefreshSurfaceCells.h"
#include"APKExtractVariables.h"

#include "APKCalcInitConjVarResidue2d.h"
#include "CgMethod/APKCalcLHS.h"

#include "APKEvaluateEnergysource.h"
#include"APKinitialize_variable.h"

double Conj_CalcConjugate(Base *base,int );

#endif


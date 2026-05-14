#ifndef __APKSolverElectric_h
#define __APKSolverElectric_h

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>

#include "../Common/APKData_Structures.h"


//#include "APKSolverElectricCoupled.h"
//#include "APKSolverElectricCoupled2d.h"
#include "APKSolverElectricSemiCoupled.h"
#include "APKSolverElectricSemiCoupled2d.h"
//#include "APKSolverElectricSegregated.h"
#include "APKElectricFieldXYZ.h"

#include "APKWriteCgnsFile.h"


int SolverElectric(Root CgnsRoot);

#endif

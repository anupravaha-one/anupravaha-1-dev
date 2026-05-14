#ifndef __APKCheckConservation_h
#define __APKCheckConservation_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKTimeDepend.h"

#include"APKData_Structures.h"
#include"APKprintsolver.h"
#include"APKprintslabs.h"

#include"APKWriteCgnsFile.h"
#include"APKFindCellCentre.h"
#include"APKRadiationSolverWithFlow2D.h"

int CheckConservation(Root CgnsRoot,int tstep);

#endif



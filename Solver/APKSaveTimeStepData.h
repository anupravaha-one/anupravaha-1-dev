#ifndef __APKSaveTimeStepData_h
#define __APKSaveTimeStepData_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKprintslabs.h"
#include"APKprintsolver.h"

#include"APKWriteCgnsFile.h"
#include"APKWriteCgnsFile_vertex.h"

#include"APKCalcVorticity.h"
#include "APKCalcCurrentDensity.h"
#include "APKElectricFieldXYZ.h"

#include"APK_File_Open.h"
#include"APKWrite_Solution.h"

int SaveTimeStepData(Root CgnsRoot,int cntTimeStep,double time);
char *FindFileNamePtr(Root CgnsRoot,double time);

#endif


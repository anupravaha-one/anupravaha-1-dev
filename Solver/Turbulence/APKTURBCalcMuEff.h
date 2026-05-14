#ifndef __APKTURBCalcMuEff_h
#define __APKTURBCalcMuEff_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFillerInterfaceVariables.h"
#include "APKFindZonePtr.h"
#include "APKHomogeneousNeumannIC.h"

#include"APK_File_Open.h"

int CalcMuEff(Base *base, int vartag);

#endif


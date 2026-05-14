#ifndef __APKExtract_Var_Coeff_PorousMedia_h
#define __APKExtract_Var_Coeff_PorousMedia_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKFillerInterfaceVariables.h"
#include "APKUpdateBCslabsKunsat.h"  
#include "APKupdateBCslabsVar_Coeff_PorousMedia.h"
#include "APKFindZonePtr.h"

int Extract_Var_Coeff_PorousMedia(int variableTag,Base *base,Zone *zone);
#endif

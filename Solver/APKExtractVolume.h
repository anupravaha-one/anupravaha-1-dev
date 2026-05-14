#ifndef __APKExtractVolume_h
#define __APKExtractVolume_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKFillerConjugateBoundary.h"
#include "APKFillerInterface.h"
#include "APKFillerBoundary.h"

#include "APKFindZonePtr.h"

int ExtractVolume(Base *base,Zone *zone);

#endif


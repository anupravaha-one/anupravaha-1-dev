#ifndef __APKCalcNewVolumefluxRestart_h
#define __APKCalcNewVolumefluxRestart_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"

#include "APKDiffusion.h"
#include "APKDiffusion2d.h"

int CalcNewVolumefluxRestart(Base *base);

#endif


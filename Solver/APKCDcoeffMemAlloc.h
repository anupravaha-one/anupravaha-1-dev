#ifndef __APKCDcoeffMemAlloc_h
#define __APKCDcoeffMemAlloc_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKCountCells.h"

#include "APKFineGridCoefficients1.h"
#include "APKFineGridDefCDCoefficients1.h"
#include "APKFineGridDefQVKCoefficients1.h"
#include "APKFineGridCoefficients2.h"
#include "APKFineGridDefCDCoefficients2.h"
#include "APKFineGridDefQVKCoefficients2.h"
#include "APKCoarseGridCoefficients.h"

int CDcoeffMemAlloc(Base *base);

#endif


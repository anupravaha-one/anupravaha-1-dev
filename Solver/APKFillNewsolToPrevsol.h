#ifndef __APKFillNewsolToPrevsol_h
#define __APKFillNewsolToPrevsol_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKCNterm.h"
#include"APKUpdateBoundarySlabs.h"
#include"APKUpdateOldZSlabs.h"
#include"APKTURBUpdateOldZSlabsKEpsilon.h"

int FillNewsolToPrevsol(Base *base);
int FillNewsolToPrevsolKepsilon(Base *base);
int FillNewsolToPrevsol_T(Base *base,int vartag);

#endif

#ifndef __APKboundary_d_t_sTurb_h
#define __APKboundary_d_t_sTurb_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

#include "APKMakeSlabPtrArray.h"
#include "APKDiffusionCommon.h"
#include "APKCalcDiffusionCoefCommon.h"
#include "APKFindSlabPtrBC.h"
#include "APKfind_convcri.h"
#include "APKupdate_boundaryTurbDMPFV.h"
#include "APKWallFunctionDirichlet.h"

double boundary_d_t_sTurb(Base* base,ZoneBoundary *zBc,int direction,double *flowsol,int vartag,Zone *z,int size,Dataset *dset,double *dummy);

#endif


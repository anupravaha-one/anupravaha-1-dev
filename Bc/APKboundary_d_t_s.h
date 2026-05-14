#ifndef __APKboundary_d_t_s_h
#define __APKboundary_d_t_s_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"
#include "../Common/APKMakeSlabPtrArray.h"
#include "../Common/APKDiffusionCommon.h"
#include "../Common/APKCalcDiffusionCoefCommon.h"
#include "../Common/APKFindSlabPtrBC.h"
#include "APKfind_convcri.h"
#include "APKboundary_d_t_sRadiation.h"
#include "APKboundary_d_t_sRadsymmetry.h"
#include "APKboundary_d_t_sRadinsulation.h"
#include "APKboundary_d_t_sRobin.h"
#include "APKboundary_d_t_sInsulation.h"


double calcR(int i,int j,int ni,double *Ycor,double *R);
double boundary_d_t_s(Base* base ,int direction,double *flowsol,int vartag,Zone *z,int size,Dataset *dset,double *dummy);

#endif


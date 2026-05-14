#ifndef __APKCalculateGeometry_h
#define __APKCalculateGeometry_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
//#include "APKFindInterfaceBCPtr.h" //Not needed

#include "APKinitializeslab.h"

#include "APKFindCellCentre.h"
#include "APKGetSurfacePoints.h"
#include "APKFindZonePtr.h"
#include "APKFindFaceCenterBC.h"
#include "APKTransformCellIndices1.h"
#include "APKAdjoiningZonePointerXYZ.h"
#include "APKCalSurf.h"
#include "APKCheckOutward.h"
#include "APKFindBeta.h"
#include "APKCalculateDmF.h"
#include "APKCalculateNormalVector.h"

int CalculateGeometry(Base *CgnsBase);
int FindDeltaX(double center[3],double centerN[3],double delta_X[3]);
int FindDeltaInverse2D(double delta[4],double *addressDeltaInverse);
int FindDeltaInverse3D(double delta[9],double *addressDeltaInverse);
int ChangeFaceCenter(int face, double center[3],double centerN[3],double delta_X[3],double surfVector[3]);
int ModifyBeta(int face,double beta[3],char tag,int Tag);

#endif


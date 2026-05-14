#ifndef __APKUpdateVarQVK_h
#define __APKUpdateVarQVK_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
#include"APKCalcvarCoeffAndPrevCoeff.h"
#include"APKCalcSource.h"
#include"APKconvfluxQ2D.h"
#include"APKconvfluxQ3D.h"
#include"APKconvfluxC3D.h"
#include"APKconvfluxC2D.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffFlag.h"
#include"APKConvfluxIntC3D.h"
#include"APKConvfluxIntQ3D.h"
#include"APKConvfluxExtQ3D.h"
#include"APKConvfluxExt3D.h"
#include"APKConvfluxExtCorner3D.h"
#include"APKConvfluxExtCother3D.h"
#include"APKConvfluxIntC2D.h"
#include"APKConvfluxIntQ2D.h"
#include"APKConvfluxExtQ2D.h"
#include"APKConvfluxExt2D.h"
#include"APKConvfluxExtCorner2D.h"
#include"APKConvfluxExtCother2D.h"
#include"APKMatrixMult.h"


double UpdateVarQVK(Base *base,Zone *zone,int vartag);

#endif


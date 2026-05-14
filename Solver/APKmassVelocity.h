#ifndef __APKmassVelocity_h
#define __APKmassVelocity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKExtractMassVelocities.h"
#include"APKExtractVariables.h"
#include"APKfillnewtoprev.h"
#include"APKfind_convcri.h"
#include"APKfindresidual.h"
#include"APKfindresidualQVK.h"
#include"APKfindresidualQVK2d.h"
#include"APKfindresidualCD.h"
#include"APKfindresidualCD2d.h"
#include"APKMakeDomainArray.h"
#include"APKRefreshMassVelocities.h"
#include"APKRefreshSurfaceCells.h"

int massVelocity(Base *base,int variableTag);

#endif


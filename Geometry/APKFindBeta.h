#ifndef __APKFindBeta_h
#define __APKFindBeta_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

int FindBeta(int face,double centreCell[3],double centreNeighbour[3],double facePoints[][3],double surf[3],double beta[3]);

#endif


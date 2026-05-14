#ifndef __APKFindCellTransformRule_h
#define __APKFindCellTransformRule_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

int GetDonorRange(int transform[3],int pDR[6]);
int GetPointRange(int transform[3],int pDR[6]);
int FindMinIndices(int vertices[][3],int *cell);
int FindInterfaceCellVertices1(int point_s[3],int Vertices[8][3],int dir);
int FindInterfaceCellVertices2(int v1[8][3],int *trans,int v2[8][3],int *begin1,int *begin2);
//int FindInterfaceCellVertices2(int v1[8][3],int trans[3][3],int v2[8][3],int *begin1,int *begin2);
int FindAddConstant(int cell1[3],int cell2[3],int *trans,int addconst[3]);
//int FindAddConstant(int cell1[3],int cell2[3],int trans[3][3],int addconst[3]);
int FindCellTransformRule(Interface *face);

#endif


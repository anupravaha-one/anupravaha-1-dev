#include "APKMFillerInterfaceVariables2d.h"

int MFillerInterfaceVariables2d(int variableTag,Interface *facePtr,Zone *zone,Zone *ptrZone,int iface)
{
	int i,j,k,indices1[3],indices2[3],targetIndices;
	int *ptrRange;
	int totalCells;
	
	ptrRange=facePtr->Point_Range;
	
	switch(iface)
	{
		case 0:/* Imax*/
			
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+3),indices1[1]=j,indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				MfindZone2Variable2d(variableTag,facePtr,((zone->grid)->DummySlabV[1]+targetIndices),ptrZone,indices2,targetIndices);
			}
		break;
		case 1:/*Imin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+0)-1,indices1[1]=j,indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				MfindZone2Variable2d(variableTag,facePtr,((zone->grid)->DummySlabV[0]+targetIndices),ptrZone,indices2,targetIndices);
			}
		break;
		case 2:/*Jmax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				
				indices1[0]=i,indices1[1]=*(ptrRange+4),indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				
				MfindZone2Variable2d(variableTag,facePtr,((zone->grid)->DummySlabV[3]+targetIndices),ptrZone,indices2,targetIndices);
			}
		break;
		case 3:/*Jmin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=*(ptrRange+1)-1,indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				MfindZone2Variable2d(variableTag,facePtr,((zone->grid)->DummySlabV[2]+targetIndices),ptrZone,indices2,targetIndices);
			}
		break;
	}
	return 0;
}

#include "APKFillerInterface.h"

int FillerInterface(Interface *facePtr,Zone *ptrZone,int iface)
{
	int i,j,k,indices1[3],indices2[3];
	int *ptrRange;
	int totalCells;
	double variable;
	 
	ptrRange=facePtr->Point_Range;
	//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
	switch(iface)
	{
		case 0:/* Imax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			facePtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+3),indices1[1]=j,indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				//printf("indices2= %4d%4d%4d\n",indices2[0],indices2[1],indices2[2]);
				variable=findZone2Volumes(ptrZone,indices2);
				*(facePtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1)))= variable;      
			}
		break;
		case 1:/*Imin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			facePtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+0)-1,indices1[1]=j,indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				//printf("indices2= %4d%4d%4d\n",indices2[0],indices2[1],indices2[2]);
				variable=findZone2Volumes(ptrZone,indices2);
				*(facePtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1)))= variable;      
			}
		break;
		case 2:/*Jmax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			facePtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=*(ptrRange+4),indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				//printf("indices2= %4d%4d%4d\n",indices2[0],indices2[1],indices2[2]);
				variable=findZone2Volumes(ptrZone,indices2);
				*(facePtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))= variable;      
			}
		break;
		case 3:/*Jmin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			facePtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=*(ptrRange+1)-1,indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				//printf("indices2= %4d%4d%4d\n",indices2[0],indices2[1],indices2[2]);
				variable=findZone2Volumes(ptrZone,indices2);
				*(facePtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))= variable;      
			}
		break;
		case 4:/*Kmax*/
			totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			facePtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+5);
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				//printf("indices2= %4d%4d%4d\n",indices2[0],indices2[1],indices2[2]);
				variable=findZone2Volumes(ptrZone,indices2);
				*(facePtr->ZSlab->Volume+(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))=  variable;      
			}
		break;
		case 5:/*Kmin*/
			totalCells=(*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			facePtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+2)-1;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				//printf("indices2= %4d%4d%4d\n",indices2[0],indices2[1],indices2[2]);
				variable=findZone2Volumes(ptrZone,indices2);
				*(facePtr->ZSlab->Volume+(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))= variable;      
			}
		break;
	}
	return 0;
}

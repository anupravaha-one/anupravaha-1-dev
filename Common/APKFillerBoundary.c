/**
	This function makes  fictitious cell volumes of zone boundary slabs to zeroes. This is the assumption in our solver that
	fictitious boundary slab volumes are zeroes.
*/
#include "APKFillerBoundary.h"

int FillerBoundary(ZoneBoundary *bcPtr,int iface)
{
	int i,j,k;
	int *ptrRange;
	int totalCells;
	ptrRange=bcPtr->Point_Range;
	//printf(" \t\t%d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
	switch(iface)
	{
		case 0:/* Imax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			
			bcPtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
						
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
				*(bcPtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1)))=0.0;      
		break;
		case 1:/*Imin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			bcPtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
				*(bcPtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1)))=0.0;      
		break;
		case 2:/*Jmax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			bcPtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
				*(bcPtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))=0.0;      
		break;
		case 3:/*Jmin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			bcPtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
				*(bcPtr->ZSlab->Volume+(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))=0.0;      
		break;
		case 4:/*Kmax*/
			totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			bcPtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
				*(bcPtr->ZSlab->Volume+(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))=0.0;      
		break;
		case 5:/*Kmin*/
			totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			bcPtr->ZSlab->Volume=(double *)malloc(totalCells*sizeof(double));
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
				*(bcPtr->ZSlab->Volume+(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0)))=0.0;      
		break;
	}
	return 0;
}

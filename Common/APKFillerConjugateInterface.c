#include "APKFillerConjugateInterface.h"

int FillerConjugateInterface(int variableTag,Interface *facePtr,int iface)
{
	int i,j,k,targetIndices;
	int *ptrRange;
	int totalCells;
	double *slabPointer;//,*slabPointerMass;
	ptrRange=facePtr->Point_Range;
	//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
	slabPointer=FindSlabPtrIface(facePtr,variableTag);
	//slabPointerMass=FindSlabPtrIface(facePtr,variableTag*100);
	
	switch(iface)
	{
		case 0:/* Imax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				*(slabPointer+targetIndices)= 0.0;
				//*(slabPointerMass+targetIndices)= 0.0;
			}
		break;
		case 1:/*Imin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				*(slabPointer+targetIndices)= 0.0;
				//*(slabPointerMass+targetIndices)= 0.0;
			}
		break;
		case 2:/*Jmax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices)= 0.0;
				//*(slabPointerMass+targetIndices)= 0.0;
			}
		break;
		case 3:/*Jmin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices)= 0.0;
				//*(slabPointerMass+targetIndices)= 0.0;
			}
		break;
		case 4:/*Kmax*/
			totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices)= 0.0;
				//*(slabPointerMass+targetIndices)= 0.0;
			}
		break;
		case 5:/*Kmin*/
			totalCells=(*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices)= 0.0;
				//*(slabPointerMass+targetIndices)= 0.0;
			}
		break;
	}
	
	facePtr->ConjugateRms = 0;//<<<< Narsimha
	return 0;
}

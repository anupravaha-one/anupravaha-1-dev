#include "APKFillerIfaceConjugatePressure.h"

int FillerIfaceConjugatePressure(Zone *zone,Interface *facePtr,int iface,int vartag)
{
	//printf("i am in pressure update\n");
	//printf("iface =%d\n",iface);
	int i,j,k,indexP,targetIndices;
	int ni,nj,nk;
	int *ptrRange;
	double *slabPointer,*domainPointer;
	ptrRange=facePtr->Point_Range;
	ni=*(zone->Size+0);
	nj=*(zone->Size+1);
	nk=*(zone->Size+2);
	
	slabPointer=FindSlabPtrIface(facePtr,vartag);
	domainPointer=MakeDomainArray(zone,vartag);
	//printf("iface =%d\n",iface);
	switch (iface)
	{
		//printf("iface =%d\n",iface);
		
		case 0 :/* Imax*/
			//printf("i am in pressure update 0\n");
			
			for (k=0;k<nk-1;k++) 
			for (j=0;j<nj-1;j++)
			{
				//printf("i am in case 0\n");
				i=ni-2;
				//targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				targetIndices = (nj-1)*k+j;
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				//*((facePtr->ZSlab->P)+targetIndices)=*((zone->Flowsol->P)+indexP);
				///printf("0 Slab Pressure %lf \n",*(slabPointer+targetIndices));
			}
		break;

		case 1 :/*Imin*/
		//printf("i am in pressure update 1\n");
		
			for (k=0;k<nk-1;k++) 
			for (j=0;j<nj-1;j++)
			{
				i=0;
				//printf("i am in case 1\n");
				//targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				targetIndices = (nj-1)*k+j;
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				//*((facePtr->ZSlab->P)+targetIndices)=*((zone->Flowsol->P)+indexP);
				//printf("1 Slab Pressure %lf \n",*(slabPointer+targetIndices));
			}
		break;

		case 2 :/*Jmax*/
		//printf("i am in pressure update 2\n");
		
			for (k=0;k<nk-1;k++) 
			for (i=0;i<ni-1;i++)
			{
				j=nj-2;
				//printf("i am in case 2\n");
				//targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				targetIndices = (ni-1)*k+i;
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				//*((facePtr->ZSlab->P)+targetIndices)=*((zone->Flowsol->P)+indexP);
				//printf("2 Slab Pressure %12.8e \n",*(slabPointer+targetIndices));
			}
		break;

		case 3 :/*Jmin*/
		//printf("i am in pressure update 3\n");
		
			for (k=0;k<nk-1;k++) 
			for (i=0;i<ni-1;i++)
			{
				j=0;
				//printf("i am in case 3\n");
				//targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				targetIndices = (ni-1)*k+i;
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				//*((facePtr->ZSlab->P)+targetIndices)=*((zone->Flowsol->P)+indexP);
				//printf("3 Slab Pressure %12.8e \n",*(slabPointer+targetIndices));
			}
		break;

		case 4 :/*Kmax*/
		//printf("i am in pressure update 4\n");
		
			for (j=0;j<nj-1;j++) 
			for (i=0;i<ni-1;i++)
			{
				k=nk-2;
				//printf("i am in case 4\n");
				//targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				targetIndices = (ni-1)*j+i;
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				//*((facePtr->ZSlab->P)+targetIndices)=*((zone->Flowsol->P)+indexP);
				//printf("4 Slab Pressure %lf \n",*(slabPointer+targetIndices));
			}
		break;

		case 5 :/*Kmin*/
		//printf("i am in pressure update 5\n");
		
			for (j=0;j<nj-1;j++) 
			for (i=0;i<ni-1;i++)
			{
				k=0;
				//printf("i am in case 5\n");
				//targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				targetIndices = (ni-1)*j+i;
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				//*((facePtr->ZSlab->P)+targetIndices)=*((zone->Flowsol->P)+indexP);
				//printf("5 Slab Pressure %lf \n",*(slabPointer+targetIndices));
			}
		break;
	}
	
	facePtr->ConjugateRms = 0.0;
	return (0);
}

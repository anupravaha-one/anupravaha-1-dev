#include "APKFillerInterfaceRadiationVariables.h"

/** This function put the interface cell center value to Interface patch*/

int FillerInterfaceRadiationVariables(Base *base, int variableTag,Interface *facePtr,Zone *ptrZone,int iface)
{
	int nt,np,NI,NJ,NK,i,j,k,I,J,K,l,m,indices1[3],indices2[3],targetIndices,*sizePtrNeighbor,*ptrRange,totalCells;
	double *slabPointer,*domainPointer;
	
	ptrRange=facePtr->Point_Range;
	nt=base->BaseParameters->ntheta;
	np=base->BaseParameters->nphi;
	
	slabPointer=FindSlabPtrIface(facePtr,variableTag);
	domainPointer=MakeDomainArray(ptrZone,variableTag);
	
	sizePtrNeighbor=ptrZone->Size;
	NI=*(sizePtrNeighbor+0);
	NJ=*(sizePtrNeighbor+1);
	NK=*(sizePtrNeighbor+2);
	
	switch(iface)
	{
		case 1:/* Imax*/
			
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			
	for(l=0;l<nt;l++)
	for(m=0;m<np;m++)
		for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
		for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
		{
			indices1[0]=*(ptrRange+3),indices1[1]=j,indices1[2]=k;
			
			TransformCellIndicesCommon(indices1,facePtr,indices2);
			I=indices2[0]-1;
			J=indices2[1]-1;
			K=indices2[2]-1;
			
			targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
			
			*(slabPointer+targetIndices*nt*np+np*l+m)= *(domainPointer+((NI-1)*(NJ-1)*K+(NI-1)*J+I)*nt*np+np*l+m);
				
			}
		break;
		
		case 2:/*Imin*/
		
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			
		for(l=0;l<nt;l++)
		for(m=0;m<np;m++)
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+0)-1,indices1[1]=j,indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				I=indices2[0]-1;
				J=indices2[1]-1;
				K=indices2[2]-1;
				
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				
				*(slabPointer+targetIndices*nt*np+np*l+m)= *(domainPointer+((NI-1)*(NJ-1)*K+(NI-1)*J+I)*nt*np+np*l+m);
			}
		break;
		
		case 3:/*Jmax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
		
		for(l=0;l<nt;l++)
		for(m=0;m<np;m++)	
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				
				indices1[0]=i,indices1[1]=*(ptrRange+4),indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				I=indices2[0]-1;
				J=indices2[1]-1;
				K=indices2[2]-1;
				
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				
				*(slabPointer+targetIndices*nt*np+np*l+m)= *(domainPointer+((NI-1)*(NJ-1)*K+(NI-1)*J+I)*nt*np+np*l+m);
			}
		break;
		
		case 4:/*Jmin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			
		for(l=0;l<nt;l++)
		for(m=0;m<np;m++)
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=*(ptrRange+1)-1,indices1[2]=k;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				I=indices2[0]-1;
				J=indices2[1]-1;
				K=indices2[2]-1;
				
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices*nt*np+np*l+m)= *(domainPointer+((NI-1)*(NJ-1)*K+(NI-1)*J+I)*nt*np+np*l+m);
			}
		break;
		case 5:/*Kmax*/
			totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
		
		for(l=0;l<nt;l++)
		for(m=0;m<np;m++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+5);
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				I=indices2[0]-1;
				J=indices2[1]-1;
				K=indices2[2]-1;
				
				targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices*nt*np+np*l+m)= *(domainPointer+((NI-1)*(NJ-1)*K+(NI-1)*J+I)*nt*np+np*l+m);
			}
		break;
		case 6:/*Kmin*/
			totalCells=(*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0));
		
		for(l=0;l<nt;l++)
		for(m=0;m<np;m++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+2)-1;
				
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				I=indices2[0]-1;
				J=indices2[1]-1;
				K=indices2[2]-1;
				
				targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				*(slabPointer+targetIndices*nt*np+np*l+m)= *(domainPointer+((NI-1)*(NJ-1)*K+(NI-1)*J+I)*nt*np+np*l+m);
			}
		break;
	}
	return 0;
}

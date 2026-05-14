#include "APKupdateBCDeltaInverse.h"

/**It will store value of Deltainverse(geometry parameters), update the interface and boundary slab value for Deltainverse in variable Delta.Based on whether domain is 2d and 3d it will store value*/


int updateBCDeltaInverse(int c,Zone *zone,ZoneBoundary *bcPtr,int dir)
{
	int k,j,i,NI,NJ,NK,m;
	double *variable,*domain;
	int *ptrRange,index;
	
	
	ptrRange = bcPtr->Point_Range;
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	
	
	variable=bcPtr->ZSlab->Delta;
	
	domain=zone->ZGeom->DeltaInverse;
		switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					for(m=0;m<c;m++)	
					*(variable+c*index+m)=*(domain+c*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+m);
			
					}
					
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						for(m=0;m<c;m++)	
					*(variable+c*index+m)=*(domain+c*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+m);
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						for(m=0;m<c;m++)	
					*(variable+c*index+m)=*(domain+c*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+m);
						
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						for(m=0;m<c;m++)	
					*(variable+c*index+m)=*(domain+c*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+m);	
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					for(m=0;m<c;m++)	
					*(variable+c*index+m)=*(domain+c*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+m);
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						
					for(m=0;m<c;m++)	
					*(variable+c*index+m)=*(domain+c*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+m);
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from updateBCDeltaInverse\n");
					exit(1);
				break;
			}
		
			return 0;
		}
		
		


int FillerInterfaceDelta(int c,Interface *facePtr,Zone *ptrZone,int iface)
{
	int i,j,k,indices1[3],indices2[3],targetIndices;
	int *ptrRange;
	int totalCells;
	
	ptrRange=facePtr->Point_Range;
	//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
	
	switch(iface)
	{
		case 0:/* Imax*/
			
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+3),indices1[1]=j,indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				findInterfaceZoneDelta(c,facePtr,ptrZone,indices2,targetIndices);
			}
		break;
		case 1:/*Imin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			{
				indices1[0]=*(ptrRange+0)-1,indices1[1]=j,indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
				findInterfaceZoneDelta(c,facePtr,ptrZone,indices2,targetIndices);
			}
		break;
		case 2:/*Jmax*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				
				indices1[0]=i,indices1[1]=*(ptrRange+4),indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				
				findInterfaceZoneDelta(c,facePtr,ptrZone,indices2,targetIndices);
			}
		break;
		case 3:/*Jmin*/
			totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=*(ptrRange+1)-1,indices1[2]=k;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				findInterfaceZoneDelta(c,facePtr,ptrZone,indices2,targetIndices);
			}
		break;
		case 4:
			totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+5);
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				findInterfaceZoneDelta(c,facePtr,ptrZone,indices2,targetIndices);
			}
		break;
		case 5:
			totalCells=(*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0));
			//printf("\t\t\ttotalCells=%d\n",totalCells);
			for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
			for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
			{
				indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+2)-1;
				//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
				TransformCellIndicesCommon(indices1,facePtr,indices2);
				targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				findInterfaceZoneDelta(c,facePtr,ptrZone,indices2,targetIndices);
			}
		break;
	}
	return 0;
}
			
int findInterfaceZoneDelta(int c,Interface *facePtr,Zone *zone,int index2[3],int indexofSlabs)
{
	int k,j,i,NI,NJ,NK,m;
	int *sizePtrNeighbor;
	double *slabPointer,*domainPointer;
	
	sizePtrNeighbor=zone->Size;
	i=index2[0]-1;
	j=index2[1]-1;
	k=index2[2]-1;
	
	NI=*(sizePtrNeighbor+0);
	NJ=*(sizePtrNeighbor+1);
	NK=*(sizePtrNeighbor+2);
	
	slabPointer=facePtr->ZSlab->Delta;
	 
	domainPointer=zone->ZGeom->DeltaInverse;
	for(m=0;m<c;m++)
	*(slabPointer+c*indexofSlabs+m)= *(domainPointer+c*((NI-1)*(NJ-1)*k+(NI-1)*j+i)+m);
	
	
	return 0;
}













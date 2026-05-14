#include "APKCalcIFACE_data.h"
/**

This function calculates the indices of the cells which are adjacent to the interface in the language of the neighbouring zone.  These calculated indices are stored in to the zone->Slab data structure.  While transferring the data across the interface we can directly use these stored indices rather than finding them again.

Modified interface updation:
Every time for filling the interfaces, the transformcellindices function is called to find out the target indices of the adjacent zone.  This can be avoided by just adding a slab for indices in the zone slab data structure.  This has improved the speed of Filler interface variables.  This change has forced a change to be made in ExtractVariables.c such that mass velocities gets updated in interfaces.  It is kept in mind that Boundary conditions for massvelocites need not be applied in ExtractVariables.c as we call separate function in findresidual*.c

This function is added by BHANS.

Refer Profiling and Optimization manual for further details.
**/

int CalcIFACE_data(Base * base)
{
	int i,j,k,indices1[3],indices2[3],targetIndices,counter;
	int *ptrRange,*Neighbour;
	int totalCells,iface;
	Interface *facePtr;
	ZoneBoundary *ZoneBc;
			
	Zone *ptrZone,*zone;
	
	int numZones,tag,NI,NJ,NK;
	
	
	for(numZones=0;numZones<base->Number_of_zones;numZones++)
	{
		zone = base->CgnsZone+numZones;
		for(iface=0;iface<6;iface++)
		{
			tag = zone->BC_IFACE_TAG[iface];
			
			if(tag==2)
			{
				facePtr = zone->ZONE_FACES[iface];
				ptrZone=base->CgnsZone+zone->ZONE_NEIGHBOUR_NUM[iface];
				
				ptrRange=facePtr->Point_Range;
				
				Neighbour =ptrZone->Size;
				NI=*(Neighbour+0);
				NJ=*(Neighbour+1);
				NK=*(Neighbour+2); 
				//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
				printf("IFACE = %d\n",iface);
				//exit(1);
				switch(iface)
				{
					case 0:/* Imax*/
						
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
						//printf("\t\t\ttotalCells=%d\n",totalCells);
						//printf("counter=%d\n",counter);
						facePtr->ZSlab->Total_Cells = totalCells;
						
						facePtr->ZSlab->Target_Indices = (int *)malloc(totalCells*sizeof(int));
						
						counter=0;
						for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
						for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
						{
							indices1[0]=*(ptrRange+3),indices1[1]=j,indices1[2]=k;
							//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
							//printf("counter=%d\n",counter);
							TransformCellIndicesCommon(indices1,facePtr,indices2);
							targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
							*(facePtr->ZSlab->Target_Indices+counter) = (NI-1)*(NJ-1)*(indices2[2]-1)+(NI-1)*(indices2[1]-1)+(indices2[0]-1);
							counter++;
							
							//findZone2Variable(variableTag,facePtr,ptrZone,indices2,targetIndices);
						}
					break;
					case 1:/*Imin*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
						//printf("\t\t\ttotalCells=%d\n",totalCells);
						
						facePtr->ZSlab->Total_Cells = totalCells;
						facePtr->ZSlab->Target_Indices = (int *)malloc(totalCells*sizeof(int));
						counter=0;
						
						for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
						for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
						{
							indices1[0]=*(ptrRange+0)-1,indices1[1]=j,indices1[2]=k;
							//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
							TransformCellIndicesCommon(indices1,facePtr,indices2);
							targetIndices=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
							*(facePtr->ZSlab->Target_Indices+counter) = (NI-1)*(NJ-1)*(indices2[2]-1)+(NI-1)*(indices2[1]-1)+(indices2[0]-1);
							counter++;
							//findZone2Variable(variableTag,facePtr,ptrZone,indices2,targetIndices);
						}
					break;
					case 2:/*Jmax*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
						//printf("\t\t\ttotalCells=%d\n",totalCells);
						
						facePtr->ZSlab->Total_Cells = totalCells;
						facePtr->ZSlab->Target_Indices = (int *)malloc(totalCells*sizeof(int));
						counter=0;
						for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
						for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
						{
							
							indices1[0]=i,indices1[1]=*(ptrRange+4),indices1[2]=k;
							//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
							TransformCellIndicesCommon(indices1,facePtr,indices2);
							targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							*(facePtr->ZSlab->Target_Indices+counter) = (NI-1)*(NJ-1)*(indices2[2]-1)+(NI-1)*(indices2[1]-1)+(indices2[0]-1);
							counter++;
							//findZone2Variable(variableTag,facePtr,ptrZone,indices2,targetIndices);
						}
					break;
					case 3:/*Jmin*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
						//printf("\t\t\ttotalCells=%d\n",totalCells);
						
						facePtr->ZSlab->Total_Cells = totalCells;
						facePtr->ZSlab->Target_Indices = (int *)malloc(totalCells*sizeof(int));
						counter=0;
						for(k= *(ptrRange+2);k < *(ptrRange+5);k++)
						for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
						{
							indices1[0]=i,indices1[1]=*(ptrRange+1)-1,indices1[2]=k;
							//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
							TransformCellIndicesCommon(indices1,facePtr,indices2);
							targetIndices=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							*(facePtr->ZSlab->Target_Indices+counter) = (NI-1)*(NJ-1)*(indices2[2]-1)+(NI-1)*(indices2[1]-1)+(indices2[0]-1);
							counter++;
							//findZone2Variable(variableTag,facePtr,ptrZone,indices2,targetIndices);
						}
					break;
					case 4:/*Kmax*/
						totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
						//printf("\t\t\ttotalCells=%d\n",totalCells);
						
						facePtr->ZSlab->Total_Cells = totalCells;
						facePtr->ZSlab->Target_Indices = (int *)malloc(totalCells*sizeof(int));
						counter=0;
						for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
						for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
						{
							indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+5);
							//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
							TransformCellIndicesCommon(indices1,facePtr,indices2);
							targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							*(facePtr->ZSlab->Target_Indices+counter) = (NI-1)*(NJ-1)*(indices2[2]-1)+(NI-1)*(indices2[1]-1)+(indices2[0]-1);
							counter++;
							//findZone2Variable(variableTag,facePtr,ptrZone,indices2,targetIndices);
						}
					break;
					case 5:/*Kmin*/
						totalCells=(*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0));
						//printf("\t\t\ttotalCells=%d\n",totalCells);
						
						facePtr->ZSlab->Total_Cells = totalCells;
						facePtr->ZSlab->Target_Indices = (int *)malloc(totalCells*sizeof(int));
						counter=0;
						for(j= *(ptrRange+1);j < *(ptrRange+4);j++)
						for(i= *(ptrRange+0);i < *(ptrRange+3);i++)
						{
							indices1[0]=i,indices1[1]=j,indices1[2]=*(ptrRange+2)-1;
							//printf("indices1= %4d%4d%4d\n",indices1[0],indices1[1],indices1[2]);
							TransformCellIndicesCommon(indices1,facePtr,indices2);
							targetIndices=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							*(facePtr->ZSlab->Target_Indices+counter) = (NI-1)*(NJ-1)*(indices2[2]-1)+(NI-1)*(indices2[1]-1)+(indices2[0]-1);
							counter++;
							//findZone2Variable(variableTag,facePtr,ptrZone,indices2,targetIndices);
						}
					break;
					
					default:
						printf("Direction not found in CalcIFACE_data.c\nExiting......\n");
						exit(12);
					break;
						
				}
			}
			else
			{
				ZoneBc = zone->ZONE_FACES[iface];
				
				ptrRange=ZoneBc->Point_Range;
				
				switch(iface)
				{
					case 0:/* Imax*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
							
						ZoneBc->ZSlab->Total_Cells = totalCells;
						ZoneBc->ZSlab->Target_Indices = NULL;
					break;
					case 1:/*Imin*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1)));
							
						ZoneBc->ZSlab->Total_Cells = totalCells;
						ZoneBc->ZSlab->Target_Indices = NULL;
					break;
					case 2:/*Jmax*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
	
						ZoneBc->ZSlab->Total_Cells = totalCells;
						ZoneBc->ZSlab->Target_Indices = NULL;
					break;
					case 3:/*Jmin*/
						totalCells=((*(ptrRange+5)-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0)));
						
						ZoneBc->ZSlab->Total_Cells = totalCells;
						ZoneBc->ZSlab->Target_Indices = NULL;
					break;
					case 4:/*Kmax*/
						totalCells=((*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0)));
			
						ZoneBc->ZSlab->Total_Cells = totalCells;
						ZoneBc->ZSlab->Target_Indices = NULL;
						
					break;
					case 5:/*Kmin*/
						totalCells=(*(ptrRange+4)-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0));
			
						ZoneBc->ZSlab->Total_Cells = totalCells;
						ZoneBc->ZSlab->Target_Indices = NULL;
					break;
					default:
						printf("Direction not found in CalcIFACE_data.c\nExiting......\n");
						exit(12);
					break;
				}
			}
		}
	}
	return 0;
}

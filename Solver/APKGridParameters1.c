#include "APKGridParameters1.h"

/**
This function is used for Multigrid method. It calculates the coarsest grid index (CGI) based on the least value among
NI,NJ,NK of all the zones. The least value is halved till a grid size of just less than 8 is reached. All the zones are
assigned with the same CGI because all the zones should necessarily have the same number of grid levels.
*/

int GridParameters1(Root CgnsRoot)
{
	Base *base;
	base=CgnsRoot.CgnsBase;
	Zone *zone;
	
	int i,NI,NJ,NK,N,VN,CGIndex;
	int *ptrSize,iZone,numZones;
	
	numZones=base->Number_of_zones;

	CGIndex = 0;
	
	if(base->BaseParameters->flag2d == 1)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			ptrSize = zone->Size ;
			
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			
			if (iZone == 0)
				N = NI;
				
			if(N >= NJ)
				N = NJ;
				
			if(N >= NI)
				N = NI;
		}
	}
	else if(base->BaseParameters->flag2d == 0)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			ptrSize = zone->Size ;
			
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			
			if (iZone == 0)
				N = NI;
				
			if(N >= NJ)
				N = NJ;
				
			if(N >= NI)
				N = NI;
				
			if(N >= NK)
				N = NK;
		}
	}
	
	do
	{
		VN=N/2;
		N=VN;
		CGIndex = CGIndex+1;
	}while(VN>8);
	
	base->CoarseGridIndex = CGIndex;
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		
		zone->zonalCoarseGridIndex = CGIndex;
	}
	
	printf("coarsest grid index is %d\n",CGIndex);
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		zone->grid=(Grid *)malloc(((zone->zonalCoarseGridIndex)+1)*sizeof(Grid));
		if((zone->grid)==NULL)
		{
			printf("Available memory is not sufficient\n");
			exit(1);
		}
		
		for(i=0;i<=(zone->zonalCoarseGridIndex);i++)
		{
			(zone->grid+i)->GSize = (int *)malloc(3*sizeof(int));
			if((zone->grid+i)->GSize == NULL)
			{
				printf("Available memory is not sufficient\n");
				exit(1);
			}
		}
		for(i=0;i<3;i++)
		{
			*(zone->grid->GSize+i) = *(zone->Size+i);
		}
	}
	
	if(base->BaseParameters->flag2d == 1)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			printf("Zone = %d\n",iZone);
			
			printf("Grid = 0... %d %d %d\n",*(zone->grid->GSize+0),*(zone->grid->GSize+1),*(zone->grid->GSize+2));

			for(i=1;i<=(zone->zonalCoarseGridIndex);i++)
			{
				CalcGridSize2d(zone,i);
				printf("Grid = %d... %d %d %d\n",i,*((zone->grid+i)->GSize+0),*((zone->grid+i)->GSize+1),*((zone->grid+i)->GSize+2));
			}
		}
	}
	else if(base->BaseParameters->flag2d == 0)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;

			printf("Zone = %d\n",iZone);
			
			printf("Grid = 0... %d %d %d\n",*(zone->grid->GSize+0),*(zone->grid->GSize+1),*(zone->grid->GSize+2));

			for(i=1;i<=(zone->zonalCoarseGridIndex);i++)
			{
				CalcGridSize3d(zone,i);
				printf("Grid = %d... %d %d %d\n",i,*((zone->grid+i)->GSize+0),*((zone->grid+i)->GSize+1),*((zone->grid+i)->GSize+2));
			}
		}
	}
	return 0;
}


int CalcGridSize3d(Zone *zone,int g)
{
	int NI,NJ,NK;
	
	NI = *((zone->grid+g-1)->GSize+0);
	NJ = *((zone->grid+g-1)->GSize+1);
	NK = *((zone->grid+g-1)->GSize+2);
	
	if((NI-1)%2 == 0)
		*((zone->grid+g)->GSize+0) = (NI-1)/2 + 1;
	else
		*((zone->grid+g)->GSize+0) = (NI-2)/2 + 1 + 1;
	
	if((NJ-1)%2 == 0)
		*((zone->grid+g)->GSize+1) = (NJ-1)/2 + 1;
	else
		*((zone->grid+g)->GSize+1) = (NJ-2)/2 + 1 + 1;
		
	if((NK-1)%2 == 0)
		*((zone->grid+g)->GSize+2) = (NK-1)/2 + 1;
	else
		*((zone->grid+g)->GSize+2) = (NK-2)/2 + 1 + 1;
		
	return 0;
}


int CalcGridSize2d(Zone *zone,int g)
{
	int NI,NJ,NK;
	
	NI = *((zone->grid+g-1)->GSize+0);
	NJ = *((zone->grid+g-1)->GSize+1);
	NK = *((zone->grid+g-1)->GSize+2);
	
	if((NI-1)%2 == 0)
		*((zone->grid+g)->GSize+0) = (NI-1)/2 + 1;
	else
		*((zone->grid+g)->GSize+0) = (NI-2)/2 + 1 + 1;
	
	if((NJ-1)%2 == 0)
		*((zone->grid+g)->GSize+1) = (NJ-1)/2 + 1;
	else
		*((zone->grid+g)->GSize+1) = (NJ-2)/2 + 1 + 1;
		
	*((zone->grid+g)->GSize+2) = NK;
	return 0;
}



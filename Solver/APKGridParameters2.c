#include "APKGridParameters2.h"

/**
This function is used for Single grid method i.e., our normal method. It just assigns,
					*(zone->grid->GSize+i) = *(zone->Size+i);
The reason for doing this is that in some functions the "zone size" has been accessed as a grid structure member i.e.,
"zone grid size".
*/

int GridParameters2(Root CgnsRoot)
{
	Base *base;
	Zone *zone;
	int i,iZone,numZones;
	
	base=CgnsRoot.CgnsBase;
	numZones=base->Number_of_zones;
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		zone->grid=(Grid *)malloc(1*sizeof(Grid));
		if((zone->grid)==NULL)
		{
			printf("Available memory is not sufficient, exiting from APKGridParameters2.c\n");
			exit(1);
		}
		
		(zone->grid)->GSize = (int *)malloc(3*sizeof(int));
		if((zone->grid)->GSize == NULL)
		{
			printf("Available memory is not sufficient, exiting from APKGridParameters2.c\n");
			exit(1);
		}
		
		for(i=0;i<3;i++)
		{
			*(zone->grid->GSize+i) = *(zone->Size+i);
		}
	}	
	
	return 0;
}


#include "APKRefreshSurfaceVolumes.h"

/**
This function is called only once in the "main" function. This function allocates memory for the 
zoneslabs of interfaces and boundaries. This function inturn calls "ExtractVolume" which will
store the adjoining zone volumes in the interface slab, and zeroes in the boundary slab (because
at boundary we have a zero volume cell)
*/

int RefreshSurfaceVolumes(Base *base)
{
	int cntZone,numZones,numiface,numBc,zBc,cntface;
	//printf("base=%u\n",base);
	
	numZones=base->Number_of_zones;
	for(cntZone=0;cntZone < numZones;cntZone++)
	{
		/* Allocation of memory for Zone slabs*/
		numiface=(base->CgnsZone+cntZone)->Number_of_Interfaces;
		numBc=(base->CgnsZone+cntZone)->Number_of_BC_Patches;
		
		for(cntface=0;cntface<numiface;cntface++)
			((base->CgnsZone+cntZone)->Iface+cntface)->ZSlab=(ZoneSlab *)malloc(sizeof(ZoneSlab));
			
		for(zBc=0;zBc<numBc;zBc++)
			((base->CgnsZone+cntZone)->ZoneBc+zBc)->ZSlab=(ZoneSlab *)malloc(sizeof(ZoneSlab));
		//printf("\n i m in zone number %d\n",cntZone+1);
		ExtractVolume(base,base->CgnsZone+cntZone);
	}
	return 0;
}

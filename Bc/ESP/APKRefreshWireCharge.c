#include "APKRefreshWireCharge.h"
/**
This function RefreshWireCharge() is similar to RefreshSurfaceCells and is used
for applying boundary condition to update the value of charge density on the wire
This function inturn calls the GetBoundaryFaces() that will actually implement 
BC for the wire boundary.
*/
int RefreshWireCharge(int variableTag,Base *base)
{
	
	int cntZone,numZones;
	
	numZones=base->Number_of_zones;
	for(cntZone=0;cntZone < numZones;cntZone++)
	{
		if ((base->CgnsZone+cntZone)->ZoneWireTag == 0)
			GetBoundaryFaces(variableTag,base,base->CgnsZone+cntZone);
	}
	return 0;
}

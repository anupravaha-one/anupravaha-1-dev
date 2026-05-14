#include "APKGetBoundaryFaces.h"
/**
The function GetBoundaryFaces() is similar to ExtractVariables() and updates
the wire boundary value of charge density. When it encounters wire boundary in
the Zone (wire is specially recognized by giving it a special tag in GridZ, which
is "Wall Viscous Boundary", having a tag "22"). then it calls the function 
Update_Boundary_Charge() which implements the boundary condition for charge density
on wire.
*/
int GetBoundaryFaces(int variableTag,Base *base,Zone *zone)
{
	int iface; 
	char tag,name[33];
	void *ptr;

	ZoneBoundary *bcPtr;
	//Zone *ptrZoneNeighbor;
	
	//Zone *FindZonePtr();
	
	
	for(iface=0;iface < 6;iface++)
	{
		//printf("\tiface=%d\n",iface);
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);//Get the pointer and tag of interface/boundary
		if(tag=='b')
		{
			//printf("\t\tthis is a boundary\n");
			bcPtr=(ZoneBoundary *)ptr;
			if (bcPtr->Boundary_Tag == 22)
				Update_Boundary_Charge(base,variableTag,zone,bcPtr);//Infolder Boundary cond 
		}
	}
	return 0;
}

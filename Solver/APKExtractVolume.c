#include "APKExtractVolume.h"
/**
The function ExtractVolume updates all the interface slab values in the current zone for volume.
It takes base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
Interface pointers, and then with the help of function FillerInterface() it updates the interface 
slabs, when it encounters boundary in the Zone then it calls the function  FillerBoundary().
*/
int ExtractVolume(Base *base,Zone *zone)
{
	int iface,Tag=1;
	char tag,name[33];
	void *ptr;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *ptrZoneNeighbor;
	for(iface=0;iface < 6;iface++)
	{
		//printf("\tiface=%d\n",iface);
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			facePtr=(Interface *)ptr;
			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);/* working properly*/	
			Tag=strcmp((zone->SolidFluidName),(ptrZoneNeighbor->SolidFluidName));
			if(Tag!=0 )
			{
				//printf("i am filling zero volume slab\n");
				FillerConjugateBoundary(facePtr,iface);
			}
			else
				FillerInterface(facePtr,ptrZoneNeighbor,iface);	
		}
		else
		if(tag=='b')
		{
			//printf("\t\tthis is a boundary\n");
			bcPtr=(ZoneBoundary *)ptr;
			FillerBoundary(bcPtr,iface);
		}
	}
	return 0;
}

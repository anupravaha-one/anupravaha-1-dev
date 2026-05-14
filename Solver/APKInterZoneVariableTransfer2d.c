#include "APKInterZoneVariableTransfer2d.h"

/**
The function ExtractVariables updates all the interface slab values in the current zone for the input 
variable, the variables here are Velocities, temperature, Scalars and so on. 
It takes variable tag, Base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
Interface pointers, and then with the help of function FillerInterfaceVariables() it updates the interface 
slabs. It also checks whether the adjoining zone is also of the same type as the current zone, i.e; whether the adjoining
zone is  a solid or fluid if the current zone is a fluid etc, and calls FillerConjugateInterface() to update and implement
interface conditions if it is a conjugate heat transfer case. When it encounters boundary in the Zone then it calls
the function  update_boundary() or update_boundaryTurb(),according to whether it is a laminar or turbulent flow problem,
which implement the boundary conditions.
*/
int InterZoneVariableTransfer2d(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface;
	char tag,name[33];
	void *ptr;
	Interface *facePtr;
	Zone *ptrZoneNeighbor;
	
	for(iface=0;iface < 4;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			facePtr=(Interface *)ptr;
			
			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);
			
			MFillerInterfaceVariables2d(variableTag,facePtr,zone,ptrZoneNeighbor,iface);
		}
	}
	return 0;
}

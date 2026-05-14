#include "APKExtractMassVelocities.h"
/**
The function ExtractMassVelocities updates all the interface slab values in the current zone for the input 
variable, the variables here are Velocities, temperature, Scalars and so on. 
It takes variable tag, Base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
Interface pointers, and then with the help of function FillerInterfaceVariables() it updates the interface 
slabs. It also checks whether the adjoining zone is also of the same type as the current zone, i.e; whether the adjoining
zone is  a solid or fluid if the current zone is a fluid etc, and calls appropriate functions such as HomogeneousDirichletIC(), FillerConjugateInterfaceTurb(), FillerConjugateInterface() etc., to update and implement
interface conditions if it is a conjugate heat transfer case. When it encounters boundary in the Zone then it calls
the function  velocitymassBc() which, in turn, calls update_boundary() or update_boundaryTurb(), according to whether it is a laminar or turbulent flow problem,
which implement the boundary conditions.

NOTE: variableTag supplied to this function is 100 for Umass, 200 for Vmass and 300 for Wmass.
 */

//variableTag: 100->Umass, 200->Vmass and 300->Wmass
int ExtractMassVelocities(int variableTag,Base *base,Zone *zone)
{
	int index,Tag=1,numBC,numIface;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *ptrZoneNeighbor;

	numBC = zone->Number_of_BC_Patches;
	numIface = zone->Number_of_Interfaces;

	if(variableTag%100 != 0)
		variableTag = variableTag * 100;
	
	//Updating the boundaries
	for(index=0;index < numBC;index++)
	{
		bcPtr = zone->ZoneBc+index;
		//printf("%s Dir=%d Implementing BC for vartag=%d\n", zone->zonename,bcPtr->Direction, variableTag);
		velocitymassBc(base,variableTag/100,zone,bcPtr);
	}
	
	//Updating the interfaces
	for(index=0;index < numIface;index++)
	{
		facePtr = zone->Iface+index;
		ptrZoneNeighbor = (Zone *) facePtr->AdjoiningZonePtr;
		//ptrZoneNeighbor = (Zone *)FindZonePtr(facePtr->Adjoining_Zone_Name,base);
		
		Tag=strcmp(zone->SolidFluidName, ptrZoneNeighbor->SolidFluidName);
		
		if(Tag!=0) //Two zones are of different material
		{
			if((zone->flagSolidFluid - ptrZoneNeighbor->flagSolidFluid) != 0)
			{
				if(zone->flagSolidFluid == 0 && (base->flagTurbModel ==1 || base->flagTurbModel ==3))
					FillerConjugateInterfaceTurb(base,variableTag,facePtr,facePtr->Direction,zone);
				else
					HomogeneousDirichletIC(variableTag, facePtr);
			}
			else 
			{
				//for Slip BC Implementation
				FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,facePtr->Direction-1);
			}
		}
		else
		{
			FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,facePtr->Direction-1);
		}
	}
	
	return 0;
}

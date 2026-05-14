#include"APKExtractVariableCoefficients.h"

/**
The function ExtractVariableCoefficients updates all the interface slab values in the current zone for the input 
variable, the variables here are Density, Dynamic viscocity, Thermal Conductivity, Specific heat and so on. 
It takes variable tag, Base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
Interface pointers, and then with the help of function FillerInterfaceVariables() it updates the interface 
slabs, when it encounters boundary in the Zone then it calls the function  updateBCslabsVarCoeff().
*/
int  ExtractVariableCoefficients(int variableTag,Base *base,Zone *zone)
{
	int iface;
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
			
			//printf("\t\tthis is an interface\n");
			facePtr=(Interface *)ptr;
			
			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);/* working properly*/	
			
			FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
		}
		else /*Note: For Coefficients of variables like rho,viscosity etc. there is no boundary condition implementation*/
		if(tag=='b')
		{
			//printf("\t\tthis is a boundary\n");
			bcPtr=(ZoneBoundary *)ptr;
			//printf("updating boundary direction=%d\n",iface);
			updateBCslabsVarCoeff(variableTag,zone,bcPtr,iface+1,base);
		}
			
	}
	return 0;
}

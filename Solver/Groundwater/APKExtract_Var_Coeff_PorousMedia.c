#include"APKExtract_Var_Coeff_PorousMedia.h"


/**This one is similar to ExtractVariables.c .. It fill the interface and boundary value for groundwater variables and diffusion coeffcents as permeability , or dsipserion coeffceints*/

int  Extract_Var_Coeff_PorousMedia(int variableTag,Base *base,Zone *zone)
{
	int iface;
	char tag,name[33];
	void *ptr;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *ptrZoneNeighbor;
			       
	for(iface=0;iface < 6;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{	facePtr=(Interface *)ptr;
			
			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);
			
			FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
		}
		else if(tag=='b')
		{         					
			bcPtr=(ZoneBoundary *)ptr;
			
			if(variableTag==719)
			{ 
			UpdateBCslabsKunsat(base,zone,bcPtr,19,variableTag,iface+1);
			}
			else
	        	updateBCslabsVar_Coeff_PorousMedia(base,variableTag,zone,bcPtr,iface+1);
			
			
		}
			
	}
	
	return 0;
}

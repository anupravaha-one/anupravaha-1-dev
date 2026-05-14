#include "APKsavemassflux.h"

/**
This function updates the old values of massflux with the current values of the mass flux for each
and every face of the cell.
*/

int savemassflux(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k,iFace;
	Zone *zone;
	
	numzones= base->Number_of_zones;
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			for(k=0;k<NK-1;k++)   /*Storing current Massflux in PrevFlowsol Massflux*/
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			for(iFace=0;iFace<6;iFace++)
			{
				index = ((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6+iFace;
				*(zone->PrevFlowsol->Massflux+index) =*(zone->Flowsol->Massflux+index);
			}
		}
	}
	return 0;
}


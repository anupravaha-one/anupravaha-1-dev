#include "APKsavevolumeflux.h"

/**
This function is exactly same as savemassflux.
This function updates the old values of volumeflux with the current values of the volumeflux for each
and every face of the cell.
*/

int savevolumeflux(Base *base)
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
			index=0;
			for(k=0;k<NK-1;k++)   /*Storing current Massflux in PrevFlowsol Massflux*/
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			for(iFace=0;iFace<6;iFace++)
			{
				//index = ((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6+iFace;
				*(zone->PrevFlowsol->Massflux+index) =*(zone->Flowsol->Massflux+index);
				
			       if((base->flagParticle==1)||(base->flagDroplet==1))
				{
				*(zone->PrevFlowsol->MassfluxDisperse+index) =*(zone->Flowsol->MassfluxDisperse+index);
				} 
				index++;
			}
		}
	}
	return 0;
}

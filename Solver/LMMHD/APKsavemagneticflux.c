#include "APKsavemagneticflux.h"


int savemagneticflux(Base *base)//rajendra
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
				*(zone->PrevFlowsol->Magneticflux+index) =*(zone->Flowsol->Magneticflux+index);
				
				
			}
		}
	}
	return 0;
}



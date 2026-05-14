#include "APKDragMPFV.h"
#define EPSILON 1e-9
/**
This function calculates the drag force for the case of disperse phase flow. The drag is calculated and stored in the data structure under the source term. This function is called once 
per time step
*/
double MPFV(Base *base)
{
	
       
	int cntZ,NI,NJ,NK,*ptrSize,numZones,i,j,k;
	Zone *zone;
	//int flagDen,flagMu,flagDia;//Flags
	int totalcells,index;/*,flagrhoDisperse=0;*/
	numZones= base->Number_of_zones;

	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		
		index=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			

			*((zone->Flowsol)->Momentumofphasex[0]+index) = (*((zone->Flowsol)->Momentumx+index))-(*((zone->Flowsol)->Momentumofphasex[1]+index));
			*((zone->Flowsol)->Momentumofphasey[0]+index) = (*((zone->Flowsol)->Momentumy+index))-(*((zone->Flowsol)->Momentumofphasey[1]+index));
			
			if((zone->PrevFlowsol->Momentumofphasez[1])!=NULL)
				(*((zone->Flowsol)->Momentumofphasez[0]+index)) = (*((zone->Flowsol)->Momentumz+index))-(*((zone->Flowsol)->Momentumofphasez[1]+index));
			
			
			*((zone->Flowsol)->VarthetaPhase[0]+index) =1.0- (*((zone->Flowsol)->VarthetaPhase[1]+index));
			
			
			
			index++;
		}	
		
	
        }

	return 0;

}



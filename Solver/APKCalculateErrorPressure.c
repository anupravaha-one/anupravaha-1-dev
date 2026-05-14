#include "APKCalculateErrorPressure.h"
/**
This function CalculateErrorPressure() computes the rms value of the difference between current pressure and
previous pressure values for the entire flow domain.
*/
double CalculateErrorPressure(Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k,totalcells=0;
	double rmsP,maxerror;
	double errorP=0.0;
	Zone *zone;
	
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = totalcells+(NI-1)*(NJ-1)*(NK-1);
	
		for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			errorP = errorP + (*(zone->Flowsol->P+index) - *(zone->PrevFlowsol->P+index))*(*(zone->Flowsol->P+index) - *(zone->PrevFlowsol->P+index));
		}
	}
	
	rmsP = sqrt(errorP/totalcells);
	maxerror=rmsP;
	return(maxerror);
}

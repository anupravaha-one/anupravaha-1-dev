#include "APKEvaluateEnergysource.h"
/**
The function EvaluateEnergysource() loops over all the zones and provides necessary data to the user defined function 
userfunc_Energysource() to compute heat generation term arising the energy equation.
*/
int EvaluateEnergysource(Base *base)
{
	int i,j,k,NI,NJ,NK,numZones,cntZ,*ptrSize;
	Zone *zone;
	double *variable,time;
	int totalcells,index;
	double *coordX,*coordY,*coordZ;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3];
	numZones= base->Number_of_zones;
	time = base->BaseParameters->Time;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if((zone->flagHeatGeneration==1) && (zone->flagHeat_con_var==1))
		{
			
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			totalcells = (NI-1)*(NJ-1)*(NK-1);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Energysource;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;
				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfunc_Energysource(centroid[0],centroid[1],centroid[2],time,zone);
			}
		}
	}
return 0;
}

		

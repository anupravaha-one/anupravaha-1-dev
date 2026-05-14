#include "APKCalcVarDropDia.h"
/**
This function CalcVarDropDia() is called in the beginning of the time step from fuction solverSemiCoupled2d() or solverSemiCoupled() if droplet phase flow is choosen. This fuction computes and stores the droplet diameter by the use of previous time step volume fraction and number density values.  
*/
double CalcVarDropDia(Base *base)
{
	int i,j,k,NI,NJ,NK,numZones,cntZ,*ptrSize,flag;
	Zone *zone;
	int totalcells,index;
	double A,B;
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
		
 		flag=(base->BaseParameters->flagVariableDia);
		  
		A=(3.141592654)*(*(zone->PrevFlowsol->NumDensity+index));
		
		B = 6.0*(*(zone->PrevFlowsol->VarthetaDisperse+index));
		
		*(zone->PrevFlowsol->ParticleDia+index*flag) = pow(B/A,1.0/3.0);
		   
		index++;
		}

	     }
	   

return 0;
}	

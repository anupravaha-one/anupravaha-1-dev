#include "APKEvaporate.h"

double Evaporate(Base *base)
{
	double rho,dia;
        double LembdaSt,Mdotv;
	//double NMV;
	int cntZ,NI,NJ,NK,numZones,i,j,k;
	int *ptrSize,flag,flagrho;
 	Zone *zone;
 	int totalcells,index;/*,flagrhoDisperse=0;*/
 	numZones= base->Number_of_zones;
 	//shash
	double EPSILON;
	EPSILON = base->BaseParameters->EPSILONDrop;
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
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
		
		LembdaSt=(*(zone->PrevFlowsol->EvapoConst+index));
		
		flagrho=base->BaseParameters->flagRhoDisperse;
		flag=base->BaseParameters->flagVariableDia;
                
		rho = *(zone->Flowsol->RhoDisperse+index*flagrho);
		dia = *(zone->PrevFlowsol->ParticleDia+index*flag);  
		Mdotv=LembdaSt*(0.7854)*(rho*dia); 
	        if(*(zone->PrevFlowsol->VarthetaDisperse+index)>EPSILON)
		*(zone->PrevFlowsol->Evaporate+index) = Mdotv*(*(zone->PrevFlowsol->NumDensity+index));
		else
		*(zone->PrevFlowsol->Evaporate+index) = 0.0;
		
		}
	
	}
	return 0;
}




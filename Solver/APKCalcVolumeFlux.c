#include "APKCalcVolumeFlux.h"

/* This function is for writing the volume flux of the gas phase at the inlet and outlet of the channel in the files volfluxinlet.dat and volfluxoutlet.dat.
*/

int CalcVolumeFlux(Base *base)
{
		
	int i,j,k,NI,NJ,NK,numZones,cntZ,*ptrSize,flag;
	Zone *zone;
	int index,indexP;
	double dia,A,B;
	double *varthetaDisperse_f;

	double suminlet=0.0,sumoutlet=0.0;
	FILE *fp1,*fp2;
	
	fp1 = fopen("SOLUTIONS/volfluxinlet.dat","a");
	fp2 = fopen("SOLUTIONS/volfluxoutlet.dat","a");

	numZones= base->Number_of_zones;
 	   
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		indexP= 0;
		k=0;
		indexP  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;

		CalcVolfracFaceValue(base, zone, (zone->PrevFlowsol->VarthetaDisperse));
		varthetaDisperse_f = zone->PrevFlowsol->VarthetaD_face;


		for(j=0;j<NJ-1;j++)
		{
			i=0;//inlet as k is already initialized as 0
			indexP  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index = (indexP*6) + 0; //west face of inlet
			suminlet =suminlet+ *(zone->PrevFlowsol->Massflux+index)*(1-(*(varthetaDisperse_f+index)*base->BaseParameters->InitialVolumeFraction));  // + *(zone->PrevFlowsol->MassfluxDisperse+index); //total volume flux at inlet //
			
			//outlet
			i=NI-2;
			indexP  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index = (indexP*6)+1;//east face of outlet
			sumoutlet =sumoutlet+*(zone->PrevFlowsol->Massflux+index)*(1-(*(varthetaDisperse_f+index)*base->BaseParameters->InitialVolumeFraction));  // + *(zone->PrevFlowsol->MassfluxDisperse+index); //total volume flux at outlet //
		}
		fprintf(fp1,"%lf\n",suminlet);
		fprintf(fp2,"%lf\n",sumoutlet);
	}
	   
	fclose(fp1);
	fclose(fp2);

	return 0;

}



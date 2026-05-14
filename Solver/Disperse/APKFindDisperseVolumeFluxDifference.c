#include "APKFindDisperseVolumeFluxDifference.h"

/**
This function calculates and returns the massflux residual for the disperse phase (using current and previous
iteration values).
*/
double FindDisperseVolumeFluxDifference(Base *base)
{
	int numZones,cntZone,i,j,k,NI,NJ,NK,iFace,index,*ptrSize;
	Zone *zone;
	double *massFluxCurrentgas,*massFluxOldgas;
	double *massFluxCurrentdrop,*massFluxOlddrop;
	double sumDifference, temp1, temp2;

	double VolFrac;

	numZones= base->Number_of_zones;
	int totalCells;
	int flagRho;

	sumDifference=0.0;
	totalCells=base->TotalCellFluid;

	for(cntZone=0;cntZone<numZones;cntZone++)
	{
		zone = (base->CgnsZone+cntZone);
		flagRho =  zone->Flowsol->flagRho;
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			massFluxCurrentgas=zone->Flowsol->Massflux;
			massFluxOldgas=zone->PrevFlowsol->Massflux;

			massFluxCurrentdrop=zone->Flowsol->MassfluxDisperse;
			massFluxOlddrop=zone->PrevFlowsol->MassfluxDisperse;

			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			for(iFace=0;iFace<6;iFace++)
			{
				index = ((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6+iFace;

				VolFrac = *(zone->PrevFlowsol->VarthetaDisperse+index)*base->BaseParameters->InitialVolumeFraction;
				temp1 = (*(massFluxCurrentgas+index)*(*(zone->PrevFlowsol->Rho+index*flagRho)*(1-VolFrac))) + (*(massFluxCurrentdrop+index)*(*(zone->PrevFlowsol->RhoDisperse)*VolFrac)); //

				temp2 = (*(massFluxOldgas+index)*(*(zone->PrevFlowsol->Rho+index*flagRho)*(1-VolFrac))) + (*(massFluxOlddrop+index)*(*(zone->PrevFlowsol->RhoDisperse)*VolFrac)); //

				sumDifference=sumDifference+(temp1-temp2)*(temp1-temp2);
			}
		}
	}
	//printf("\nsum difference = %e",sumDifference);
	if (totalCells==0)
	totalCells=base->TotalCellSolid;

	sumDifference=sqrt(sumDifference/(totalCells*6));
	//printf("\nrms = %e",sumDifference);
	return(sumDifference);
}

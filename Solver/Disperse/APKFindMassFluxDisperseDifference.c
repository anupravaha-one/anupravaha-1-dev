#include "APKFindMassFluxDisperseDifference.h"
/**
This function calculates and returns the massflux residual for the disperse phase (using current and previous
iteration values).
*/
double FindMassFluxDisperseDifference(Base *base)
{
	int numZones,cntZone,i,j,k,NI,NJ,NK,iFace,index,*ptrSize;
	Zone *zone;
	double *massFluxCurrent,*massFluxOld;
	double sumDifference;
	numZones= base->Number_of_zones;
	int totalCells;
	sumDifference=0.0;
	totalCells=base->TotalCellFluid;
	for(cntZone=0;cntZone<numZones;cntZone++)
	{
		zone = (base->CgnsZone+cntZone);
		if(zone->flagSolidFluid==0)
		{/*Done By Ashish for Conjugate Heat Transfer Case*/
			//totalCells=totalCells+CountCells(zone);
			ptrSize = zone->Size;
			massFluxCurrent=zone->Flowsol->MassfluxDisperse;
			massFluxOld=zone->PrevFlowsol->MassfluxDisperse;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			for(iFace=0;iFace<6;iFace++)
			{
				index = ((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6+iFace;
				sumDifference=sumDifference+(*(massFluxCurrent+index)-*(massFluxOld+index))*(*(massFluxCurrent+index)-*(massFluxOld+index));
			}
		}
	}
	sumDifference=sqrt(sumDifference/(totalCells*6));
	return(sumDifference);
}

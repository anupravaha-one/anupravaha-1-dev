#include "APKCalcExtraTermsOfVarDensityMassflux.h"

/**
In variable density case, the pressure Poisson equation has some extra terms which add up to diagonal and to right hand
side of the equation. These extra terms are calculated by this function which are used in UpdatePressureVarDensityMassFlux(). It
is called only for 3D case with Two phase flow boiling and combustion problems with variable density.
*/

double CalcExtraTermsOfVarDensityMassflux(Base *base,Zone *zone)
{
	double sumMassFlux,*massFlux,prevRho,currentRho,currentP;
	double CONSTANT,volumeP,prevT,prevP;
	int NI,NJ,NK,i,j,k,index[3],*ptrSize,flagrho,indexP,iFace;
	double delta_t,*domainRho,*slabRho[6],*newMassFlux;
	
	delta_t = base->BaseParameters->TimeStep;
	ptrSize = zone->Size ;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	flagrho = zone->Flowsol->flagRho;
	
	newMassFlux = (double *)malloc(6*(NI-1)*(NJ-1)*(NK-1)*sizeof(double));
	CalcMassFlux(newMassFlux,zone);
        domainRho=MakeDomainArray(zone,710);/*Density*/
        MakeSlabPtrArray(zone,slabRho,710);//Density


	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		

		currentRho = *(zone->Flowsol->Rho+indexP*flagrho);
		prevRho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
		currentP = *(zone->Flowsol->P+indexP);//Contains pressure of current cell
		prevP = *(zone->PrevFlowsol->P+indexP);//Contains Prevpressure of current cell
		volumeP = *(zone->ZGeom->Volume+indexP);
		if(zone->Flowsol->T!=NULL)
		{
			if(base->flagtypeVariableDensity==1)
			prevT = 273.0;
			else
			prevT = *(zone->PrevFlowsol->T+indexP);
		}
		else
		prevT = 1.0;
		
		massFlux=newMassFlux+indexP*6;
		
		CONSTANT = userEnterconstant(zone,indexP,base);//See userDefinedVarCoefficients.c function
		// A constant linking density to pressure; constant = MolWt./ Unviversal Gas Constant.
		
		sumMassFlux=0.0;
		for(iFace=0;iFace<6;iFace++)
		{
			sumMassFlux=sumMassFlux +(*(massFlux+iFace));
		}
//		sumMassFlux=(*(massFlux+0))+(*(massFlux+1))+(*(massFlux+2))+(*(massFlux+3));
		
		*((zone->grid+0)->GRHS+indexP) = sumMassFlux - (CONSTANT*volumeP/(delta_t*prevT))*prevP;
		
		*(zone->Flowsol->VarDensityCoeff+indexP) = -(CONSTANT*volumeP/(delta_t*prevT));
	}
	
	free(newMassFlux);
	return(0);
}


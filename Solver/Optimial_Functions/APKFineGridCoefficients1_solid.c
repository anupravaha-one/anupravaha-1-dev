#include "APKFineGridCoefficients1_solid.h"


/**
This function is used for calculating unsteady-diffusion equation coefficients in case of BTCS/CN time stepping scheme.
*/

double FineGridCoefficients1_solid(Base *base,Zone *zone,int g,int vartag)
{
	int NI,NJ,NK,i,j,k,indexP,*ptrSize,flagrho,scheme,flagSpheat;
	double *Ddiffdia,*DCoeffN[18],delta_t;
	double *beta,*volumeDomain;
	double rho,convectivecoeff;
	double *CoeffN[18],*diffdia;
	double Volg,Temp;

	
	int Timetag;
	
	Timetag = base->BaseParameters->TimeStepingSchemeTag;
	
	delta_t = base->BaseParameters->TimeStep;
	
	if(vartag==600 || vartag==601)
		delta_t = base->BaseParameters->Turb->DualTimeStep;
	
	diffdia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
	}
	Ddiffdia = (zone->grid+g)->DiffusionCoeffPhi_P;
	for(i=0;i<18;i++)
	{
		DCoeffN[i] = (zone->grid+g)->DiffusionCoeffPhi_nb[i];
	}
	
	
	ptrSize = zone->Size ;
	beta = zone->ZGeom->Beta;
	volumeDomain = zone->ZGeom->Volume;
	
	
	scheme = base->SolverScheme;
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	
	flagrho = zone->Flowsol->flagRho;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	indexP=0;
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		if(vartag==9)
			rho = 1.0;
		else
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
				
		Volg = *(volumeDomain+indexP);
		
		Temp = (delta_t/rho);
		
		if(Timetag==1)
		{
			Temp = Temp/2.0;
		}
		
		*(diffdia + indexP) = convectivecoeff*Volg - Temp*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) =  - Temp*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) =  - Temp*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) =  - Temp*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) =  - Temp*(*(DCoeffN[3]+indexP));
		
		*(CoeffN[4]+indexP) =  - Temp*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) =  - Temp*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) =  - Temp*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) =  - Temp*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) =  - Temp*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) =  - Temp*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = - Temp*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = - Temp*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = - Temp*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = - Temp*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = - Temp*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = - Temp*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = - Temp*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = - Temp*(*(DCoeffN[17]+indexP));
		
		indexP++;
	}
	
	return(0);
}

/*
double Greaterof(double A,double B)
{
	if(A>B)
		return A;
	else
		return B;
}
*/

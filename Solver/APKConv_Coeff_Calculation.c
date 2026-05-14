#include "APKConv_Coeff_Calculation.h"


/**
This function is used for calculating convection-diffusion equation coefficients in case of BTCS time stepping scheme.
*/

double Conv_Coeff_Calculation(Base *base,int g,int vartag)
{
	int NI,NJ,NK,i,j,k,indexP,*ptrSize,scheme,flagSpheat;
	double delta_t;
	double *F;
	double Temp,*CoeffN[18],*diffdia;
	
	Zone *zone;
	int iZone,numZones;
	
	scheme = base->SolverScheme;
	delta_t = base->BaseParameters->TimeStep;
	
	if(vartag==600||vartag==601)///Other vartag in the dual time step loop are to be added here.
		delta_t = base->BaseParameters->Turb->DualTimeStep;//delta_t/10.0;
	
	numZones=base->Number_of_zones;
	
	for (iZone=0;iZone<numZones;iZone++)
	{
		zone = base->CgnsZone+iZone;	
		
		diffdia = (zone->grid+g)->Conv_Coeff_P;
		for(i=0;i<6;i++)
		{
			CoeffN[i] = (zone->grid+g)->Conv_Coeff_nb[i];
		}
	
		ptrSize = zone->Size ;
		
		flagSpheat = zone->Flowsol->flagSpecificHeat;
		
		if(scheme==0)
		{
			if (vartag==9)
				F=zone->Flowsol->Electricflux;
			else if(vartag==62||vartag==63||vartag==64||vartag==35||vartag==65) //||vartag==67
 				F=zone->Flowsol->MassfluxDisperse;
			else if(vartag==69)
				F=zone->Flowsol->MassfluxDisperse;
			else if (vartag==68||vartag==710)
		        	F=zone->Flowsol->Massflux;
			else
				F=zone->Flowsol->Massflux;
			
		}
		else
		{
			if (vartag==9)
				F=zone->Flowsol->Electricflux;
			else if(vartag==62||vartag==63||vartag==64||vartag==65||vartag==35) //||vartag==67
 				F=zone->Flowsol->MassfluxDisperse;
			else if(vartag==69)
				F=zone->Flowsol->MassfluxDisperse;
			else if (vartag==68)
		        	F=zone->Flowsol->Massflux;
			else if(vartag==710)//shash
		                F=zone->Flowsol->Massflux;
			else
				F=zone->PrevFlowsol->Massflux;
		}

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		if(base->BaseParameters->flag2d == 1)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
			
				*(diffdia + indexP) = delta_t*(Temp);
				
				*(CoeffN[0]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+0),0.0));
				
				*(CoeffN[1]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+1),0.0));
				
				*(CoeffN[2]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+2),0.0));
				
				*(CoeffN[3]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+3),0.0));
				
				*(CoeffN[4]+indexP) = 0.0;
			
				*(CoeffN[5]+indexP) = 0.0;
				
				indexP++;
			}
		}
		else
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				Temp = 
				(Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
				*(diffdia + indexP) = delta_t*(Temp);
			
				*(CoeffN[0]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+0),0.0));
				
				*(CoeffN[1]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+1),0.0));
				
				*(CoeffN[2]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+2),0.0));
				
				*(CoeffN[3]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+3),0.0));
				
				*(CoeffN[4]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+4),0.0));
			
				*(CoeffN[5]+indexP) = delta_t*(Greaterof(-*(F+indexP*6+5),0.0));
				
				indexP++;
			}
		}
	}
	return(0);
}





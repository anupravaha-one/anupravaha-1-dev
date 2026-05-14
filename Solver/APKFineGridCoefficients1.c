#include "APKFineGridCoefficients1.h"


/**
This function is used for calculating convection-diffusion equation coefficients in case of BTCS time stepping scheme.
*/

double FineGridCoefficients1(Base *base,Zone *zone,int g,int vartag)
{
	int NI,NJ,NK,i,j,k,index[3],indexP,*ptrSize,flagrho,scheme,flagSpheat,flaginduced;//rajendra
	double *Ddiffdia,*DCoeffN[18],*slabVolume[6],delta_t;
	double *beta,*volumeDomain;
	double rho,*F,convectivecoeff;
	double Temp,varCoefficient,*CoeffN[18],*diffdia;
	double Volg;
	
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
	if(scheme==0)
	{
		if (vartag==9)
			F=zone->Flowsol->Electricflux;
		
		else
			F=zone->Flowsol->Massflux;
	}
	else
	{
		if (vartag==9)
			F=zone->Flowsol->Electricflux;
		
		else
			F=zone->PrevFlowsol->Massflux;
	}
	flagrho = zone->Flowsol->flagRho;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	
	if(base->BaseParameters->flag2d == 1)
	{
//--------------------------------------------- Boundary cells ----------------------------------(1)---------
		k = 0; j = 0; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
		
		//*(diffdia + indexP) = (delta_t/rho)*(*(Ddiffdia + indexP));
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		 
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		
//---------------------------------------------------------------------------------------------(2)---------
		k = 0; j = 0;
		
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
					
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-------------------------------------------------------------------------------------------------(3)---------
		k = 0; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
					
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------------------------------------------------------------------(4)-------
		k = 0; j = NJ-2; i = NI-2; 
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		 
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//----------------------------------------------------------------------------------------------------(5)------
		k = 0; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		 
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//------------------------------------------------------------------------------------------------------(6)-------
		k = 0; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------------------------------------------------------------------------(7)-----
		k = 0; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-------------------------------------------------------------------------------------------------------(8)---
		k = 0; i = 0; j = NJ-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		 
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//--------------------------------------------Interior Cells------------------------------------------------------------(9)------
		k = 0;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -(delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -(delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//---------------------------------------End of Interior Cells------------------------------------------
	}
	else if(base->BaseParameters->flag2d == 0)
	{
//--------------------------------------------- Boundary cells ----------------------------------(1)---------
		k = 0; j = 0; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		 
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//---------------------------------------------------------------------------------------------(2)---------
		k = 0; j = 0;
		
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-------------------------------------------------------------------------------------------------(3)---------
		k = 0; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------------------------------------------------------------------(4)-------
		j = NJ-2;
		for(k=1;k<NK-2;k++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//--------------------------------------------------------------------------------------------------(5)--------
		k = 0; j = NJ-2; i = NI-2; 
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//----------------------------------------------------------------------------------------------------(6)------
		j = NJ-2; i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//------------------------------------------------------------------------------------------------------(7)-------
		k = NK-2; j = 0;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//------------------------------------------------------------------------------------------------------(8)-------
		k = 0; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//------------------------------------------------------------------------------------------------------(9)-------
		j = 0; i = 0;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-----------------------------------------------------------------------------------------------------(10)--------
		j = 0; i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-----------------------------------------------------------------------------------------------------(11)-----
		j = NJ-2; i = 0;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//------------------------------------------------------------------------------------------------------(12)-------
		j = 0; i = 0; k = NK-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//------------------------------------------------------------------------------------------------------(13)----
		k = 0; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------------------------------------------------------------------------(14)-----
		k = NK-2; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------------------------------------------------------------------------(15)-----
		k = 0; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-------------------------------------------------------------------------------------------------------(16)---
		k = 0; i = 0; j = NJ-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//--------------------------------------------------------------------------------------------------------(17)------
		j = 0;
		for(k=1;k<NK-2;k++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
				
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//------------------------------------------------------------------------------------------------------(18)------
		k = 0;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//-------------------------------------------------------------------------------------------------------------(19)---
		i=0;
		for(k=1;k<NK-2;k++)
		{
			for(j=1;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//----------------------------------------------------------------------------------------------------------(20)-----
		k = NK-2; j = NJ-2; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//--------------------------------------------------------------------------------------------------------(21)-----
		k = NK-2; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//-------------------------------------------------------------------------------------------------------(22)------
		k = NK-2;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//-------------------------------------------------------------------------------------------------------(23)------
		k = NK-2; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//---------------------------------------------------------------------------------------------------------(24)-----
		k = NK-2; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------------------------------------------------------------------------(25)-----
		i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			for(j=1;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				
				if(vartag==9||flaginduced==1)//rajendra
					rho = 1.0;
				else
			
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+indexP);
				
				Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
				
				*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
				
				*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
					
				*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
					
				*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
					
				*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
					
				*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
				
				*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
				
				*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
				
				*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
				
				*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
				
				*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
				
				*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
				
				*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
				
				*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
				
				*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
				
				*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
				
				*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
				
				*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
				
				*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
			}
		}
//---------------------------------------------------------------------------------------------------------(26)-------
		k = NK-2; j = NJ-2; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		
		if(vartag==9||flaginduced==1)//rajendra
			rho = 1.0;
		else
			
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+indexP);
		
		Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
		
		*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
		
		*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
			
		*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
			
		*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
			
		*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
			
		*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
		
		*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
		
		*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
		
		*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
		
		*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
		
		*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
		
		*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
		
		*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
		
		*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
		
		*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
		
		*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
		
		*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
		
		*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
		
		*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
//------------------------------------------- End of Boundary cells ----------------------------------------
	
//------------------------------------------- Interior cells ----------------------------------------------(27)------
		for(k=1;k<NK-2;k++)
		for(j=1;j<NJ-2;j++)
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			
			if(vartag==9||flaginduced==1)//rajendra
				rho = 1.0;
			else
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			varCoefficient = CalcvarCoefficient(base,zone,vartag,indexP);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+indexP);
			
			Temp = (Greaterof(*(F+indexP*6+1),0.0)) + (Greaterof(*(F+indexP*6+0),0.0)) + (Greaterof(*(F+indexP*6+3),0.0)) + (Greaterof(*(F+indexP*6+2),0.0)) + (Greaterof(*(F+indexP*6+5),0.0)) + (Greaterof(*(F+indexP*6+4),0.0));
			
			*(diffdia + indexP) = varCoefficient*Volg + delta_t*convectivecoeff*(Temp) - (delta_t/rho)*(*(Ddiffdia + indexP));
			
			*(CoeffN[0]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+0),0.0)) - (delta_t/rho)*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+1),0.0)) - (delta_t/rho)*(*(DCoeffN[1]+indexP));
				
			*(CoeffN[2]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+2),0.0)) - (delta_t/rho)*(*(DCoeffN[2]+indexP));
				
			*(CoeffN[3]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+3),0.0)) - (delta_t/rho)*(*(DCoeffN[3]+indexP));
				
			*(CoeffN[4]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+4),0.0)) - (delta_t/rho)*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*delta_t*(Greaterof(-*(F+indexP*6+5),0.0)) - (delta_t/rho)*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -(delta_t/rho)*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -(delta_t/rho)*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -(delta_t/rho)*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -(delta_t/rho)*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -(delta_t/rho)*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -(delta_t/rho)*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -(delta_t/rho)*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -(delta_t/rho)*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -(delta_t/rho)*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -(delta_t/rho)*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -(delta_t/rho)*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -(delta_t/rho)*(*(DCoeffN[17]+indexP));
		}
//---------------------------------------End of Interior Cells------------------------------------------
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

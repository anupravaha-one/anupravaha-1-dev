#include "APKFineGridDefQVKCoefficients2.h"

/**
This function is used for calculating deffered correction coefficients in case of Crank Nicholson time stepping scheme and
Quick convection scheme.
*/

double FineGridDefQVKCoefficients2(Base *base,Zone *zone,int vartag)
{
	int NI,NJ,NK,i,j,k,index[3],indexP,*ptrSize,flagrho,scheme,flagSpheat;
	double *slabVolume[6],delta_t;
	double *beta,*volumeDomain,*CoeffN[12],*diffdia;
	double Volg,Vol[6],rho,*F,convectivecoeff;
	double Ve,Vw,Vn,Vs,Vt,Vb,Vee,Vww,Vnn,Vss,Vtt,Vbb,alpha;
	
	ptrSize = zone->Size ;
	beta = zone->ZGeom->Beta;
	volumeDomain = zone->ZGeom->Volume;
	delta_t = base->BaseParameters->TimeStep;
	scheme = base->SolverScheme;
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	alpha = base->BaseParameters->RatioCenUP;
	
	if(vartag == 600 || vartag == 601)
		delta_t = base->BaseParameters->Turb->DualTimeStep;
	
	diffdia = zone->DeferredCorrCoeffPhi_P;
	for(i=0;i<12;i++)
	{
		CoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
	}
	
	if(scheme==0)
	{
		F=zone->Flowsol->Massflux;
	}
	else
	{
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
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
//---------------------------------------------------------------------------------------------(2)---------
		k = 0; j = 0;
		
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//-------------------------------------------------------------------------------------------------(3)---------
		k = 0; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//---------------------------------------------------------------------------------------------------(4)-------
		k = 0; j = NJ-2; i = NI-2; 
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
//----------------------------------------------------------------------------------------------------(5)------
		k = 0; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
//------------------------------------------------------------------------------------------------------(6)-------
		k = 0; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//---------------------------------------------------------------------------------------------------------(7)-----
		k = 0; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//-------------------------------------------------------------------------------------------------------(8)---
		k = 0; i = 0; j = NJ-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
//--------------------------------------------Interior Cells------------------------------------------------------------(9)------
//----------------------------------------------------------------------------------------------------Q1----------------------
		k = 0;i = 1;
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
		
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									    +((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									    +((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									    +((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     -((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									     -(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     +(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     -((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									     -(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     +(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     -((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									     -(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     +(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     -((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									     -(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     +(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//--------------------------------------------------------------------------------------------------------Q2----------------------
		k = 0;i = 1;j = 1;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
	
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
//---------------------------------------------------------------------------------------------------------Q3----------------
		k = 0;i = 1;j = NJ-3;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
	
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(slabVolume[2]+(NI-1)*k+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
//---------------------------------------------------------------------------------------------------------Q4----------------
		k = 0;i = NI-3;
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
		
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									    +((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									    +((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									    +((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     -((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									     -(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     +(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     -((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									     -(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     +(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     -((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									     -(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     +(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     -((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									     -(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     +(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//------------------------------------------------------------------------------------------------------------Q5------------------
		k = 0;i = NI-3;j = 1;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
	
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
//----------------------------------------------------------------------------------------------------------Q6---------------
		k = 0;i = NI-3;j = NJ-3;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
	
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(slabVolume[2]+(NI-1)*k+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
//-------------------------------------------------------------------------------------------------------------Q7------------
		k = 0;j = 1;
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
		
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									    +((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									    +((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									    +((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     -((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									     -(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     +(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     -((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									     -(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     +(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     -((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									     -(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     +(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     -((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									     -(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     +(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//--------------------------------------------------------------------------------------------------------Q8-----------------
		k = 0;j = NJ-3;
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
		
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(slabVolume[2]+(NI-1)*k+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									    +((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									    +((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									    +((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     -((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									     -(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     +(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									     -((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									     -(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									     +(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     -((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									     -(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     +(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									     -((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									     -(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									     +(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		}
//-----------------------------------------------------------------------------------------------------------Q9--------------
		k = 0;
		for(j=2;j<NJ-3;j++)
		{
			for(i=2;i<NI-3;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
										+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
										+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
										+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0))));
				//west
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
										-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
										-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
										+(Greaterof(-*(F+indexP*6+0),0.0)));
				//west to west
				*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
				//east
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
										-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
										-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
										+(Greaterof(-*(F+indexP*6+1),0.0)));
				//east to east
				*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
				//south
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
										-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
										-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
										+(Greaterof(-*(F+indexP*6+2),0.0)));
				//south to south
				*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
				//north
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
										-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
										-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
										+(Greaterof(-*(F+indexP*6+3),0.0)));
				//north to north
				*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
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
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//---------------------------------------------------------------------------------------------(2)---------
		k = 0; j = 0;
		
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-------------------------------------------------------------------------------------------------(3)---------
		k = 0; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------(4)-------
		j = NJ-2;
		for(k=1;k<NK-2;k++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				Vt = Volg + Vol[5];
				Vb = Volg + Vol[4];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
										+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
										+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
										+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
										+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
										+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
				
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
					
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
					
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
					
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
					
				*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
				
				*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
			}
		}
//--------------------------------------------------------------------------------------------------(5)--------
		k = 0; j = NJ-2; i = NI-2; 
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//----------------------------------------------------------------------------------------------------(6)------
		j = NJ-2; i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//------------------------------------------------------------------------------------------------------(7)-------
		k = NK-2; j = 0;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//------------------------------------------------------------------------------------------------------(8)-------
		k = 0; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//------------------------------------------------------------------------------------------------------(9)-------
		j = 0; i = 0;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-----------------------------------------------------------------------------------------------------(10)--------
		j = 0; i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-----------------------------------------------------------------------------------------------------(11)-----
		j = NJ-2; i = 0;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//------------------------------------------------------------------------------------------------------(12)-------
		j = 0; i = 0; k = NK-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//------------------------------------------------------------------------------------------------------(13)----
		k = 0; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------------(14)-----
		k = NK-2; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------------(15)-----
		k = 0; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-------------------------------------------------------------------------------------------------------(16)---
		k = 0; i = 0; j = NJ-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//--------------------------------------------------------------------------------------------------------(17)------
		j = 0;
		for(k=1;k<NK-2;k++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				Vt = Volg + Vol[5];
				Vb = Volg + Vol[4];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
										+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
										+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
										+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
										+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
										+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
				
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
					
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
					
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
					
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
					
				*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
				
				*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
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
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				Vt = Volg + Vol[5];
				Vb = Volg + Vol[4];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
										+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
										+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
										+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
										+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
										+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
				
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
					
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
					
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
					
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
					
				*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
				
				*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
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
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				Vt = Volg + Vol[5];
				Vb = Volg + Vol[4];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
										+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
										+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
										+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
										+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
										+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
				
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
					
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
					
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
					
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
					
				*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
				
				*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
			}
		}
//----------------------------------------------------------------------------------------------------------(20)-----
		k = NK-2; j = NJ-2; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//--------------------------------------------------------------------------------------------------------(21)-----
		k = NK-2; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-------------------------------------------------------------------------------------------------------(22)------
		k = NK-2;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				Vt = Volg + Vol[5];
				Vb = Volg + Vol[4];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
										+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
										+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
										+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
										+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
										+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
				
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
					
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
					
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
					
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
					
				*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
				
				*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
			}
		}
//-------------------------------------------------------------------------------------------------------(23)------
		k = NK-2; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//---------------------------------------------------------------------------------------------------------(24)-----
		k = NK-2; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
									+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
									+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
									+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
									+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
									+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
			
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
				
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
				
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
				
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
				
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
			
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------------(25)-----
		i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			for(j=1;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;
			
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				
				Ve = Volg + Vol[1];
				Vw = Volg + Vol[0];
				Vn = Volg + Vol[3];
				Vs = Volg + Vol[2];
				Vt = Volg + Vol[5];
				Vb = Volg + Vol[4];
				
				*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
										+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
										+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
										+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
										+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
										+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
				
				*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
					
				*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
					
				*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
					
				*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
					
				*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
				
				*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
			}
		}
//---------------------------------------------------------------------------------------------------------(26)-------
		k = NK-2; j = NJ-2; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
			
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*((((*(F+indexP*6+1))*Vol[1]/Ve)-(Greaterof(*(F+indexP*6+1),0.0)))
								+(((*(F+indexP*6+0))*Vol[0]/Vw)-(Greaterof(*(F+indexP*6+0),0.0)))
								+(((*(F+indexP*6+3))*Vol[3]/Vn)-(Greaterof(*(F+indexP*6+3),0.0)))
								+(((*(F+indexP*6+2))*Vol[2]/Vs)-(Greaterof(*(F+indexP*6+2),0.0)))
								+(((*(F+indexP*6+5))*Vol[5]/Vt)-(Greaterof(*(F+indexP*6+5),0.0)))
								+(((*(F+indexP*6+4))*Vol[4]/Vb)-(Greaterof(*(F+indexP*6+4),0.0))));
		 
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+0))*Volg/Vw)+(Greaterof(-*(F+indexP*6+0),0.0)));
			
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+1))*Volg/Ve)+(Greaterof(-*(F+indexP*6+1),0.0)));
			
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+2))*Volg/Vs)+(Greaterof(-*(F+indexP*6+2),0.0)));
			
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+3))*Volg/Vn)+(Greaterof(-*(F+indexP*6+3),0.0)));
			
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+4))*Volg/Vb)+(Greaterof(-*(F+indexP*6+4),0.0)));
		
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(((*(F+indexP*6+5))*Volg/Vt)+(Greaterof(-*(F+indexP*6+5),0.0)));
//------------------------------------------- End of Boundary cells ----------------------------------------
	
//------------------------------------------- Interior cells ----------------------------------------------(27)---------
//-------------------------------------------------------------------------------------------------------Q1-----------------
		k = 1;j = 1;i = 1;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//------------------------------------------------------------------------------------------------------Q2----------------
		j = 1;i = 1;
		for(k=2;k<NK-3;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------Q3-------------------
		k = NK-3;j = 1;i = 1;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//--------------------------------------------------------------------------------------------------Q4--------------------
		k = 1;j = NJ-3;i = 1;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//-----------------------------------------------------------------------------------------------------Q5-----------------
		j = NJ-3;i = 1;
		for(k=2;k<NK-3;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//--------------------------------------------------------------------------------------------------Q6--------------------
		k = NK-3;j = NJ-3;i = 1;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//--------------------------------------------------------------------------------------------------Q7--------------------
		k = 1;j = 1;i = NI-3;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//------------------------------------------------------------------------------------------------------Q8----------------
		j = 1;i = NI-3;
		for(k=2;k<NK-3;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------Q9-------------------
		k = NK-3;j = 1;i = NI-3;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//--------------------------------------------------------------------------------------------------Q10--------------------
		k = 1;j = NJ-3;i = NI-3;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//-----------------------------------------------------------------------------------------------------Q11-----------------
		j = NJ-3;i = NI-3;
		for(k=2;k<NK-3;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//--------------------------------------------------------------------------------------------------Q12--------------------
		k = NK-3;j = NJ-3;i = NI-3;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		if(vartag==6)
			convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
		else
			convectivecoeff = 1.0;
		
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vee = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vnn = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vtt = *(slabVolume[4]+(NI-1)*j+i);
		
		Ve = Volg + Vol[1];
		Vw = Volg + Vol[0];
		Vn = Volg + Vol[3];
		Vs = Volg + Vol[2];
		Vt = Volg + Vol[5];
		Vb = Volg + Vol[4];
		
		*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
								+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
								+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
								+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
								+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
								+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
		//west
		*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
								-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								+(Greaterof(-*(F+indexP*6+0),0.0)));
		//west to west
		*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
		//east
		*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
								-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
								-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
								+(Greaterof(-*(F+indexP*6+1),0.0)));
		//east to east
		*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
		//south
		*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
								-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								+(Greaterof(-*(F+indexP*6+2),0.0)));
		//south to south
		*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
		//north
		*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
								-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
								-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
								+(Greaterof(-*(F+indexP*6+3),0.0)));
		//north to north
		*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
		//bottom
		*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
								-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								+(Greaterof(-*(F+indexP*6+4),0.0)));
		//bottom to bottom
		*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
		//top
		*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
								-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
								-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
								+(Greaterof(-*(F+indexP*6+5),0.0)));
		//top to top
		*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
//--------------------------------------------------------------------------------------------------------Q13--------------
		i = 1;k = 1;
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------Q14--------------------
		i = 1;k = NK-3;
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------Q15--------------------
		i = NI-3;k = 1;
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-----------------------------------------------------------------------------------------------------------Q16------------
		i = NI-3;k = NK-3;
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//------------------------------------------------------------------------------------------------------Q17-----------------
		i = 1;
		for(k=2;k<NK-3;k++)
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-------------------------------------------------------------------------------------------------------Q18----------------
		i = NI-3;
		for(k=2;k<NK-3;k++)
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//----------------------------------------------------------------------------------------------------Q19-------------------
		k = 1;
		for(i=2;i<NI-3;i++)
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//------------------------------------------------------------------------------------------------------Q20-----------------
		k = NK-3;
		for(i=2;i<NI-3;i++)
		for(j=2;j<NJ-3;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------------------------------------------------------------------Q21--------------------
		j = 1;
		for(k=2;k<NK-3;k++)
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-----------------------------------------------------------------------------------------------------Q22------------------
		j = NJ-3;
		for(k=2;k<NK-3;k++)
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//--------------------------------------------------------------------------------------------------Q23---------------------
		k = 1;
		j = 1;
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-----------------------------------------------------------------------------------------------------------Q24----------
		k = NK-3;
		j = 1;
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//------------------------------------------------------------------------------------------------------------Q25----------
		k = 1;
		j = NJ-3;
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//-----------------------------------------------------------------------------------------------------------Q26----------
		k = NK-3;
		j = NJ-3;
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(slabVolume[4]+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//----------------------------------------------------------------------------------------------------Q27---------------
		for(k=2;k<NK-3;k++)
		for(j=2;j<NJ-3;j++)
		for(i=2;i<NI-3;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else
				convectivecoeff = 1.0;
			
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vww = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-2);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vee = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+2);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vss = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-2)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vnn = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+2)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vbb = *(volumeDomain+(NI-1)*(NJ-1)*(k-2)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vtt = *(volumeDomain+(NI-1)*(NJ-1)*(k+2)+(NI-1)*j+i);
			
			Ve = Volg + Vol[1];
			Vw = Volg + Vol[0];
			Vn = Volg + Vol[3];
			Vs = Volg + Vol[2];
			Vt = Volg + Vol[5];
			Vb = Volg + Vol[4];
			
			*(diffdia + indexP) = 0.5*delta_t*convectivecoeff*alpha*(((Vol[1]*(2*Volg+Vol[0])/(Ve*Vw))*(Greaterof(*(F+indexP*6+1),0.0))-(Vol[1]*(2*Vol[1]+Vee)/(Ve*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))-(Greaterof(*(F+indexP*6+1),0.0)))
									+((Vol[0]*(2*Volg+Vol[1])/(Ve*Vw))*(Greaterof(*(F+indexP*6+0),0.0))-(Vol[0]*(2*Vol[0]+Vww)/(Vw*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))-(Greaterof(*(F+indexP*6+0),0.0)))
									+((Vol[3]*(2*Volg+Vol[2])/(Vn*Vs))*(Greaterof(*(F+indexP*6+3),0.0))-(Vol[3]*(2*Vol[3]+Vnn)/(Vn*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))-(Greaterof(*(F+indexP*6+3),0.0)))
									+((Vol[2]*(2*Volg+Vol[3])/(Vn*Vs))*(Greaterof(*(F+indexP*6+2),0.0))-(Vol[2]*(2*Vol[2]+Vss)/(Vs*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))-(Greaterof(*(F+indexP*6+2),0.0)))
									+((Vol[5]*(2*Volg+Vol[4])/(Vt*Vb))*(Greaterof(*(F+indexP*6+5),0.0))-(Vol[5]*(2*Vol[5]+Vtt)/(Vt*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))-(Greaterof(*(F+indexP*6+5),0.0)))
									+((Vol[4]*(2*Volg+Vol[5])/(Vt*Vb))*(Greaterof(*(F+indexP*6+4),0.0))-(Vol[4]*(2*Vol[4]+Vbb)/(Vb*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))-(Greaterof(*(F+indexP*6+4),0.0))));
			//west
			*(CoeffN[0]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[1])*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									-((2*Vol[0]+Vww)*Volg/(Vw*(Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0))
									-(Vol[1]*Volg/(Vw*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									+(Greaterof(-*(F+indexP*6+0),0.0)));
			//west to west
			*(CoeffN[6]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[0]*Volg/((Vol[0]+Vww)*(Volg+2*Vol[0]+Vww)))*(Greaterof(-*(F+indexP*6+0),0.0)));
			//east
			*(CoeffN[1]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[0])*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+1),0.0))
									-((2*Vol[1]+Vee)*Volg/(Ve*(Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0))
									-(Vol[0]*Volg/(Ve*(Vol[1]+2*Volg+Vol[0])))*(Greaterof(*(F+indexP*6+0),0.0))
									+(Greaterof(-*(F+indexP*6+1),0.0)));
			//east to east
			*(CoeffN[7]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[1]*Volg/((Vol[1]+Vee)*(Volg+2*Vol[1]+Vee)))*(Greaterof(-*(F+indexP*6+1),0.0)));
			//south
			*(CoeffN[2]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[3])*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									-((2*Vol[2]+Vss)*Volg/(Vs*(Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0))
									-(Vol[3]*Volg/(Vs*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									+(Greaterof(-*(F+indexP*6+2),0.0)));
			//south to south
			*(CoeffN[8]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[2]*Volg/((Vol[2]+Vss)*(Volg+2*Vol[2]+Vss)))*(Greaterof(-*(F+indexP*6+2),0.0)));
			//north
			*(CoeffN[3]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[2])*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+3),0.0))
									-((2*Vol[3]+Vnn)*Volg/(Vn*(Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0))
									-(Vol[2]*Volg/(Vn*(Vol[3]+2*Volg+Vol[2])))*(Greaterof(*(F+indexP*6+2),0.0))
									+(Greaterof(-*(F+indexP*6+3),0.0)));
			//north to north
			*(CoeffN[9]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[3]*Volg/((Vol[3]+Vnn)*(Volg+2*Vol[3]+Vnn)))*(Greaterof(-*(F+indexP*6+3),0.0)));
			//bottom
			*(CoeffN[4]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[5])*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									-((2*Vol[4]+Vbb)*Volg/(Vb*(Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0))
									-(Vol[5]*Volg/(Vb*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									+(Greaterof(-*(F+indexP*6+4),0.0)));
			//bottom to bottom
			*(CoeffN[10]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[4]*Volg/((Vol[4]+Vbb)*(Volg+2*Vol[4]+Vbb)))*(Greaterof(-*(F+indexP*6+4),0.0)));
			//top
			*(CoeffN[5]+indexP) = 0.5*delta_t*convectivecoeff*alpha*(+((2*Volg+Vol[4])*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+5),0.0))
									-((2*Vol[5]+Vtt)*Volg/(Vt*(Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0))
									-(Vol[4]*Volg/(Vt*(Vol[5]+2*Volg+Vol[4])))*(Greaterof(*(F+indexP*6+4),0.0))
									+(Greaterof(-*(F+indexP*6+5),0.0)));
			//top to top
			*(CoeffN[11]+indexP) = 0.5*delta_t*convectivecoeff*alpha*((Vol[5]*Volg/((Vol[5]+Vtt)*(Volg+2*Vol[5]+Vtt)))*(Greaterof(-*(F+indexP*6+5),0.0)));
		}
//---------------------------------------End of Interior Cells------------------------------------------
	}
	
	return(0);
}


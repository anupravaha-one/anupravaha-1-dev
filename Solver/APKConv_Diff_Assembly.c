#include "APKConv_Diff_Assembly.h"


/**
This function is used for calculating convection-diffusion equation coefficients in case of BTCS time stepping scheme.
*/

double Conv_Diff_Assembly(Base *base,Zone *zone,int g,int vartag)
{
	int NI,NJ,NK,i,j,k,indexP,*ptrSize,flagrho,flagSpheat,flagCPDrop;//rajendra
	double *Ddiffdia,*DCoeffN[18],delta_t;
	double *volumeDomain;
	double rho,convectivecoeff;
	double Temp,*CoeffN[18],*diffdia;
	double Volg;
	double *conv,*conv_nb[6],*RHO,*CP,*CPDrop;//shash;
	
	int Timetag;
	int VarFlag=0;//This keeps track of which variable it is, and later on only by checking the value of this one variable suitable action can be taken. This is done to avoid lengthy "if" statements where a lot of conditions are being tested
	
	delta_t = base->BaseParameters->TimeStep;
	
	if(vartag==600||vartag==601)
		delta_t = base->BaseParameters->Turb->DualTimeStep;//delta_t/10.0;
	
	Timetag = base->BaseParameters->TimeStepingSchemeTag;
	
	RHO=zone->Flowsol->Rho;
	CP =(zone->Flowsol->SpecificHeat);
	
	CPDrop =(zone->Flowsol->CPdrop);//shash
	
	diffdia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
	}
	
	if((vartag>=41)&&(vartag<=55))
	{
		Ddiffdia = zone->DeferredCorrCoeffDummy_P;
		for(i=0;i<18;i++)
		{
			DCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
		}
	}

	
	//disperse phase....adhi+shash
	else if((vartag==62)||(vartag==63)||(vartag==64)||(vartag==65)||(vartag==710)||(vartag==69))//||(vartag==67)
	{
		if(base->flagTurb==0)//chandrasekhar--as in case of turbulence diffusion coefficients are there.
		{
			Ddiffdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<18;i++)
			{
			DCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
		}
		else if(base->flagTurb==1)
		{
			Ddiffdia = (zone->grid+g)->DiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
			DCoeffN[i] = (zone->grid+g)->DiffusionCoeffPhi_nb[i];
			}
		}
		VarFlag=1;
	}
	else if(vartag==35)
	{
		Ddiffdia = (zone->grid+g)->DiffusionCoeffPhi_P;
		for(i=0;i<18;i++)
		{
			DCoeffN[i] = (zone->grid+g)->DiffusionCoeffPhi_nb[i];
		}
		VarFlag=1;
	}
	else
	{
		Ddiffdia = (zone->grid+g)->DiffusionCoeffPhi_P;
		for(i=0;i<18;i++)
		{
			DCoeffN[i] = (zone->grid+g)->DiffusionCoeffPhi_nb[i];
		}
	}
	
	conv = (zone->grid+g)->Conv_Coeff_P;
	for(i=0;i<6;i++)
	{
		conv_nb[i] = (zone->grid+g)->Conv_Coeff_nb[i];
	}
	
	if(vartag==9||vartag==56||vartag==57||vartag==58)//raj
	{
		VarFlag=1;
	}
	if(base->flagMPFV==1)
	{
		VarFlag=1;
	}
	
	ptrSize = zone->Size ;
	volumeDomain = zone->ZGeom->Volume;
	
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	flagCPDrop = zone->Flowsol->flagCPdrop;
	flagrho = zone->Flowsol->flagRho;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	if(base->BaseParameters->flag2d == 1)
	{
//--------------------------------------------- Boundary cells ----------------------------------(1)---------
		indexP=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			if(VarFlag==1)//rajendra//adhi+shash//sushil
				rho = 1.0;
			else
				rho = *(RHO+indexP*flagrho);
			
			Temp=(delta_t/rho);
			
			if(vartag==6)
			        convectivecoeff = *(CP+flagSpheat*indexP);   //setting this =0 gave problem for 2D case - try 2D lid-driven cavity with energy eqn
			else if(vartag==69)
				convectivecoeff = *(CPDrop+flagCPDrop*indexP);
			else
				convectivecoeff = 1.0;
			
			
			Volg = *(volumeDomain+indexP);
			
			if(Timetag==1)
			{
				convectivecoeff = convectivecoeff/2.0;
				Temp = Temp/2.0;
				*(diffdia + indexP) = convectivecoeff*Volg*2.0 + convectivecoeff*(*(conv+indexP)) - Temp*(*(Ddiffdia + indexP));
			}
			else
			{
				*(diffdia + indexP) = convectivecoeff*Volg + convectivecoeff*(*(conv+indexP)) - Temp*(*(Ddiffdia + indexP));
				
				
			}
			
			*(CoeffN[0]+indexP) = -convectivecoeff*(*(conv_nb[0]+indexP))- Temp*(*(DCoeffN[0]+indexP));
			
			*(CoeffN[1]+indexP) = -convectivecoeff*(*(conv_nb[1]+indexP))- Temp*(*(DCoeffN[1]+indexP));
			
			*(CoeffN[2]+indexP) = -convectivecoeff*(*(conv_nb[2]+indexP))- Temp*(*(DCoeffN[2]+indexP));
			
			*(CoeffN[3]+indexP) = -convectivecoeff*(*(conv_nb[3]+indexP))- Temp*(*(DCoeffN[3]+indexP));
			
			*(CoeffN[6]+indexP) = -Temp*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -Temp*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[10]+indexP) = -Temp*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -Temp*(*(DCoeffN[11]+indexP));
			indexP++;
		}
	}
	else if(base->BaseParameters->flag2d == 0)
	{
		indexP=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			if(VarFlag==1)//rajendra
				rho = 1.0;
			else
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			Temp=(delta_t/rho);
			
			if(vartag==6)
				convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
			else if(vartag==69)
				convectivecoeff = *(CPDrop+flagCPDrop*indexP);
			else
				convectivecoeff = 1.0;
				
			Volg = *(volumeDomain+indexP);
			if(Timetag==1)
			{
				convectivecoeff = convectivecoeff/2.0;
				Temp = Temp/2.0;
				*(diffdia + indexP) = convectivecoeff*Volg*2.0 + convectivecoeff*(*(conv+indexP)) - Temp*(*(Ddiffdia + indexP));
			}
			else
			{
				*(diffdia + indexP) = convectivecoeff*Volg + convectivecoeff*(*(conv+indexP)) - Temp*(*(Ddiffdia + indexP));
			}
			
			*(CoeffN[0]+indexP) = -convectivecoeff*(*(conv_nb[0]+indexP))- Temp*(*(DCoeffN[0]+indexP));
				
			*(CoeffN[1]+indexP) = -convectivecoeff*(*(conv_nb[1]+indexP))- Temp*(*(DCoeffN[1]+indexP));
			
			*(CoeffN[2]+indexP) = -convectivecoeff*(*(conv_nb[2]+indexP))- Temp*(*(DCoeffN[2]+indexP));
			
			*(CoeffN[3]+indexP) = -convectivecoeff*(*(conv_nb[3]+indexP))- Temp*(*(DCoeffN[3]+indexP));
			
			*(CoeffN[4]+indexP) = -convectivecoeff*(*(conv_nb[4]+indexP))- Temp*(*(DCoeffN[4]+indexP));
			
			*(CoeffN[5]+indexP) = -convectivecoeff*(*(conv_nb[5]+indexP))- Temp*(*(DCoeffN[5]+indexP));
			
			*(CoeffN[6]+indexP) = -Temp*(*(DCoeffN[6]+indexP));
			
			*(CoeffN[7]+indexP) = -Temp*(*(DCoeffN[7]+indexP));
			
			*(CoeffN[8]+indexP) = -Temp*(*(DCoeffN[8]+indexP));
			
			*(CoeffN[9]+indexP) = -Temp*(*(DCoeffN[9]+indexP));
			
			*(CoeffN[10]+indexP) = -Temp*(*(DCoeffN[10]+indexP));
			
			*(CoeffN[11]+indexP) = -Temp*(*(DCoeffN[11]+indexP));
			
			*(CoeffN[12]+indexP) = -Temp*(*(DCoeffN[12]+indexP));
			
			*(CoeffN[13]+indexP) = -Temp*(*(DCoeffN[13]+indexP));
			
			*(CoeffN[14]+indexP) = -Temp*(*(DCoeffN[14]+indexP));
			
			*(CoeffN[15]+indexP) = -Temp*(*(DCoeffN[15]+indexP));
			
			*(CoeffN[16]+indexP) = -Temp*(*(DCoeffN[16]+indexP));
			
			*(CoeffN[17]+indexP) = -Temp*(*(DCoeffN[17]+indexP));
			
			indexP++;
		}
	}
	else
	{
		printf("Wrong tag\nExiting in APKConv_Diff_Assembly.c\n");
		exit(12);
	}

	return(0);
}





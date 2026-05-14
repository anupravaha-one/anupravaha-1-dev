#include "APKBC_userDefinedVarCoefficients.h"

/**
In this function is used to specify the boundary face center value of density as a function of flow variables.
User can write the boundary face center value of density in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncDensity(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	double f,Rho_s,Rho_f,E_fm,E_sm;
        double P,Temp,TotMolWt; //Shyam
	int nop,nophase=0,flagrho=0.0;
// 	double omegaA,omegaB,MolWtA,MolWtB;
// 	double xA,xB;
	
	nop=base->nphase;
		
	//printf("updating density\n");
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		//Do not Disturb this if condition.
		Rho_s = base->BaseParameters->Rho_s;
		Rho_f = base->BaseParameters->Rho_f;
		E_sm = base->BaseParameters->E_sm;
		E_fm = base->BaseParameters->E_fm;
		
		if(base->flagSolidification==1)
		{
			*((bcPtr->ZSlab)->f+faceCenterindex) = *((bcPtr->ZSlab)->Latent_Heat_Content+faceCenterindex)/(E_fm - E_sm);
			f = *((bcPtr->ZSlab)->f+faceCenterindex);
			//This statement has been written keeping in mind that variation of variable
			//is due to phase change only.  Others conditions can also be written in
			//this block if needed.
			variable = Rho_s + (Rho_f - Rho_s) *f;
			//printf("variable = %lf\n",variable);
			//*((z->Flowsol)->Rho+i) = Rho_s + (Rho_f - Rho_s) *f;	
		}
		else
		{
			*((bcPtr->ZSlab)->f+faceCenterindex) = (*((bcPtr->ZSlab)->enthalpy+faceCenterindex)- E_fm)/(E_sm - E_fm);
			f = *((bcPtr->ZSlab)->f+faceCenterindex);
			variable =  Rho_f /(1.0 - f + (Rho_f/Rho_s)*f);
		}
		
	}
	else if(base->flagTurbComb==1&&base->flagVariableDensityProblem==1) //Shyam
	{
		P = *((bcPtr->ZSlab)->P+faceCenterindex);
		Temp = *((bcPtr->ZSlab)->T+faceCenterindex);
		TotMolWt = *((bcPtr->ZSlab)->MolWteq+faceCenterindex);
		variable = ((P+101325.0)*TotMolWt)/(8314.0*Temp);
	}
	else if(base->flagMPFV==1)// Added by sushil...for multiphase flow....
	{	
		*((bcPtr->ZSlab)->Rho+faceCenterindex)=0.0;
// 		printf("\\n\nIf density of each phase is variable then user have to write separate function to calculate density of each phase\n\n");
		for(nophase=0;nophase<nop;nophase++)
		{
//  			flagrho=((bcPtr->ZSlab)->flagRhoPhase[nophase]);
			flagrho=base->BaseParameters->flagRhoPhase[nophase];
			
			if(flagrho==1)
			{
				*((bcPtr->ZSlab)->RhoofPhase[nophase]+faceCenterindex)=1.0;//Default value
				*((bcPtr->ZSlab)->Rho+faceCenterindex) =*((bcPtr->ZSlab)->Rho+faceCenterindex) + *((bcPtr->ZSlab)->RhoofPhase[nophase]+(faceCenterindex*flagrho))*(*((bcPtr->ZSlab)->VarthetaPhase[nophase]+faceCenterindex));
				
			}
			
			else
			{	
				*((bcPtr->ZSlab)->Rho+faceCenterindex) =*((bcPtr->ZSlab)->Rho+faceCenterindex) + *((bcPtr->ZSlab)->RhoofPhase[nophase]+faceCenterindex*flagrho)*(*((bcPtr->ZSlab)->VarthetaPhase[nophase]+faceCenterindex));
			}
		}	
	
		variable=*((bcPtr->ZSlab)->Rho+faceCenterindex);
	}
	else if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
	{
		variable=*((bcPtr->ZSlab)->Rho+faceCenterindex);
		//we are solving equation for the density of carrier phase, which is done just before disperse phase equation.
		//Updation of boundary is done as a solvable variable.
	}
	else
	{
		//Write User defined functions here >>>>>
		
		/*double omegaA,omegaB,MolWtA,MolWtB,TotMolWt;
		double xA,xB;
		
		omegaA = (*(bcPtr->ZSlab->Scal[0]+faceCenterindex));
		omegaB = (*(bcPtr->ZSlab->Scal[1]+faceCenterindex));
		//MolWtA = 39.944;
		//MolWtB = 2.016;
		MolWtA = 34.29548;
		MolWtB = 25.536;
		
		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
		
		TotMolWt = MolWtA*xA + MolWtB*xB;
		variable = (101325.0 * TotMolWt)/(8314.0*273.0);
		//variable = (101325.0 * TotMolWt)/(8314.0*(*(bcPtr->ZSlab->T+faceCenterindex)));
		//printf("slab value of temp = %lf Density = %lf\n",*(bcPtr->ZSlab->T+faceCenterindex),variable);*/
		//variable = (100000.0 * 0.001)/(*(bcPtr->ZSlab->T+faceCenterindex));
		//variable = 1.0;
		//variable = (100000.0* (28.9/8314.0))/(*(bcPtr->ZSlab->T+faceCenterindex));
		
		//For Saclay Test Problem
// 		omegaA = (*(bcPtr->ZSlab->Scal[2]+faceCenterindex));
// 		omegaB = (*(bcPtr->ZSlab->Scal[3]+faceCenterindex));
// 		
// 		MolWtA = 28.84;
// 		MolWtB = 4.0;
// 	
// 		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 	
// 		TotMolWt = MolWtA*xA + MolWtB*xB;
// 		variable = (101325.0 * TotMolWt)/(8314.0*298.0);
		
	}
	
	return(variable);
}
/**
In this function is used to specify the boundary face center value of viscosity as a function of flow variables.
User can write the boundary face center value of viscosity in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncViscosity(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable = 0.01;	
// 	double omegaA,omegaB,MolWtA,MolWtB;
// 	double xA,xB,muA,muB;
// 	double phiAA,phiAB,phiBB,phiBA;
	double Mu_s,Mu_f,f;
	int nop,nophase=0,flagvisc=0.0;
	
	nop=base->nphase;
	
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		//Do not Disturb this if condition.
		Mu_s = base->BaseParameters->Mu_s;
		Mu_f = base->BaseParameters->Mu_f;
		f = *((bcPtr->ZSlab)->f+faceCenterindex);
		//This statement has been written keeping in mind that variation of variable
		//is due to phase change only.  Others conditions can also be written in
		//this block if needed.
		if(base->flagSolidification==1)
		variable = Mu_s + (Mu_f - Mu_s) *f;
		else
		variable = f*Mu_s + (1.0-f)*Mu_f;
	}
// 	printf("\n\nUser have to write a separate function to evaluate viscosity of each phase\n\n");
	else if(base->flagMPFV==1)//Added by sushil for multiphase flow problems.....
	{
		*((bcPtr->ZSlab)->Viscoeff+faceCenterindex)=0.0;
// 		printf("\\n\nIf viscosity of each phase is variable then user have to write separate function to calculate viscosity\n\n");
		for(nophase=0;nophase<nop;nophase++)
		{
			flagvisc=base->BaseParameters->flagViscoeffPhase[nophase];
			if(flagvisc==1)
			{
				*((bcPtr->ZSlab)->ViscosityofPhase[nophase]+faceCenterindex)=1.0;//Default value
				
				*((bcPtr->ZSlab)->ViscoeffDummy[nophase]+faceCenterindex) =(*((bcPtr->ZSlab)->ViscosityofPhase[nophase]+faceCenterindex*flagvisc)*(*((bcPtr->ZSlab)->VarthetaPhase[nophase]+faceCenterindex)));
				
				
				*((bcPtr->ZSlab)->Viscoeff+faceCenterindex) =*((bcPtr->ZSlab)->Viscoeff+faceCenterindex)+ *((bcPtr->ZSlab)->ViscosityofPhase[nophase]+faceCenterindex*flagvisc)*(*((bcPtr->ZSlab)->VarthetaPhase[nophase]+faceCenterindex));
				
			}
			else
			{
				
				*((bcPtr->ZSlab)->ViscoeffDummy[nophase]+faceCenterindex) =(*((bcPtr->ZSlab)->ViscosityofPhase[nophase]+faceCenterindex*flagvisc)*(*((bcPtr->ZSlab)->VarthetaPhase[nophase]+faceCenterindex)));
				
				*((bcPtr->ZSlab)->Viscoeff+faceCenterindex) =*((bcPtr->ZSlab)->Viscoeff+faceCenterindex)+ *((bcPtr->ZSlab)->ViscosityofPhase[nophase]+faceCenterindex*flagvisc)*(*((bcPtr->ZSlab)->VarthetaPhase[nophase]+faceCenterindex));
			}
			//printf("The updated dummy viscosity of the phase %d at %s is %lf\n",nophase,bcPtr->BCName,*((bcPtr->ZSlab)->ViscoeffDummy[nophase]+faceCenterindex));
		}
		
	
		variable=*((bcPtr->ZSlab)->Viscoeff+faceCenterindex);
		//printf("\n index=%d\tvisco=%e",faceCenterindex,variable);
	}
	else
	{
		//Write User defined functions here >>>>>
		
// 		omegaA = (*(bcPtr->ZSlab->Scal[0]+faceCenterindex));
// 		omegaB = (*(bcPtr->ZSlab->Scal[1]+faceCenterindex));
// 		//MolWtA = 39.94;
// 		//MolWtB = 2.016;
// 		MolWtA = 34.29548;
// 		MolWtB = 25.536;
// 		muA = 0.01255;
// 		muB= 0.01;
// 
// 		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 
// 		variable = (xA*muA*sqrt(MolWtA) + xB*muB*sqrt(MolWtB))/(xA*sqrt(MolWtA) + xB*sqrt(MolWtB));
		
// 		phiAA = 1.0;
// 		phiBB = 1.0;
// 		phiAB = (1.0/sqrt(8.0))* pow((1.0+(MolWtA/MolWtB)),-0.5) * pow((1.0+pow((muA/muB),0.5)*pow((MolWtB/MolWtA),0.25)),2.0);
// 		phiBA = (1.0/sqrt(8.0))* pow((1.0+(MolWtB/MolWtA)),-0.5) * pow((1.0+pow((muB/muA),0.5)*pow((MolWtA/MolWtB),0.25)),2.0);
// 	
// 		variable = (xA*muA/(xA*phiAA+xA*phiAB))+(xB*muB/(xB*phiBA+xB*phiBB));
	
		//For Saclay Test Problem
// 		omegaA = (*(bcPtr->ZSlab->Scal[2]+faceCenterindex));
// 		omegaB = (*(bcPtr->ZSlab->Scal[3]+faceCenterindex));
// 		
// 		MolWtA = 28.84;
// 		MolWtB = 4.0;
// 		muA = 0.00001855;
// 		muB=  0.00001989;
// 
// 		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 
// 		variable = (xA*muA*sqrt(MolWtA) + xB*muB*sqrt(MolWtB))/(xA*sqrt(MolWtA) + xB*sqrt(MolWtB));
	}
	
	return(variable);
}



/**
In this function is used to specify the boundary face center value of ThermalConductivity as a function of flow variables.
User can write the boundary face center value of ThermalConductivity in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncThermalConductivity(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable = 1.0;
	//variable = 2.40 *X;
	//double omegaA,omegaB,MolWtA,MolWtB,variable;
	//double xA,xB,muA,muB,kA,kB;
	//double phiAA,phiAB,phiBB,phiBA;
	double K_s,K_f,f,CP_s,CP_f;
	//for evaporation of droplet phase
	double Ka,Kv,massfrac;//shash
	
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		K_s = base->BaseParameters->K_s;
		K_f = base->BaseParameters->K_f;
		CP_s = base->BaseParameters->CP_s;
		CP_f = base->BaseParameters->CP_f;
		f = *((bcPtr->ZSlab)->f+faceCenterindex);
		//This statement has been written keeping in mind that variation of variable
		//is due to phase change only.  Others conditions can also be written in
		//this block if needed.
		if(base->flagSolidification==1)
		variable = (K_s + (K_f - K_s) *f) / (CP_s + (CP_f - CP_s) *f);
		else
		variable = (f *K_s + (1.0-f) *K_f)/(f * CP_s + (1.0-f)*CP_f);
	}
	else if(base->flagDroplet==1)
	{
	 	Ka = 0.0242;
	 	Kv = 0.0128;
	 	massfrac = 0.0;
		//if mass fraction is not being solved, the mixture properties will become the properties of pure gas phase.
		if(base->BaseParameters->flagDropFuelMassFrac==1)
	 	{
	 	massfrac = *(bcPtr->ZSlab->DropFuelMassFrac+faceCenterindex);
	 	variable = massfrac*Kv+(1.0-massfrac)*Ka;
		}
		else
		{
		variable = massfrac*Kv+(1.0-massfrac)*Ka;
		}
	}
	else
	{
		/*omegaA = (*(bcPtr->ZSlab->Scal[0]+faceCenterindex));
		omegaB = (*(bcPtr->ZSlab->Scal[1]+faceCenterindex));
		MolWtA = 39.94;
		MolWtB = 2.016;
		muA= 0.01;
		muB= 0.01;
		kA = 0.01;
		kB= 0.01;

		phiAA = 1.0;
		phiBB = 1.0;
		phiAB = (1.0/sqrt(8.0))* pow((1.0+(MolWtA/MolWtB)),-0.5) * pow((1.0+pow((muA/muB),0.5)*pow((MolWtB/MolWtA),0.25)),2.0);
		phiBA = (1.0/sqrt(8.0))* pow((1.0+(MolWtB/MolWtA)),-0.5) * pow((1.0+pow((muB/muA),0.5)*pow((MolWtA/MolWtB),0.25)),2.0);

		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));

		variable = (xA*kA/(xA*phiAA+xA*phiAB))+(xB*kB/(xB*phiBA+xB*phiBB));*/
	}
	return(variable);
}
/**
In this function is used to specify the boundary face center value of SpecificHeat as a function of flow variables.
User can write the boundary face center value of SpecificHeat in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncSpecificHeat(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable = 1.0;
	double CP_s,CP_f,f;
	//for evaporation droplet phase flow
	double Cpa,Cpv,massfrac;
	
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		CP_s = base->BaseParameters->CP_s;
		CP_f = base->BaseParameters->CP_f;
		f = *((bcPtr->ZSlab)->f+faceCenterindex);
		//This statement has been written keeping in mind that variation of variable
		//is due to phase change only.  Others conditions can also be written in
		//this block if needed.
		if(base->flagSolidification==1)
		variable = CP_s + (CP_f - CP_s) *f;
		else
		variable = f * CP_s + (1.0-f)*CP_f ;
	}
	else if(base->flagDroplet==1)//shash
	{
	 	Cpa = 1006.43;
	 	Cpv = 2471.0;
	 	massfrac = 0.0;
		//if mass fraction is not being solved, the mixture properties will become the properties of pure gas phase.
		if(base->BaseParameters->flagDropFuelMassFrac==1)
	 	{
	 	massfrac = *(bcPtr->ZSlab->DropFuelMassFrac+faceCenterindex);
	 	variable = massfrac*Cpv+(1.0-massfrac)*Cpa;
		}
		else
		{
		variable = massfrac*Cpv+(1.0-massfrac)*Cpa;
		}
	}
	else
	{
		/*double omegaA,omegaB,MolWtA,MolWtB,variable;
		double xA,xB,CpA,CpB,muA,muB;

		omegaA = (*(bcPtr->ZSlab->Scal[0]+faceCenterindex));
		omegaB = (*(bcPtr->ZSlab->Scal[1]+faceCenterindex));
		MolWtA = 39.94;
		MolWtB = 2.016;
		CpA = 1.0;
		CpB = 1.0;
		muA= 0.01;
		muB= 0.01;

		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));

		variable = (xA*CpA + xB*CpB);*/
	}
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff1 as a function of flow variables.
User can write the boundary face center value of ScalCoeff1 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff1(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 1)
	{
		printf("Exiting from BCuserfuncScalCoeff1 \n");
		exit(88);
	}
	variable = 0.0;
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff2 as a function of flow variables.
User can write the boundary face center value of ScalCoeff2 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff2(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 2)
	{
		printf("Exiting from BCuserfuncScalCoeff2 \n");
		exit(88);
	}
	variable = 0.0;
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff3 as a function of flow variables.
User can write the boundary face center value of ScalCoeff3 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff3(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 3)
	{
		printf("Exiting from BCuserfuncScalCoeff3 \n");
		exit(88);
	}
	//For Species/Sclars in Turbulent flow used for Saclay Test Problem.
// 	variable = 7.0e-05 + *(bcPtr->ZSlab->TViscoeff+faceCenterindex)/0.7;
	variable = 0.0;
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff4 as a function of flow variables.
User can write the boundary face center value of ScalCoeff4 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff4(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 4)
	{
		printf("Exiting from BCuserfuncScalCoeff4 \n");
		exit(88);
	}
	//For Species/Sclars in Turbulent flow used for Saclay Test Problem.
// 	variable = 7.0e-05 + *(bcPtr->ZSlab->TViscoeff+faceCenterindex)/0.7;
	variable = 0.0;
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff5 as a function of flow variables.
User can write the boundary face center value of ScalCoeff5 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff5(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//double Tx,Ty,Tz;
	/*Tx=(X+Y+Z)/sqrt(3.0);
	Ty=(Y-Z)/sqrt(2.0);
	Tz=(2.0*X-Y-Z)/sqrt(6.0);
	variable=2*Ty+3*Tz+4;*/
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff6 as a function of flow variables.
User can write the boundary face center value of ScalCoeff6 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff6(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff7 as a function of flow variables.
User can write the boundary face center value of ScalCoeff7 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff7(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff8 as a function of flow variables.
User can write the boundary face center value of ScalCoeff8 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff8(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff9 as a function of flow variables.
User can write the boundary face center value of ScalCoeff9 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff9(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function is used to specify the boundary face center value of ScalCoeff10 as a function of flow variables.
User can write the boundary face center value of ScalCoeff10 in the form of functions. Face center coordinates X,Y,Z,
zone boundary pointer and face center index is given as input.
*/
double BCuserfuncScalCoeff10(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}

double BCuserfuncChargeDiffusion(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=5.32e-06;
	return(variable);
}
/**
In this function is used to specify the boundary face center value of SpecificHeat of liquid droplet as a function of flow variables like droplet temperature. User can write the boundary face center value of SpecificHeat in the form of functions. Face center coordinates X, Y , Z, zone boundary pointer and face center index is given as input.
*/
double BCuserfuncCPdrop(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base)
{
	double variable = 1.0;
	
	double Cpf,Tsurf;
	
	if(base->flagDroplet==1)//shash
	{
	 Tsurf = *((bcPtr->ZSlab)->Tdrop+faceCenterindex);
	 //liquid droplet specific heat for n-Heptane.
	 //Cpf=13058.45044066-126.5095282565*Tsurf+0.5279613848638*Tsurf*Tsurf-0.0009457386295687*Tsurf*Tsurf*Tsurf+6.369853422618e-7*Tsurf*Tsurf*Tsurf*Tsurf;
	 Cpf = 2219.0;
	 variable = Cpf;
	}
	else
	{
	variable = 1.0;
	}
	return(variable);
}

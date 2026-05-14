#include "APKuserDefinedVarCoefficients.h"

double userEnterconstant(Zone *z,int cellindex,Base *base)
{
	double constant;
	double *MolWt;
	//double omegaA,omegaB,MolWtA,MolWtB,TotMolWt;
	//double xA,xB;

	if(base->flagTwoPhase==1)//calculation of del(Rho)/del(P) for the water-steam mixture
	{
		double f,constant_L,constant_V;// constant_L and constant_V are del(Rho)/del(P) values for water and steam from thermodynamic propety tables
		constant_L = base->BaseParameters->Pc_l;
		constant_V = base->BaseParameters->Pc_v;
		f=*(z->Flowsol->f+cellindex);
		constant=(f*constant_V)+((1-f)*constant_L);
	}
	else
	{
		if(z->Flowsol->flagTurbComb==1)
		{
			MolWt = z->Flowsol->Tcoeff->MolWteq;

			constant = (*(MolWt + cellindex))/8314.0;
		}
		else
		{
			//Write User defined functions here >>>>>
			
// 			constant = 28.9/8314.0;
			
// 			omegaA = *(z->Flowsol->Scal[0]+cellindex);
// 			omegaB = *(z->Flowsol->Scal[1]+cellindex);
// 			//MolWtA = 39.944;
// 			//MolWtB = 2.016;
// 			MolWtA = 34.29548;
// 			MolWtB = 25.536;
// 		
// 			xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 			xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 		
// 			TotMolWt = MolWtA*xA + MolWtB*xB;
// 		
// 			constant = (TotMolWt/8314.0);//0.001;
			
			//For Saclay Test Problem
// 			omegaA = *(z->Flowsol->Scal[2]+cellindex);
// 			omegaB = *(z->Flowsol->Scal[3]+cellindex);
// 			MolWtA = 28.84;
// 			MolWtB = 4.0;
// 	
// 			xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 			xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 	
// 			TotMolWt = MolWtA*xA + MolWtB*xB;
// 	
// 			constant = (TotMolWt/8314.0);
		}

	}
	return(constant);
}
/**
In this function user can write the Density in the form of functions. Coordinate X,Y,Z, Zone pointer and
cell index is given as input. This function can be used if density is a function of flow variables like scalar concentrations etc.
*/
double userfuncDensity(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable=1.0;
// 	double omegaA,omegaB,MolWtA,MolWtB;
// 	double xA,xB;
        double Temp,P,TotMolWt;
	double f,Rho_s,Rho_f,E_fm,E_sm;
	
	int nop,nophase=0,flagrho=0.0;
	
	nop=base->nphase;
	
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		//Do not Disturb this if condition.
		Rho_s = base->BaseParameters->Rho_s;
		Rho_f = base->BaseParameters->Rho_f;
		E_sm = base->BaseParameters->E_sm;
		E_fm = base->BaseParameters->E_fm;
		
		if(base->flagSolidification==1)
		{
			*((z->Flowsol)->f+cellindex) = *((z->Flowsol)->Latent_Heat_Content+cellindex)/(E_fm - E_sm);
			f = *((z->Flowsol)->f+cellindex);
			//This statement has been written keeping in mind that variation of variable
			//is due to phase change only.  Others conditions can also be written in
			//this block if needed.
			variable = Rho_s + (Rho_f - Rho_s) *f;
			//*((z->Flowsol)->Rho+i) = Rho_s + (Rho_f - Rho_s) *f;
		}
		else
		{
			*((z->Flowsol)->f+cellindex) =(*((z->Flowsol)->enthalpy+cellindex)- E_fm)/(E_sm - E_fm);
			f = *((z->Flowsol)->f+cellindex);
			variable =  Rho_f /(1.0 - f + (Rho_f/Rho_s)*f);
		}
	}
	else if(base->flagTurbComb==1&&base->flagVariableDensityProblem == 1) //Shyam
	{
		Temp = *((z->Flowsol)->T+cellindex);
		P = *((z->Flowsol)->P+cellindex);
		TotMolWt = *(((z->Flowsol)->Tcoeff)->MolWteq + cellindex);
		variable = ((P+101325.0)*TotMolWt)/(8314.0*Temp);
	}
	else if(base->flagMPFV==1)// Added by sushil...for multiphase flow....
	{	
		*((z->Flowsol)->Rho+cellindex)=0.0;
		for(nophase=0;nophase<nop;nophase++)
		{
			flagrho=(base->BaseParameters->flagRhoPhase[nophase]);
			if(flagrho==1)
			{
				*((z->Flowsol)->RhoofPhase[nophase]+cellindex)=1.0;//Default value
				
				*((z->Flowsol)->Rho+cellindex) =*((z->Flowsol)->Rho+cellindex) + *((z->Flowsol)->RhoofPhase[nophase]+cellindex*flagrho)*(*((z->Flowsol)->VarthetaPhase[nophase]+cellindex));
				
			}
			else
			{
				*((z->Flowsol)->Rho+cellindex) =*((z->Flowsol)->Rho+cellindex) + *((z->Flowsol)->RhoofPhase[nophase]+cellindex*flagrho)*(*((z->Flowsol)->VarthetaPhase[nophase]+cellindex));
				
			}
		
		}
		
		variable=*((z->Flowsol)->Rho+cellindex);
		
	}
	else if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
	{
		variable=*((z->Flowsol)->Rho+cellindex);
		//we are solving equation for the density of carrier phase, which is done just before disperse phase equation.
	}	
	else
	{
		//Write User defined functions here >>>>>
		
		/*variable =  (100000.0 * 0.001)/(*(z->Flowsol->T+cellindex));
		omegaA = *(z->Flowsol->Scal[0]+cellindex);
		omegaB = *(z->Flowsol->Scal[1]+cellindex);
		//MolWtA = 39.944;
		//MolWtB = 2.016;
		MolWtA = 34.29548;
		MolWtB = 25.536;

		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));

		TotMolWt = MolWtA*xA + MolWtB*xB;
		variable = (101325.0 * TotMolWt)/(8314.0*273.0);*/
                //variable = (101325.0 * TotMolWt)/(8314.0*(*(z->Flowsol->T+cellindex)));
		//variable = (100000.0* 0.001)/(*(z->Flowsol->T+cellindex));
		//variable = (100000.0* (28.9/8314.0))/(*(z->Flowsol->T+cellindex));
		
		//For Saclay Test Problem
// 		omegaA = *(z->Flowsol->Scal[2]+cellindex);
// 		omegaB = *(z->Flowsol->Scal[3]+cellindex);
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
In this function user can write the Viscosity in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if viscosity is a function of flow variables like
scalar concentrations etc.
*/
double userfuncViscosity(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable=  0.01;
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
		f = *((z->Flowsol)->f+cellindex);
		//This statement has been written keeping in mind that variation of variable
		//is due to phase change only.  Others conditions can also be written in
		//this block if needed.
		if(base->flagSolidification==1)
		variable = Mu_s + (Mu_f - Mu_s) *f;
		else
		variable = f*Mu_s + (1.0-f)*Mu_f;
	}
	else if(base->flagMPFV==1)//Added by sushil for multiphase flow problems.....
	{	
		*((z->Flowsol)->Viscoeff+cellindex)=0.0;
		for(nophase=0;nophase<nop;nophase++)
		{
			flagvisc=(base->BaseParameters->flagViscoeffPhase[nophase]);
			if(flagvisc==1)
			{
				*((z->Flowsol)->ViscosityofPhase[nophase]+cellindex)=1.0;//Default value
				
				*((z->Flowsol)->ViscoeffDummy[nophase]+cellindex) = (*((z->Flowsol)->ViscosityofPhase[nophase]+cellindex*flagvisc)*(*((z->Flowsol)->VarthetaPhase[nophase]+cellindex)));
				
				
				*((z->Flowsol)->Viscoeff+cellindex) =*((z->Flowsol)->Viscoeff+cellindex)+ *((z->Flowsol)->ViscosityofPhase[nophase]+cellindex*flagvisc)*(*((z->Flowsol)->VarthetaPhase[nophase]+cellindex));
			}
			else
			{
				
				*((z->Flowsol)->ViscoeffDummy[nophase]+cellindex) = (*((z->Flowsol)->ViscosityofPhase[nophase]+cellindex*flagvisc)*(*((z->Flowsol)->VarthetaPhase[nophase]+cellindex)));
				
				*((z->Flowsol)->Viscoeff+cellindex) =*((z->Flowsol)->Viscoeff+cellindex)+ *((z->Flowsol)->ViscosityofPhase[nophase]+cellindex*flagvisc)*(*((z->Flowsol)->VarthetaPhase[nophase]+cellindex));
		
			}
			//printf("\n zone---%s and dummy[%d]=%lf",z->zonename,nophase,*((z->Flowsol)->ViscoeffDummy[nophase]+cellindex));
		}
		
	
		variable=*((z->Flowsol)->Viscoeff+cellindex);
	}
	else
	{
		//Write User defined functions here >>>>>
		
// 		omegaA = *(z->Flowsol->Scal[0]+cellindex);
// 		omegaB = *(z->Flowsol->Scal[1]+cellindex);
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
		
// 		muA= 0.01;
// 		muB= 0.01;
// 		
// 		phiAA = 1.0;
// 		phiBB = 1.0;
// 		phiAB = (1.0/sqrt(8.0))* pow((1.0+(MolWtA/MolWtB)),-0.5) * pow((1.0+pow((muA/muB),0.5)*pow((MolWtB/MolWtA),0.25)),2.0);
// 		phiBA = (1.0/sqrt(8.0))* pow((1.0+(MolWtB/MolWtA)),-0.5) * pow((1.0+pow((muB/muA),0.5)*pow((MolWtA/MolWtB),0.25)),2.0);
// 		
// 		variable = (xA*muA/(xA*phiAA+xA*phiAB))+(xB*muB/(xB*phiBA+xB*phiBB));
		
		//Saclay Test Problem
// 		omegaA = *(z->Flowsol->Scal[2]+cellindex);
// 		omegaB = *(z->Flowsol->Scal[3]+cellindex);
// 		MolWtA = 28.84;
// 		MolWtB = 4.0;
// 		muA = 0.00001855;
// 		muB = 0.00001989;
// 
// 		xA = (omegaA/MolWtA)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 		xB = (omegaB/MolWtB)/((omegaA/MolWtA)+(omegaB/MolWtB));
// 
// 		variable = (xA*muA*sqrt(MolWtA) + xB*muB*sqrt(MolWtB))/(xA*sqrt(MolWtA) + xB*sqrt(MolWtB));
	}
	
	return(variable);
}
/**
In this function user can write the ThermalConductivity in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ThermalConductivity is a function of flow variables like
scalar concentrations etc.
*/
double userfuncThermalConductivity(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable = 0.1;
	//double omegaA,omegaB,MolWtA,MolWtB;
	//double xA,xB,muA,muB,kA,kB;
	//double phiAA,phiAB,phiBB,phiBA;
	double K_s,K_f,f,CP_s,CP_f;
	//for evaporation droplet phase flow
	double Ka,Kv,massfrac;
	
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		K_s = base->BaseParameters->K_s;
		K_f = base->BaseParameters->K_f;
		CP_s = base->BaseParameters->CP_s;
		CP_f = base->BaseParameters->CP_f;
		f = *((z->Flowsol)->f+cellindex);
		//This statement has been written keeping in mind that variation of variable
		//is due to phase change only.  Others conditions can also be written in
		//this block if needed.
		if(base->flagSolidification==1)
		variable = (K_s + (K_f - K_s) *f) / (CP_s + (CP_f - CP_s) *f);
		else
		variable = (f *K_s + (1.0-f) *K_f)/(f * CP_s + (1.0-f)*CP_f);
	}
	else if(base->flagDroplet==1)//shash
	{
		Ka = 0.0242;
		Kv = 0.0128;
	 	massfrac = 0.0;
		//if mass fraction is not being solved, the mixture properties will become the properties of pure gas phase.
		if(base->BaseParameters->flagDropFuelMassFrac==1)
	 	{
	 	massfrac = *(z->Flowsol->DropFuelMassFrac+cellindex);
		variable = massfrac*Kv+(1.0-massfrac)*Ka;
	  	}
		else
		{
		variable = massfrac*Kv+(1.0-massfrac)*Ka;
		}	
	}
	else
	{
		/*omegaA = *(z->Flowsol->Scal[0]+cellindex);
		omegaB = *(z->Flowsol->Scal[1]+cellindex);
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
In this function user can write the SpecificHeat in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if SpecificHeat is a function of flow variables like
scalar concentrations etc.
*/
double userfuncSpecificHeat(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable = 1.0;
	//double omegaA,omegaB,MolWtA,MolWtB;
	//double xA,xB,CpA,CpB,muA,muB;	
	double CP_s,CP_f,f;
	//for evaporation droplet phase flow
	double Cpa,Cpv,massfrac;
	
	if((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		CP_s = base->BaseParameters->CP_s;
		CP_f = base->BaseParameters->CP_f;
		f = *((z->Flowsol)->f+cellindex);
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
	 	massfrac = *(z->Flowsol->DropFuelMassFrac+cellindex);
	 	variable = massfrac*Cpv+(1.0-massfrac)*Cpa;
		}
		else
		{
		variable = massfrac*Cpv+(1.0-massfrac)*Cpa;
		}
	}
	else
	{	
		/*omegaA = *(z->Flowsol->Scal[0]+cellindex);
		omegaB = *(z->Flowsol->Scal[1]+cellindex);
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
In this function user can write the specific heat of liquid droplet for energy equation of the droplet phase in the form of functions. Coordinate X,Y,Z, Zone pointer and cell index is given as input. This function can be used if specific heat of liquid droplet is a function of flow variables like droplet temperature etc.
*/
double userfuncCPdrop(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	
	double Cpf,Tsurf;
	
	
	if(base->flagDroplet==1)//shash
	{
	 Tsurf = *(z->Flowsol->Tdrop+cellindex);
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
double userDefinedCD(double Re)
{
	double Cd;
// 	//Moukalled et al. 2003
	if(Re>1e-9&&Re<0.49)
	{
		Cd=16.0/Re;
				//Cd=(24.0*(1+0.15*pow(Re,0.687)))/Re;
	}
	else if(Re>0.49&&Re<100.0)
		Cd=20.0/(pow(Re,0.643));
	else if(Re>100.0)
		Cd = 6.3/(pow(Re,0.385));
	else
		Cd=0.0;
// 	
	printf("\nMoukalled Cd=%e\n",Cd);
	
	//Schiller and Naumann drag relation
// 	if(Re<=1000.0&&Re!=0.000)
// 		Cd=(24.0*(1+0.15*pow(Re,0.687)))/Re;
// 	else if(Re==0.00)
// 		Cd=0.0000;
// 	else 
// 		Cd=0.44;
// 
// 	printf("\nSchiller and Naumann Cd=%e\n",Cd);
	
// 	//White drag relation
// 	if(Re>0.000)
// 	Cd=0.5+(24.0/Re)+(6.0/(1+sqrt(Re)));
// 	else
// 	Cd=0.000;	
// 	//printf("\nWhilte Cd=%e\n",Cd);
	
// 	//Ishii and Zuber relation
// 	if(Re>0.000)
// 	Cd=(24.0*(1+0.1*pow(Re,0.75)))/Re;
// 	else
// 	Cd=0.000;
// 	
// 	//printf("\nIshii and Zuber Cd=%e\n",Cd);
	
	return(Cd);
		
}
double userDefinedCDDisperse(double Re)
{
	double Cd;
	
	if(Re<=200000.0)
	{
		if(Re < 0.00000001)
			Cd=0;
		else if((Re > 0.00000001 )&&(Re<=1.0))
			Cd=24.0/Re;
		else if((Re > 1.0 )&&(Re<=1000.0))
			Cd=(24.0/Re)*(1+0.15*pow(Re,0.687));
		else if((Re>1000.0)&&(Re<=200000.0))
        		Cd=0.44;
       	}
	else
		Cd=0.09;
	
	return(Cd);
		
}

double userDefinedDropPfYf(double Tsurf,double P)
{
	double Pfs,Yfs,MwtA,MwtF;
	double Tbn;
	Tbn = 371.4; 
	MwtA=28.97;//air 
	MwtF=100.204;//n-heptane
	if(Tsurf>Tbn)
	{
	Pfs = exp(14.2146-(3151.68/(Tsurf-43.0)));
	}
	else
	{
	Pfs = exp(14.3896-(3209.45/(Tsurf-43.0)));
	}
	Yfs=pow((1+(((P/Pfs)-1)*(MwtA/MwtF))),-1);
	
	return(Yfs);
}

double userDefinedDropCpa(double Tref)
{
	double Cpa,c1,c2,c3,c4,c5;
	double MwtA;
	MwtA=28.97;
	
	//for air
	c1=0.2896*1e5;
	c2=0.09390*1e5;
	c3=3.0120*1e3;
	c4=0.0758*1e5;
	c5=1484.0;
	
	Cpa = c1+c2*((c3/Tref)/(sinh(c3/Tref)))*((c3/Tref)/(sinh(c3/Tref)))+c4*((c5/Tref)/(cosh(c5/Tref)))*((c5/Tref)/(cosh(c5/Tref)));
	
	Cpa=Cpa/MwtA;
	
	return(Cpa);
}

double userDefinedDropCpv(double Tref)
{
	double Cpv;
	
	Cpv = (363+0.467*Tref)*(5-0.001*883.05);

	return(Cpv);
}


double userDefinedDropKa(double Tref)
{
	double Ka;
	
	Ka=1.5207*1e-11*Tref*Tref*Tref-4.8574*1e-8*Tref*Tref+1.0184*1e-4*Tref-3.9333*1e-4; 

	return(Ka);
}
double userDefinedDropKv(double Tref,double Cpv)
{
	double Kv,Tcri,Pcri,Tr,lambda,MwtF;
	
	MwtF=100.204;//n-heptane
	Tcri=540.17;//n-heptane
	Pcri=2631.633*1000.0;//n-heptane(Pa)
	
	Tr=Tref/Tcri;
	
	lambda = pow(Tcri,(1.0/6.0))*pow(MwtF,(1.0/2.0))*pow((101.325/(Pcri/1000.0)),(2.0/3.0));
	
	Kv=1e-7*pow((14.52*Tr-5.14),(2.0/3.0))*(Cpv*MwtF/lambda);
 

	return(Kv);
}

double userDefinedDropL(double Tsurf)
{
	double L,Tcri,Tbn;
	Tcri=540.17;//n-heptane
	Tbn = 371.4;
	
	L=317.8*1000.0*pow(((Tcri-Tsurf)/(Tcri-Tbn)),-0.38);

	return(L);
}
double userDefinedDropCpf(double Tsurf)
{
	double Cpf;
	
	//n-heptane
	Cpf = 13058.45044066-126.5095282565*Tsurf+0.5279613848638*Tsurf*Tsurf-0.0009457386295687*Tsurf*Tsurf*Tsurf+6.369853422618e-7*Tsurf*Tsurf*Tsurf*Tsurf;


	return(Cpf);
}

double userDefinedDropRhoF(double Tsurf)
{
	double RhoF;
	
	//n-heptane
	if(Tsurf<=538.0)
		RhoF = -941.03+19.96181*Tsurf-0.08612051*Tsurf*Tsurf+1.579494e-4*Tsurf*Tsurf*Tsurf-1.089345e-7*Tsurf*Tsurf*Tsurf*Tsurf;
	else
		RhoF=4.195281e7-2.360524e5*Tsurf+442.7316*Tsurf*Tsurf-0.2767921*Tsurf; 


	return(RhoF);
}

/**
In this function user can write the ScalCoeff1 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff1 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff1(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 1)
	{
		printf("Exiting from userfuncScalCoeff1 \n");
		exit(88);
	}
	variable = 0.0;
	return(variable);
}
/**
In this function user can write the ScalCoeff2 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff2 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff2(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 2)
	{
		printf("Exiting from userfuncScalCoeff2 \n");
		exit(88);
	}
	variable = 0.0;
	return(variable);
}
/**
In this function user can write the ScalCoeff3 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff3 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff3(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 3)
	{
		printf("Exiting from userfuncScalCoeff3 \n");
		exit(88);
	}
	//For Species/Sclars in Turbulent flow used for Saclay Test Problem.
// 	variable = 7.0e-05 + *(z->Flowsol->Tcoeff->TViscoeff+cellindex)/0.7;
	variable = 0.0;
	return(variable);
}
/**
In this function user can write the ScalCoeff4 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff4 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff4(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	//NOTE: Keep this message as it warns, in case of turbulence, of unwanted updation of scalar diffusion coefficients
	if(base->BaseParameters->Turb->numTurbScalars >= 4)
	{
		printf("Exiting from userfuncScalCoeff4 \n");
		exit(88);
	}
	//For Species/Sclars in Turbulent flow used for Saclay Test Problem.
// 	variable = 7.0e-05 + *(z->Flowsol->Tcoeff->TViscoeff+cellindex)/0.7;
	variable = 0.0;
	return(variable);
}
/**
In this function user can write the ScalCoeff5 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff5 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff5(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
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
In this function user can write the ScalCoeff6 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff6 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff6(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function user can write the ScalCoeff7 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff7 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff7(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function user can write the ScalCoeff8 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff8 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff8(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function user can write the ScalCoeff9 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff9 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff9(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}
/**
In this function user can write the ScalCoeff10 in the form of functions. Coordinate X,Y,Z, Zone pointer and 
cell index is given as input. This function can be used if ScalCoeff10 is a function of flow variables like
scalar concentrations etc.
*/
double userfuncScalCoeff10(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=6.0*Y*(1.0-Y);
	return(variable);
}

double userfuncChargeDiffusion(Zone *z,int cellindex,double X,double Y,double Z,Base *base)
{
	double variable;
	variable=5.32e-06;
	return(variable);
}


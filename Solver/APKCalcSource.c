#include "APKCalcSource.h"
/**
The function CalcSource() is used to calculate the various source terms occuring in the
governing equations of various variables. This function, hence, takes the variable tag
as input and depending on the variable tag, either computes or calls necessary functions
to compute the source terms corresponding to that particular variable.
*/

//#define Ra 1000.0
//#define Pr 0.71

double CalcSource(double *source,int variableTag,Zone *zone,Base *base,int index[3])
{
	int NI,NJ,NK,k,j,i,indexP,face;
	int *ptrSize,flagrho,flagspecificheat;
	double diagonalCoef;
	double *slabPressure[6],*slabVolume[6],*slabrho[6],*slabKE[6],*density;
	double *SurfArea,*pressure,*ptrTemp,*ptrVolume;
	double gXcomp,gYcomp,gZcomp,energySource;//,density;
	double *massFlux,rho,velP;
	double sumMassFlux,Tref,beta;
	
	double dele_prev,dele,delt,A,liq_frac,C,C_p,P_o;

        double FMassFrac,OxMassFrac,PMassFrac,*KE,epsilon,MinFrac,Stoich,MinProd,rhoinfty; //Shyam
	double *Potential, *slabPotential[6], Charge, DiElec,volP, Scaling;//Venky
	
	double B_x, B_y, B_z, velPu, velPv, velPw, sigma, Dphi[3];//naveen
	
	double K,Kold,Omega,OmegaOld,TProduction,betastar,mue_t,Dissipation=0.0,diacoeff=0.0,gamma,beta1,vol;
	double sumGradSquare,MuEffK,MuEffEps,Eps,Epsold,Fsum1=0.0;
	double C1epsi,C2epsi,Gamma,C0epsi;
	double CrossDiffusion=0,F1,Radflux;
	double Sigma=5.67e-8,*TempP,*Kp;
	
	*source=0.0;
	diagonalCoef=0.0;
	
	ptrSize= zone->Size;
	delt = base->BaseParameters->TimeStep;
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	MakeSlabPtrArray(zone,slabrho,710);//shyam
	MakeSlabPtrArray(zone,slabKE,600);//shyam
	
	ptrVolume= zone->ZGeom->Volume;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	k=index[0],j=index[1],i=index[2];
	indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	flagrho = zone->Flowsol->flagRho;
	flagspecificheat = zone->Flowsol->flagSpecificHeat;
	rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
        density = zone->PrevFlowsol->Rho+indexP*flagrho; //Shyam
        KE = zone->PrevFlowsol->Scal[0]+indexP; //Shyam
	P_o = base->BaseParameters->P_o;//reference pressure to handle {(P/Rho).(del.U)} term
	if(base->flagLMMHD==1)
	{
		FindgradientPhiMHD(base,zone,index,Dphi);//naveen
		B_x = *(zone->Flowsol->MfieldX+indexP);
		B_y = *(zone->Flowsol->MfieldY+indexP);
		B_z = *(zone->Flowsol->MfieldZ+indexP);
		sigma = *(zone->Flowsol->Eleconductivity+indexP);
	}
		
	if((variableTag==600 || variableTag==601) && base->BaseParameters->Turb->numTurbScalars >= 2)
	{
		vol = *(zone->ZGeom->Volume+indexP);
		mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+indexP);
		delt = base->BaseParameters->Turb->DualTimeStep;//delt=delt/10.0;
		
		if(base->flagTurbModel==2)
		{
			Kold = *(zone->PrevFlowsol->Scal[0]+indexP);
			OmegaOld = *(zone->PrevFlowsol->Scal[1]+indexP);
			K        = *(zone->Flowsol->Scal[0]+indexP);
			Omega    = *(zone->Flowsol->Scal[1]+indexP);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
			betastar = (zone->Flowsol->Tcoeff->betastar);
			
			gamma = zone->Flowsol->Tcoeff->gamma;
 			beta1  = zone->Flowsol->Tcoeff->beta;
		}
		else if(base->flagTurbModel==6)
		{
			Kold = *(zone->PrevFlowsol->Scal[0]+indexP);
			OmegaOld = *(zone->PrevFlowsol->Scal[1]+indexP);
			K        = *(zone->Flowsol->Scal[0]+indexP);
			Omega    = *(zone->Flowsol->Scal[1]+indexP);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
			betastar = (zone->Flowsol->Tcoeff->betastar);
			F1= *(zone->Flowsol->Tcoeff->F1+indexP);
		}
		else
		{
			Kold   = *(zone->PrevFlowsol->Scal[0]+indexP);
			Epsold = *(zone->PrevFlowsol->Scal[1]+indexP);
			K   = *(zone->Flowsol->Scal[0]+indexP);
			Eps = *(zone->Flowsol->Scal[1]+indexP);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
			
			//(Epsold/Kold) is calculated outside the dual time step and stored in Gamma for use through the dual time step
// 			Gamma = *(zone->Flowsol->Gamma+indexP);
			Gamma = (Epsold/Kold);

// 			if(Gamma > 1.0e5)
// 				printf("Zone: %s k=%d j=%d i=%d Epsold=%e Kold=%e Gamma=%e Epsold/Kold=%e ...\n", zone->zonename,k,j,i,Epsold,Kold,Gamma,Epsold/Kold);
			
// 			if(Gamma < 0)
// 			{
// 				printf("Zone: %s k=%d j=%d i=%d Epsold=%e Kold=%e Gamma=%e limiting Gamma to 0 ...\n", zone->zonename,k,j,i,Epsold,Kold,Gamma);
// 				Gamma = 0;	//This is to stabilize the equations by ensuring the positivity of coeffients of Reaction terms (Production, Dissipation etc.)
// 			}
		}
	}

	switch(variableTag)
	{
		case 1:
		{
			SurfArea=zone->ZGeom->SurfaceArea;
			MakeSlabPtrArray(zone,slabPressure,4);/*Called for making slabs of Pressure*/
			pressure=zone->Flowsol->P;
			
			*source=(1.0/rho)*CalcPressureGrad(pressure,slabPressure,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);
			//printf("source=%e\n",*(source));

// 			if(base->BaseParameters->Turb->numTurbScalars >= 2 && base->flagTurbComb==0) //Shyam, Only for non-reactive cases
// 			{
// 				*source= *source + (2.0/3.0)*CalcPressureGrad(KE,slabKE,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				if(flagrho == 1)
// 				{
// 					*source= *source + (2.0/3.0)*(*KE/rho)*CalcPressureGrad(density,slabrho,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				}
// 			}
			
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->u+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			//printf("source=%e\t",*(source));
			if(base->flagSolidification==1)
			{
				liq_frac = *(zone->Flowsol->f+indexP);
				C = base->BaseParameters->C;
				A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
				A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
				A  = A/rho;
				//printf("A= %e\tliq_frac=%e\n",A,liq_frac);
				diagonalCoef = diagonalCoef - A;
				*(source) = *(source) + A * *(zone->Flowsol->u+indexP);
				//printf("source=%e\n",*(source));
			}
			
			if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
				*source = *source + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);
			
			if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==1)
			{
				if(base->flagTurbComb==1) //Shyam
				{
				gXcomp=base->BaseParameters->gravityX;
				rhoinfty = *(base->BaseParameters->rhoinfinity);
				*source=*source+(*(ptrVolume+indexP))*gXcomp*(rho-rhoinfty)/rho;
				}
				else
				{
				gXcomp=base->BaseParameters->gravityX;
				*source=*source+(*(ptrVolume+indexP))*gXcomp;
				//*source=*source+(*(ptrVolume+indexP))*gXcomp*density;
				}
			}
			if(base->BaseParameters->flagBuoyancy==1)
			{
				
				Tref=base->BaseParameters->TReference;
				beta=base->BaseParameters->VolumetricExpCoeff;
				gXcomp=base->BaseParameters->gravityX;
				if(base->flagSolidification ==1)
				{
					ptrTemp=zone->Flowsol->enthalpy;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
				}
				else
				{
					ptrTemp=zone->Flowsol->T;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
					//C_p = 1.0;
				}
				*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP) - Tref)/C_p);
				//*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP)/C_p) - Tref);
			}
			
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_u+indexP));
			}
			
			//This for Electric field Source Term			
			if((base->flagEfield==1) || (base->flagEfield==2))
			{
				//Calculation of Grad(Phi) just like Grad(pressure)
				
				MakeSlabPtrArray(zone,slabPotential,8);
				Potential=zone->Flowsol->Potential;
				Charge = *(zone->Flowsol->RhoIon+indexP);// Get Ion charge Density
				Scaling = base->BaseParameters->ScaleFactor;
				*source = *source + Scaling*(Charge/rho)*CalcPressureGrad(Potential,slabPotential,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);
			}
			
			//This is for LMMHD problems(naveen). The function FindgradientPhiMHD() is in solver/LMMHD
			if(base->flagLMMHD==1)
			{
				velPu = *(zone->Flowsol->u+indexP);
				velPv = *(zone->Flowsol->v+indexP);
				velPw = *(zone->Flowsol->w+indexP);
				*source = *source + (*(ptrVolume+indexP))*(sigma/rho)*((B_y*Dphi[2] - B_z*Dphi[1] + B_x*(velPv*B_y+velPw*B_z))-((B_y*B_y+B_z*B_z)*velPu));
				diagonalCoef = diagonalCoef - *(ptrVolume+indexP)*(sigma/rho) *(B_y*B_y+B_z*B_z) ;
			}
		}
		break;
		case 100:
		{
			SurfArea=zone->ZGeom->SurfaceArea;
		
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->u+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = velP*sumMassFlux;

// 			if(base->BaseParameters->Turb->numTurbScalars >= 2 && base->flagTurbComb==0) //Shyam, Only for non-reactive cases
// 			{
// 				*source= *source + (2.0/3.0)*CalcPressureGrad(KE,slabKE,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				if(flagrho == 1)
// 				{
// 					*source= *source + (2.0/3.0)*(*KE/rho)*CalcPressureGrad(density,slabrho,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				}
// 			}
			
			if(base->flagSolidification==1)
			{
				liq_frac = *(zone->Flowsol->f+indexP);
				C = base->BaseParameters->C;
				A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
				A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
				A  = A/rho;
				diagonalCoef = diagonalCoef - A;
				*(source) = *(source) + A *  *(zone->Flowsol->umass+indexP);
				
			}
			
			if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
			{
				*source = *source + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);
				//printf("findgradient value = %e\n",*source);
			}
			
			if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==1)
			{
				if(base->flagTurbComb==1) //Shyam
				{
				gXcomp=base->BaseParameters->gravityX;
				rhoinfty = *(base->BaseParameters->rhoinfinity);
				*source=*source+(*(ptrVolume+indexP))*gXcomp*(rho-rhoinfty)/rho;
				}
				else
				{
				gXcomp=base->BaseParameters->gravityX;
				*source=*source+(*(ptrVolume+indexP))*gXcomp;
				//*source=*source+(*(ptrVolume+indexP))*gXcomp*density;
				}
			}
			if(base->BaseParameters->flagBuoyancy==1)
			{
				
				Tref=base->BaseParameters->TReference;
				beta=base->BaseParameters->VolumetricExpCoeff;
				gXcomp=base->BaseParameters->gravityX;
				if(base->flagSolidification ==1)
				{
					ptrTemp=zone->Flowsol->enthalpy;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
				}
				else
				{
					ptrTemp=zone->Flowsol->T;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
					//C_p = 1.0;
				}
				*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP) - Tref)/C_p);
				//*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP)/C_p) - Tref);
			}
			
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_u+indexP));
			}
			
			//This for Electric field Source Term			
			if((base->flagEfield==1) || (base->flagEfield==2))
			{
				//Calculation of Grad(Phi) just like Grad(pressure)
				
				MakeSlabPtrArray(zone,slabPotential,8);
				Potential=zone->Flowsol->Potential;
				Charge = *(zone->Flowsol->RhoIon+indexP);// Get Ion charge Density
				Scaling = base->BaseParameters->ScaleFactor;
				*source = *source + Scaling*(Charge/rho)*CalcPressureGrad(Potential,slabPotential,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);

			}
			//This is for LMMHD problems(naveen)
			if(base->flagLMMHD==1)
			{
				velPu = *(zone->PrevFlowsol->u+indexP);
				velPv = *(zone->PrevFlowsol->v+indexP);
				velPw = *(zone->PrevFlowsol->w+indexP);
				*source = *source + (*(ptrVolume+indexP))*(sigma/rho)*(B_y*Dphi[2] - B_z*Dphi[1] + B_x*(velPv*B_y+velPw*B_z) - velPu*(B_y*B_y + B_z*B_z));
			}
		}
		break;
		case 2:
		{
			SurfArea=zone->ZGeom->SurfaceArea;
			MakeSlabPtrArray(zone,slabPressure,4);/*Called for making slabs of Pressure*/
			pressure=zone->Flowsol->P;
			
			*source=(1.0/rho)*CalcPressureGrad(pressure,slabPressure,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);

// 			if(base->BaseParameters->Turb->numTurbScalars >= 2 && base->flagTurbComb==0) //Shyam, Only for non-reactive cases
// 			{
// 				*source= *source + (2.0/3.0)*CalcPressureGrad(KE,slabKE,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				if(flagrho == 1)
// 				{
// 					*source= *source + (2.0/3.0)*(*KE/rho)*CalcPressureGrad(density,slabrho,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				}
// 			}
			
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->v+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			if(base->flagSolidification==1)
			{
				liq_frac = *(zone->Flowsol->f+indexP);
				C = base->BaseParameters->C;
				A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
				A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
				A  = A/rho;
				diagonalCoef = diagonalCoef - A;
				*(source) = *(source) + A * *(zone->Flowsol->v+indexP);
			}
			if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
				*source = *source + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);

			if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==1)
			{
				if(base->flagTurbComb==1) //Shyam
				{
				gYcomp=base->BaseParameters->gravityY;
				rhoinfty = *(base->BaseParameters->rhoinfinity);
				*source=*source+(*(ptrVolume+indexP))*gYcomp*(rho-rhoinfty)/rho;
				}
				else
				{
				gYcomp=base->BaseParameters->gravityY;
				*source=*source+(*(ptrVolume+indexP))*gYcomp;
				//*source=*source+(*(ptrVolume+indexP))*gYcomp*density;
				}
			}
			if(base->BaseParameters->flagBuoyancy==1)
			{
				
				Tref=base->BaseParameters->TReference;
				beta=base->BaseParameters->VolumetricExpCoeff;
				gYcomp=base->BaseParameters->gravityY;
				if(base->flagSolidification ==1)
				{
					ptrTemp=zone->Flowsol->enthalpy;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
				}
				else
				{
					ptrTemp=zone->Flowsol->T;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
					//C_p = 1.0;
				}
				*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP) - Tref)/C_p);
				//*source=*source-(*(ptrVolume+indexP))*beta*gYcomp*((*(ptrTemp+indexP)/C_p) - Tref);
			}
			
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_v+indexP));
			}
			
			//This for Electric field Source Term			
			if((base->flagEfield==1) || (base->flagEfield==2))
			{
				//Calculation of Grad(Phi) just like Grad(pressure)
				
				MakeSlabPtrArray(zone,slabPotential,8);
				Potential=zone->Flowsol->Potential;
				Charge = *(zone->Flowsol->RhoIon+indexP);// Get Ion charge Density
				Scaling = base->BaseParameters->ScaleFactor;
				*source = *source + Scaling*(Charge/rho)*CalcPressureGrad(Potential,slabPotential,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);

			}

			//This is for LMMHD problems(naveen)
			if(base->flagLMMHD==1)
			{
				velPu = *(zone->Flowsol->u+indexP);
				velPv = *(zone->Flowsol->v+indexP);
				velPw = *(zone->Flowsol->w+indexP);
				*source = *source + (*(ptrVolume+indexP))*(sigma/rho)*((B_z*Dphi[0] - B_x*Dphi[2] + B_y*(velPu*B_x+velPw*B_z))-((B_x*B_x+B_z*B_z)*velPv));
				diagonalCoef = diagonalCoef - *(ptrVolume+indexP)*(sigma/rho) *(B_x*B_x+B_z*B_z) ;
			}
		}
		break;
		case 200:
		{
			SurfArea=zone->ZGeom->SurfaceArea;
		
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->v+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = velP*sumMassFlux;
			
// 			if(base->BaseParameters->Turb->numTurbScalars >= 2 && base->flagTurbComb==0) //Shyam, Only for non-reactive cases
// 			{
// 				*source= *source + (2.0/3.0)*CalcPressureGrad(KE,slabKE,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				if(flagrho == 1)
// 				{
// 					*source= *source + (2.0/3.0)*(*KE/rho)*CalcPressureGrad(density,slabrho,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				}
// 			}
			
			if(base->flagSolidification==1)
			{
				liq_frac = *(zone->Flowsol->f+indexP);
				C = base->BaseParameters->C;
				A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
				A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
				A  = A/rho;
				diagonalCoef = diagonalCoef - A;
				*(source) = *(source) + A * *(zone->Flowsol->vmass+indexP);
			}
			if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
				*source = *source + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);
			
			if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==1)
			{
				if(base->flagTurbComb==1) //Shyam
				{
				gYcomp=base->BaseParameters->gravityY;
				rhoinfty = *(base->BaseParameters->rhoinfinity);
				*source=*source+(*(ptrVolume+indexP))*gYcomp*(rho-rhoinfty)/rho;
				}
				else
				{
				gYcomp=base->BaseParameters->gravityY;
				*source=*source+(*(ptrVolume+indexP))*gYcomp;
				//*source=*source+(*(ptrVolume+indexP))*gYcomp*density;
				}
			}
			if(base->BaseParameters->flagBuoyancy==1)
			{
				
				Tref=base->BaseParameters->TReference;
				beta=base->BaseParameters->VolumetricExpCoeff;
				gYcomp=base->BaseParameters->gravityY;
				if(base->flagSolidification ==1)
				{
					ptrTemp=zone->Flowsol->enthalpy;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
				}
				else
				{
					ptrTemp=zone->Flowsol->T;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
					//C_p = 1.0;
				}
				*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP) - Tref)/C_p);
				//*source=*source-(*(ptrVolume+indexP))*beta*gYcomp*((*(ptrTemp+indexP)/C_p) - Tref);
			}
			
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_v+indexP));
			}
			
			//This for Electric field Source Term
			if((base->flagEfield==1) || (base->flagEfield==2))
			{
				//Calculation of Grad(Phi) just like Grad(pressure)
				
				MakeSlabPtrArray(zone,slabPotential,8);
				Potential=zone->Flowsol->Potential;
				Charge = *(zone->Flowsol->RhoIon+indexP);// Get Ion charge Density
				Scaling = base->BaseParameters->ScaleFactor;
				*source = *source + Scaling*(Charge/rho)*CalcPressureGrad(Potential,slabPotential,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);
			}
			
			//This is for LMMHD problems(naveen)
			if(base->flagLMMHD==1)
			{
				velPu = *(zone->PrevFlowsol->u+indexP);
				velPv = *(zone->PrevFlowsol->v+indexP);
				velPw = *(zone->PrevFlowsol->w+indexP);
				*source = *source + (*(ptrVolume+indexP))*(sigma/rho)*(B_z*Dphi[0] - B_x*Dphi[2] + B_y*(velPu*B_x+velPw*B_z) - velPv*(B_x*B_x + B_z*B_z));
			}
		}
		break;
		case 3:
		{
			SurfArea=zone->ZGeom->SurfaceArea;
			MakeSlabPtrArray(zone,slabPressure,4);/*Called for making slabs of Pressure*/
			pressure=zone->Flowsol->P;
			
			*source=(1.0/rho)*CalcPressureGrad(pressure,slabPressure,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag);
			
// 			if(base->BaseParameters->Turb->numTurbScalars >= 2 && base->flagTurbComb==0) //Shyam, Only for non-reactive cases
// 			{
// 				*source= *source + (2.0/3.0)*CalcPressureGrad(KE,slabKE,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				if(flagrho == 1)
// 				{
// 					*source= *source + (2.0/3.0)*(*KE/rho)*CalcPressureGrad(density,slabrho,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				}
// 			}
			
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->w+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			
			if(base->flagSolidification==1)
			{
				liq_frac = *(zone->Flowsol->f+indexP);
				C = base->BaseParameters->C;
				A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
				A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
				A  = A/rho;
				diagonalCoef = diagonalCoef - A;
				*(source) = *(source) + A * *(zone->Flowsol->w+indexP);
			}
			if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
				*source = *source + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);

			if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==1)
			{
				if(base->flagTurbComb==1) //Shyam
				{
				gZcomp=base->BaseParameters->gravityZ;
				rhoinfty = *(base->BaseParameters->rhoinfinity);
				*source=*source+(*(ptrVolume+indexP))*gZcomp*(rho-rhoinfty)/rho;
				}
				else
				{
				gZcomp=base->BaseParameters->gravityZ;
				*source=*source+(*(ptrVolume+indexP))*gZcomp;
				//*source=*source+(*(ptrVolume+indexP))*gYcomp*density;
				}
			}
			if(base->BaseParameters->flagBuoyancy==1)
			{
				
				Tref=base->BaseParameters->TReference;
				beta=base->BaseParameters->VolumetricExpCoeff;
				gZcomp=base->BaseParameters->gravityZ;
				if(base->flagSolidification ==1)
				{
					ptrTemp=zone->Flowsol->enthalpy;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
				}
				else
				{
					ptrTemp=zone->Flowsol->T;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
					//C_p = 1.0;
				}
				*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP) - Tref)/C_p);
				//*source=*source-(*(ptrVolume+indexP))*beta*gZcomp*((*(ptrTemp+indexP)/C_p) - Tref);
			}
			
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_w+indexP));
			}
			
			//This is for LMMHD problems(naveen)
			if(base->flagLMMHD==1)
			{
				velPu = *(zone->Flowsol->u+indexP);
				velPv = *(zone->Flowsol->v+indexP);
				velPw = *(zone->Flowsol->w+indexP);
				*source = *source + (*(ptrVolume+indexP))*(sigma/rho)*((B_x*Dphi[1] - B_y*Dphi[0] + B_z*(velPu*B_x+velPv*B_y))-((B_x*B_x+B_y*B_y)*velPw));
				diagonalCoef = diagonalCoef - *(ptrVolume+indexP)*(sigma/rho) *(B_x*B_x+B_y*B_y) ;
			}
		}
		break;
		case 300:
		{
			SurfArea=zone->ZGeom->SurfaceArea;
		
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->w+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = velP*sumMassFlux;
			
// 			if(base->BaseParameters->Turb->numTurbScalars >= 2 && base->flagTurbComb==0) //Shyam, Only for non-reactive cases
// 			{
// 				*source= *source + (2.0/3.0)*CalcPressureGrad(KE,slabKE,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				if(flagrho == 1)
// 				{
// 					*source= *source + (2.0/3.0)*(*KE/rho)*CalcPressureGrad(density,slabrho,ptrVolume,slabVolume,SurfArea,index,ptrSize,variableTag); //Shyam
// 				}
// 			}
			
			if(base->flagSolidification==1)
			{
				liq_frac = *(zone->Flowsol->f+indexP);
				C = base->BaseParameters->C;
				A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
				A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
				A  = A/rho;
				diagonalCoef = diagonalCoef - A;
				*(source) = *(source) + A * *(zone->Flowsol->wmass+indexP);
			}
			
			if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
				*source = *source + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);

			if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==1)
			{
				if(base->flagTurbComb==1) //Shyam
				{
				gZcomp=base->BaseParameters->gravityZ;
				rhoinfty = *(base->BaseParameters->rhoinfinity);
				*source=*source+(*(ptrVolume+indexP))*gZcomp*(rho-rhoinfty)/rho;
				}
				else
				{
				gZcomp=base->BaseParameters->gravityZ;
				*source=*source+(*(ptrVolume+indexP))*gZcomp;
				//*source=*source+(*(ptrVolume+indexP))*gYcomp*density;
				}
			}
			if(base->BaseParameters->flagBuoyancy==1)
			{
				
				Tref=base->BaseParameters->TReference;
				beta=base->BaseParameters->VolumetricExpCoeff;
				gZcomp=base->BaseParameters->gravityZ;
				if(base->flagSolidification ==1)
				{
					ptrTemp=zone->Flowsol->enthalpy;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
				}
				else
				{
					ptrTemp=zone->Flowsol->T;
					C_p =*(zone->Flowsol->SpecificHeat+indexP);
					//C_p = 1.0;
				}
				*source=*source-(*(ptrVolume+indexP))*beta*gXcomp*((*(ptrTemp+indexP) - Tref)/C_p);
				//*source=*source-(*(ptrVolume+indexP))*beta*gZcomp*((*(ptrTemp+indexP)/C_p)- Tref);
			}
			
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_w+indexP));
			}
			
			
			//This is for LMMHD problems(naveen)
			if(base->flagLMMHD==1)
			{
				velPu = *(zone->Flowsol->u+indexP);
				velPv = *(zone->Flowsol->v+indexP);
				velPw = *(zone->Flowsol->w+indexP);
				*source = *source + (*(ptrVolume+indexP))*(sigma/rho)*(B_x*Dphi[1] - B_y*Dphi[0] + B_z*(velPu*B_x+velPv*B_y) - velPw*(B_x*B_x + B_y*B_y));
			}
		}
		break;
		case 6:
		{
			energySource=0.0;
			if(zone->flagHeatGeneration==1)
			{
				if (zone->flagHeat_con_var==0)
				energySource = *((zone->Flowsol)->Energysource);
				else if (zone->flagHeat_con_var==1)
				energySource = *((zone->Flowsol)->Energysource+indexP);
				else
				{
					printf("Wrong tag is supplied for EnergySource in CalcSource.c .....\n So exiting...");
					exit(0);
				}
			}
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->T+indexP)*(*(zone->PrevFlowsol->SpecificHeat+indexP*flagspecificheat));
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			
			if(zone->flagHeatGeneration == 1)
			{
				*source=*source + (1.0/rho)*energySource*(*(ptrVolume+indexP));
				//*source=*source + energySource*(*(ptrVolume+indexP));
			}
			if(base->BaseParameters->flagflowRadiation==1)
			{
				TempP=zone->PrevFlowsol->T;
				Kp=zone->Flowsol->kappa;
				Radflux=*(zone->Flowsol->RadEnergysource+indexP);
				//printf("var=%lf",Radflux);
				*source=*source - (1.0/rho)*Radflux*(*(ptrVolume+indexP));
				
				diagonalCoef=diagonalCoef+(1.0/rho)*(*Kp)*(4.0*Sigma*pow(*(TempP+indexP),3.0))*(*(ptrVolume+indexP));
			}
				
				
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_T+indexP));
			}
			if(base->BaseParameters->flagDissipation==1)
			{
				printf("\nPlease write the function for dissipation in CalcSourece, Function not yet written hence exiting from CalcSource.c\n");//A function to be written for Calculating Dissipation here
				exit(0);
			}
		}
		break;
		case 7:
		{
			energySource=0.0;

			if(zone->flagHeatGeneration==1)
			{
				if(zone->flagHeat_con_var==0)
					energySource = *((zone->Flowsol)->Energysource);
				else if(zone->flagHeat_con_var==1)
					energySource = *((zone->Flowsol)->Energysource+indexP);
				else
				{
					printf("Wrong tag is supplied for EnergySource in CalcSource.c .....\n So exiting...");
					exit(0);
				}
			}

			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->enthalpy+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			if(base->flagTwoPhase==1)
			{
			
			*source = *source -(P_o/rho)*sumMassFlux;
			}
			
			if(base->flagSolidification==1)
			{
				dele 	 = *(zone->Flowsol->Latent_Heat_Content+indexP);
				dele_prev= *(zone->PrevFlowsol->Latent_Heat_Content+indexP);
				*source  = *source - ((dele - dele_prev)*(*(ptrVolume+indexP))/delt);
				
			}
			if(base->BaseParameters->flagEnergySource==1 && zone->flagHeatGeneration == 1)
			{
				*source=*source + (1.0/rho)*energySource*(*(ptrVolume+indexP));
				//*source=*source + energySource*(*(ptrVolume+indexP));
			}
			if(base->BaseParameters->TimeStepingSchemeTag==1)
			{
				*source=*source+(*(zone->Flowsol->PrevTStepFluxes_T+indexP));
			}
			if(base->BaseParameters->flagDissipation==1)
			{
				printf("\nPlease write the function for dissipation in CalcSourece, Function not yet written hence exiting from CalcSource.c\n");//A function to be written for Calculating Dissipation here
				exit(0);
			}
		}
		break;
		case 8:
		{
			velP = *(zone->Flowsol->RhoIon+indexP); //Get the ion charge density
			volP = *(zone->ZGeom->Volume+indexP);   //Get the volume of current cell
			DiElec = *(zone->Flowsol->DiElecConst); //Get the di-electric constant
			*source = *source + (velP*volP)/DiElec;
			//printf("\nCharge in source = %e", velP);
		}
		break;
		case 9:
		{
			*source=0.0;
			diagonalCoef=0.0;
		}
		break;

                case 15: //Shyam
		{
	                MinFrac = 0.0;
                        Stoich = *(base->BaseParameters->Stoich);
			massFlux = zone->Flowsol->Massflux+indexP*6;
			FMassFrac = *(zone->PrevFlowsol->FuelMassFrac+indexP);
                        OxMassFrac = *(zone->PrevFlowsol->OxiMassFrac+indexP);
                        PMassFrac = *(zone->PrevFlowsol->ProdMassFrac+indexP);
                        KE = (zone->PrevFlowsol->Scal[0]+indexP);
                        epsilon = *(zone->PrevFlowsol->Scal[1]+indexP);
                        MinProd = 0.5*PMassFrac/(1.0+Stoich);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + FMassFrac*sumMassFlux;

                        MinFrac = FMassFrac;
                        if(MinFrac > (OxMassFrac/Stoich))
                        MinFrac = OxMassFrac/Stoich;

                        if(MinFrac > MinProd)
                        MinFrac = MinProd;

                        *source = *source - 4.0*(epsilon/(*KE))*MinFrac*(*(ptrVolume+indexP));
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;

		case 14: //Shyam
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->MixFrac+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;             

		case 17:
                {

                        *source = *(zone->Flowsol->PotentialMHDpoissonRHS+indexP);
                        diagonalCoef=0.0;
                }
                break;

		case 600:
		{
			if(base->flagTurbModel==8)
			{
				diagonalCoef = SpalartAllmarasSource(base,zone,source,index);
			}
			else
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(zone->PrevFlowsol->Scal[0]+indexP);
				sumMassFlux=0.0;
				for(face=0;face<6;face++)
				{
					sumMassFlux = sumMassFlux + (*(massFlux+face));
				}
				*source = *source + velP*sumMassFlux;
				
				if(base->BaseParameters->Turb->numTurbScalars >= 1)
				{
					if(base->flagTurbModel==2 || base->flagTurbModel==6)
					{
						*source = *source + mue_t * TProduction;
						diacoeff    = vol * betastar * OmegaOld;
						Dissipation = diacoeff * K;
	
						*source=*source-(Dissipation);
					}
					else
					{
						if(base->flagLogKEps == 1)
						{
							sumGradSquare = CalcGradSquare(base,zone,index,variableTag);
							MuEffK = *(zone->Flowsol->Scalcoeff[0]+indexP);
		
							*source = *source + vol * (MuEffK/rho) * sumGradSquare;
							*source = *source + mue_t * TProduction * exp(-Kold) * (1.0+Kold);
							Dissipation = vol * exp(Epsold - Kold) * (1.0+Kold);
							diacoeff = - (vol * exp(Epsold - Kold)) + mue_t * TProduction * exp(-Kold);
							Fsum1 = Fsum1 + diacoeff * K;
						}
						else
						{
							*source = *source + mue_t * TProduction;
							diacoeff    = vol * Gamma;
							Dissipation = diacoeff * K;
// 							diacoeff    = 0.0;
// 							Dissipation = vol * Epsold;
							if(zone->Flowsol->flagTurb == 4)
								*source = *source - *(zone->Flowsol->Tcoeff->PhiK+indexP);
						}
						if((base->BaseParameters->TimeStepingSchemeTag)==0)
						{
							*source=*source-(Dissipation+Fsum1);
						}
						else if((base->BaseParameters->TimeStepingSchemeTag)==1)
						{
							*source=*source-(Dissipation+Fsum1/2.0);
						}
					}
					diagonalCoef = diacoeff;
				}
			}
		}
		break;
		case 601:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[1]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;

			if(base->flagTurb == 1 && base->BaseParameters->Turb->numTurbScalars >= 2)
			{
				
				if(base->flagTurbModel==2)//k-Omega turbulence model
				{
					*source = *source + gamma * rho * TProduction;
					diacoeff    = vol * beta1 * OmegaOld;
					Dissipation = diacoeff * (2 * Omega - OmegaOld);
					diacoeff = diacoeff * 2;
					
					*source = *source - Dissipation;
				}
				else if(base->flagTurbModel==6)//SST turbulence model
				{
					gamma = F1 * zone->Flowsol->Tcoeff->Gamma1 + (1.0-F1) * zone->Flowsol->Tcoeff->Gamma2;
					beta1 = F1 * zone->Flowsol->Tcoeff->Beta1  + (1.0-F1) * zone->Flowsol->Tcoeff->Beta2;
					CrossDiffusion = *(zone->Flowsol->CrossDiff+indexP);
					
					*source = *source + gamma * rho * TProduction;
					*source = *source + CrossDiffusion;
					
					diacoeff    = vol * beta1 * OmegaOld;
					Dissipation = diacoeff * (2 * Omega - OmegaOld);
					diacoeff = diacoeff * 2;
					
					*source = *source - Dissipation;
				}
				else//k-Epsilon turbulence model (including standard, RNG, Launder-Sharma and Lam-Bremhorst models)
				{
					C0epsi = *(zone->Flowsol->Tcoeff->C0epsi+indexP);
					C1epsi = *(zone->Flowsol->Tcoeff->C1epsi+indexP);
					C2epsi = *(zone->Flowsol->Tcoeff->C2epsi+indexP);
					if(base->flagLogKEps == 1)
					{
						sumGradSquare = CalcGradSquare(base,zone,index,variableTag);
						MuEffEps = *(zone->Flowsol->Scalcoeff[1]+indexP);
	
						*source = *source + vol * (MuEffEps/rho) *sumGradSquare;
						*source = *source + C1epsi * mue_t * TProduction * exp(-K);
						Dissipation = vol * C2epsi * exp(Epsold-K) * (1.0-Epsold);
						diacoeff = vol * C2epsi * exp(Epsold-K);
						Fsum1 = Fsum1 + diacoeff * Eps;
					}
					else
					{
						*source = *source + mue_t * C1epsi * Gamma * TProduction;
// 						*source = *source + C1epsi * 0.09 * Kold * TProduction;//This boosts up the speed of convergence/divergence

						diacoeff    = vol * C2epsi * Gamma;
						Dissipation = diacoeff * (2 * Eps - Epsold);
						diacoeff    = diacoeff * 2;
// 						diacoeff    = vol * C2epsi * Gamma;
// 						Dissipation = diacoeff * Eps;

						if((zone->Flowsol->flagRNG_KEPS == 1) && (zone->Flowsol->flagTurb == 1))
						{
							Dissipation = Dissipation + vol * C0epsi * Gamma * Eps;
							diacoeff    = diacoeff    + vol * C0epsi * Gamma;
						}
						else if(zone->Flowsol->flagTurb == 4)
						{
							*source = *source + *(zone->Flowsol->Tcoeff->PhiEps+indexP);
						}
					}
					if((base->BaseParameters->TimeStepingSchemeTag)==0)
					{
						*source=*source-(Dissipation+Fsum1);
					}
					else if((base->BaseParameters->TimeStepingSchemeTag)==1)
					{
						*source=*source-(Dissipation+(Fsum1)/2.0);
					}
				}
				diagonalCoef=diacoeff;
			}
		}
		break;
		case 602:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[2]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
		}
		break;
		case 603:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[3]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;
		case 604:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[4]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;
		case 605:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[5]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;
		case 606:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[6]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;
		case 607:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[7]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/
		}
		break;
		case 608:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[8]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/		
		}
		break;
		case 609:
		{
			massFlux = zone->Flowsol->Massflux+indexP*6;
			velP = *(zone->PrevFlowsol->Scal[9]+indexP);
			sumMassFlux=0.0;
			for(face=0;face<6;face++)
			{
				sumMassFlux = sumMassFlux + (*(massFlux+face));
			}
			*source = *source + velP*sumMassFlux;
			/*An explicit function can be written in terms of constants and other scalars*/		
		}
		break;
	}
	return (diagonalCoef);
}


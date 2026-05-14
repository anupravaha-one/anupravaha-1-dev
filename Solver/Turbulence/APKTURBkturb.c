#include "APKTURBkturb.h"
/**
The fucntion kturb() performs Gauss-Siedel iterations of the k-epsilon equations for a given Zone.
*/
double kturb(Base *base,Zone *zone,double *variable,double *prevTstepVar,int vartag,double *prevIterVar)
{
	double *Fc,*Fc_n,Fdiff[6],Fsum1,Fsum2;
	double rho,vol,vel,velP,AX,RS,diacoeff=0;
	int NI,NJ,NK,index[3],INDEX,numBC,flagSpheat,flagrho,scheme;
	int k,j,i,face,*ptrSize;
	double residue,res=0.0;
	double *deferCorrFlux,*upwindFlux,*diagConvect;
	double condia,diffdia,Dissipation;
	double alpha,delta_t,omega=0.01;
	double *slabVariable[6],*slabVolume[6];
	double *vols,*beta,convectivecoeff;
	double varCoefficient,varCoefficientPrev;
	double source,DiagSource,CoeffN[6];
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];
	double sumGradSquare,MuEffK,MuEffEps,K,Eps,Kold,Epsold,mue_t,TProduction;
	double C1epsi,C2epsi,Gamma,C0epsi;//,*f1;
	//int printGamma=0;
	FILE *fp;

	//printf("I am in kturb.c\n");
	fp=fopen("Settings/TurbRelaxation.dat","r");
	fscanf(fp,"%lf", &omega);
	fclose(fp);

	/*fp = fopen("Settings/printGamma.dat","r");
	fscanf(fp,"%d", &printGamma);
	fclose(fp);*/
	
	ptrSize= zone->Size;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);

	alpha = base->BaseParameters->RatioCenUP;
	delta_t = base->BaseParameters->Turb->DualTimeStep;
	scheme=base->SolverScheme;

	deferCorrFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	diagConvect=(double *)malloc((NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	upwindFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	
	numBC= zone->Number_of_BC_Patches;
	vols= zone->ZGeom->Volume;
	beta= zone->ZGeom->Beta;
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	flagrho = zone->Flowsol->flagRho;
	//f1 = zone->Flowsol->Tcoeff->f1;
	
	MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of variable*/
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	
	if(base->BaseParameters->flag2d==1)
	{
		
		diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
		MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;

			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

			if(base->BaseParameters->ConvectionScheme==1)
			{
			convfluxQ2D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
			}
			else if(base->BaseParameters->ConvectionScheme==0)
			{
			convfluxC2D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
			}

			convectivecoeff = 1.0;

			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			vol = *(zone->ZGeom->Volume+INDEX);
			vel=*(variable+INDEX);
			velP=*(prevTstepVar+INDEX);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);

			diffdia =(-1.0/rho)*Diffusion2dWithVaryingDiffCoeff(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,Fdiff,CoeffN);

			Fc= (upwindFlux+INDEX*6);
			Fc_n= (deferCorrFlux+INDEX*6);

			condia = (convectivecoeff)*(*(diagConvect+INDEX));
			Fsum1=0.0;
			Fsum2=0.0;
			for(face=0;face<4;face++)
			{
				Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+((convectivecoeff)*(*(Fc+face)));
				Fsum2=Fsum2+(convectivecoeff)*(*(Fc_n+face));
			}

			varCoefficient = CalcvarCoefficient(base,zone,vartag,INDEX);
			varCoefficientPrev = CalcvarCoefficientPrev(base,zone,vartag,INDEX);

			DiagSource = CalcSource(&source,vartag,zone,base,index);

			Kold   = *(zone->PrevFlowsol->Scal[0]+INDEX);
			Epsold = *(zone->PrevFlowsol->Scal[1]+INDEX);
			K   = *(zone->Flowsol->Scal[0]+INDEX);
			Eps = *(zone->Flowsol->Scal[1]+INDEX);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+INDEX);
			Gamma = (Epsold/Kold);
			//Gamma = *(zone->Flowsol->Gamma+INDEX);	//(Epsold/Kold) is calculated outside the dual time step and stored in Gamma for use through the dual time step

			if(Epsold < 0 || Kold < 0)
			{
				//if(printGamma == 1)
					//if(vartag==600)
					//	printf("[%s %d %d %d %d] Gamma = %e, EpsOld = %e, kOld= %e\n", zone->zonename,k,j,i,vartag,Gamma,Epsold,Kold);
				if(Gamma < 0)
					Gamma = 0;	//This is to stabilize the equations by ensuring the positivity of coeffients of Reaction terms (Production, Dissipation etc.)
			}

			if(vartag==600)
			{
				if(base->flagLogKEps == 1)
				{
					sumGradSquare = CalcGradSquare(base,zone,index,vartag);
					MuEffK = *(zone->Flowsol->Scalcoeff[0]+INDEX);

					source = source + vol * (MuEffK/rho) * sumGradSquare;
					source = source + mue_t * TProduction * exp(-Kold) * (1.0+Kold);
					Dissipation = vol * exp(Epsold - Kold) * (1.0+Kold);
					diacoeff = - (vol * delta_t * exp(Epsold - Kold)) + delta_t * mue_t * TProduction * exp(-Kold);
					Fsum1=Fsum1 - (vol * exp(Epsold - Kold))*K + mue_t * TProduction * exp(-Kold) * K;
				}
				else
				{
					source = source + mue_t * TProduction;
					if(zone->Flowsol->flagTurb == 4)
						source = source - *(zone->Flowsol->Tcoeff->PhiK+INDEX);
					Dissipation = vol * K * Gamma;
					diacoeff = vol * delta_t * Gamma;
					//Dissipation = vol * K * (Epsold/Kold);
					//diacoeff = vol * delta_t * (Epsold/Kold) ;
				}
			}
			else if(vartag==601)
			{
				C0epsi = *(zone->Flowsol->Tcoeff->C0epsi+INDEX);
				C1epsi = *(zone->Flowsol->Tcoeff->C1epsi+INDEX);
				C2epsi = *(zone->Flowsol->Tcoeff->C2epsi+INDEX);
				if(base->flagLogKEps == 1)
				{
					sumGradSquare = CalcGradSquare(base,zone,index,vartag);
					MuEffEps = *(zone->Flowsol->Scalcoeff[1]+INDEX);

					source = source + vol * (MuEffEps/rho) *sumGradSquare;
					source = source + C1epsi * mue_t * TProduction * exp(-K);
					Dissipation = vol * C2epsi * exp(Epsold-K) * (1.0-Epsold);
					diacoeff = vol * C2epsi * delta_t * exp(Epsold-K);
					Fsum1=Fsum1+(vol * C2epsi * exp(Epsold-K) * Eps);
				}
				else
				{
					source = source + mue_t * C1epsi * Gamma * TProduction;
					Dissipation = vol * C2epsi * Eps * Gamma ;
					diacoeff =vol * C2epsi * delta_t * Gamma ;
					if((zone->Flowsol->flagRNG_KEPS == 1) && (zone->Flowsol->flagTurb == 1))
					{
						diacoeff = diacoeff + vol * C0epsi * delta_t * Gamma;
						Dissipation = Dissipation +  vol * C0epsi * Gamma* Eps;
					}
					else if(zone->Flowsol->flagTurb == 4)
						source = source + *(zone->Flowsol->Tcoeff->PhiEps+INDEX);

					//Dissipation = vol * C2epsi * (2*Eps - Epsold) * Gamma ;
					//diacoeff =2* vol * C2epsi * delta_t * Gamma ;
				}
			}

			switch(base->BaseParameters->TimeStepingSchemeTag)
			{

				case 0:
					AX = varCoefficient*vol*vel + delta_t*Fsum1 + delta_t*Dissipation ;

					RS = varCoefficientPrev*vol*velP + delta_t*source - delta_t*Fsum2;
				break;
				case 1:
					//For C-N scheme Dissipation term has to be modified.........
					//Currently not implemented for turbulence case ..........
					AX = varCoefficient*vol*vel + delta_t*(Fsum1/2.0) + delta_t*Dissipation;
					RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*(Fsum2/2.0);
				break;
				default:
					printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
					printf("exiting in kturb\n");
					exit(1);
				break;
			}

			residue=RS-AX;
			diacoeff = diacoeff + varCoefficient*vol + delta_t*(diffdia+condia+DiagSource);
			
			*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
			res=res+(residue*residue);
		}
	}
	else /*3D case with varying Diff Coeff*/
	{
		diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
		MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;

			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

			if(base->BaseParameters->ConvectionScheme==1)
			{
			convfluxQ3D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
			}
			else if(base->BaseParameters->ConvectionScheme==0)
			{
			convfluxC3D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
			}

			convectivecoeff = 1.0;

			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			vol = *(zone->ZGeom->Volume+INDEX);
			vel=*(variable+INDEX);
			velP=*(prevTstepVar+INDEX);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);

			diffdia =(-1.0/rho)*Diffusion3dWithVaryingDiffCoeff(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,Fdiff);

			Fc= (upwindFlux+INDEX*6);
			Fc_n= (deferCorrFlux+INDEX*6);

			condia = (convectivecoeff)*(*(diagConvect+INDEX));

			Fsum1=0.0;
			Fsum2=0.0;
			for(face=0;face<6;face++)
			{
				Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+((convectivecoeff)*(*(Fc+face)));
				Fsum2=Fsum2+(convectivecoeff)*(*(Fc_n+face));
			}

			varCoefficient = CalcvarCoefficient(base,zone,vartag,INDEX);
			varCoefficientPrev = CalcvarCoefficientPrev(base,zone,vartag,INDEX);

			DiagSource = CalcSource(&source,vartag,zone,base,index);

			Kold   = *(zone->PrevFlowsol->Scal[0]+INDEX);
			Epsold = *(zone->PrevFlowsol->Scal[1]+INDEX);
			K   = *(zone->Flowsol->Scal[0]+INDEX);
			Eps = *(zone->Flowsol->Scal[1]+INDEX);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+INDEX);
			Gamma = (Epsold/Kold);
			//Gamma = *(zone->Flowsol->Gamma+INDEX);	//(Epsold/Kold) is calculated outside the dual time step and stored in Gamma for use through the dual time step

			if(Epsold < 0 || Kold < 0)
			{
				//if(printGamma == 1)
					//if(vartag==600)
					//	printf("[%s %d %d %d %d] Gamma = %e, EpsOld = %e, kOld= %e\n", zone->zonename,k,j,i,vartag,Gamma,Epsold,Kold);
				if(Gamma < 0)
					Gamma = 0;	//This is to stabilize the equations by ensuring the positivity of coeffients of Reaction terms (Production, Dissipation etc.)
			}

			if(vartag==600)
			{
				if(base->flagLogKEps == 1)
				{
					sumGradSquare = CalcGradSquare(base,zone,index,vartag);
					MuEffK = *(zone->Flowsol->Scalcoeff[0]+INDEX);

					source = source + vol * (MuEffK/rho) * sumGradSquare;
					source = source + mue_t * TProduction * exp(-Kold) * (1.0+Kold);
					Dissipation = vol * exp(Epsold - Kold) * (1.0+Kold);
					diacoeff = - (vol * delta_t * exp(Epsold - Kold)) + delta_t * mue_t * TProduction * exp(-Kold);
					Fsum1=Fsum1 - (vol * exp(Epsold - Kold))*K + mue_t * TProduction * exp(-Kold) * K;
				}
				else
				{
					source = source + mue_t * TProduction;
					if(zone->Flowsol->flagTurb == 4)
						source = source - *(zone->Flowsol->Tcoeff->PhiK+INDEX);
					Dissipation = vol * K * Gamma;
					diacoeff = vol * delta_t * Gamma;
					//Dissipation = vol * K * (Epsold/Kold);
					//diacoeff = vol * delta_t * (Epsold/Kold) ;
				}
			}
			else if(vartag==601)
			{
				C0epsi = *(zone->Flowsol->Tcoeff->C0epsi+INDEX);
				C1epsi = *(zone->Flowsol->Tcoeff->C1epsi+INDEX);
				C2epsi = *(zone->Flowsol->Tcoeff->C2epsi+INDEX);
				if(base->flagLogKEps == 1)
				{
					sumGradSquare = CalcGradSquare(base,zone,index,vartag);
					MuEffEps = *(zone->Flowsol->Scalcoeff[1]+INDEX);

					source = source + vol * (MuEffEps/rho) *sumGradSquare;
					source = source + C1epsi * mue_t * TProduction * exp(-K);
					Dissipation = vol * C2epsi * exp(Epsold-K) * (1.0-Epsold);
					diacoeff = vol * C2epsi * delta_t * exp(Epsold-K);
					Fsum1=Fsum1+(vol * C2epsi * exp(Epsold-K) * Eps);
				}
				else
				{
					source = source + mue_t * C1epsi * Gamma * TProduction;
					Dissipation = vol * C2epsi * Eps * Gamma ;
					diacoeff =vol * C2epsi * delta_t * Gamma ;
					if((zone->Flowsol->flagRNG_KEPS == 1) && (zone->Flowsol->flagTurb == 1))
					{
						diacoeff = diacoeff + vol * C0epsi * delta_t * Gamma;
						Dissipation = Dissipation +  vol * C0epsi * Gamma* Eps;
					}
					else if(zone->Flowsol->flagTurb == 4)
						source = source + *(zone->Flowsol->Tcoeff->PhiEps+INDEX);

					//Dissipation = vol * C2epsi * (2*Eps - Epsold) * Gamma ;
					//diacoeff =2* vol * C2epsi * delta_t * Gamma ;
				}
			}

			switch(base->BaseParameters->TimeStepingSchemeTag)
			{
				case 0:
					AX = varCoefficient*vol*vel + delta_t*Fsum1 + delta_t*Dissipation;
					RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*Fsum2;
				break;
				case 1:
					//For C-N scheme Dissipation term has to be modified.........
					AX = varCoefficient*vol*vel + delta_t*(Fsum1/2.0) + delta_t*Dissipation;
					RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*(Fsum2/2.0);
				break;
				default:
					printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
					printf("exiting in UpdateVar\n");
					exit(1);
				break;
			}

			residue=RS-AX;
			diacoeff = diacoeff + varCoefficient*vol + delta_t*(diffdia+condia+DiagSource);
			*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
			res=res+(residue*residue);
		}
	}
	
	free(deferCorrFlux);
	free(upwindFlux);
	free(diagConvect);

	return(res);
}




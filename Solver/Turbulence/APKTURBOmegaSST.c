#include "APKTURBOmegaSST.h"
/**
The fucntion Omega() performs Gauss-Siedel iterations of the k-Omega equations for a given Zone.
*/

double OmegaSST(Base *base,Zone *zone,double *variable,double *prevTstepVar,int vartag,double *prevIterVar)
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
	double *vols,*Beta,convectivecoeff;
	double varCoefficient,varCoefficientPrev;
	double source=0,DiagSource,CoeffN[6];
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];
	double K,Omega,Kold,OmegaOld,mue_t,TProduction,CrossDiffusion=0;
	double gamma,beta,betastar=0.09,F1;
	int flagViscoeff;
	FILE *fp;

	fp=fopen("Settings/OmegaRelaxation.dat","r");
	fscanf(fp,"%lf", &omega);
	fclose(fp);
	omega=0.7;

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
	Beta= zone->ZGeom->Beta;
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	flagrho = zone->Flowsol->flagRho;
	flagViscoeff = zone->Flowsol->flagViscoeff;

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

			diffdia =(-1.0/rho)*Diffusion2dWithVaryingDiffCoeff(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,Beta,Fdiff,CoeffN);

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

			Kold     = *(zone->PrevFlowsol->Scal[0]+INDEX);
			OmegaOld = *(zone->PrevFlowsol->Scal[1]+INDEX);
			K        = *(zone->Flowsol->Scal[0]+INDEX);
			Omega    = *(zone->Flowsol->Scal[1]+INDEX);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+INDEX);
			betastar = (zone->Flowsol->Tcoeff->betastar);
			F1	 = *(zone->Flowsol->Tcoeff->F1+INDEX);

			if(vartag==600)
			{
				betastar = zone->Flowsol->Tcoeff->betastar;

				source = source + mue_t * TProduction;
				Dissipation = vol * betastar * OmegaOld * K;
				diacoeff    = vol * betastar * OmegaOld * delta_t;
			}
			else if(vartag==601)
			{
				gamma = F1 * zone->Flowsol->Tcoeff->Gamma1 + (1.0-F1) * zone->Flowsol->Tcoeff->Gamma2;
				beta  = F1 * zone->Flowsol->Tcoeff->Beta1  + (1.0-F1) * zone->Flowsol->Tcoeff->Beta2;
				CrossDiffusion = *(zone->Flowsol->CrossDiff+INDEX);

				source = source + gamma * rho * TProduction;
				source = source + CrossDiffusion;

				Dissipation = vol * beta * OmegaOld * (2 * Omega - OmegaOld);
				diacoeff    = vol * beta * OmegaOld * (2 * delta_t);
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
					printf("C-N scheme is currently not implemented for turbulence problems\n");
					printf("*** Choose BTCS scheme ... exiting the program *** \n");
					exit(10);
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

			diffdia =(-1.0/rho)*Diffusion3dWithVaryingDiffCoeff(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,Beta,Fdiff);

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

			Kold     = *(zone->PrevFlowsol->Scal[0]+INDEX);
			OmegaOld = *(zone->PrevFlowsol->Scal[1]+INDEX);
			K        = *(zone->Flowsol->Scal[0]+INDEX);
			Omega    = *(zone->Flowsol->Scal[1]+INDEX);
			TProduction = *(zone->Flowsol->Tcoeff->TProduction+INDEX);
			betastar = (zone->Flowsol->Tcoeff->betastar);
			F1	 = *(zone->Flowsol->Tcoeff->F1+INDEX);

			if(vartag==600)
			{
				betastar = zone->Flowsol->Tcoeff->betastar;

				source = source + mue_t * TProduction;
				Dissipation = vol * betastar * OmegaOld * K;
				diacoeff    = vol * betastar * OmegaOld * delta_t;
			}
			else if(vartag==601)
			{
				gamma = F1 * zone->Flowsol->Tcoeff->Gamma1 + (1.0-F1) * zone->Flowsol->Tcoeff->Gamma2;
				beta  = F1 * zone->Flowsol->Tcoeff->Beta1  + (1.0-F1) * zone->Flowsol->Tcoeff->Beta2;
				CrossDiffusion = *(zone->Flowsol->CrossDiff+INDEX);

				source = source + gamma * rho * TProduction;
				source = source + CrossDiffusion;

				Dissipation = vol * beta * OmegaOld * (2 * Omega - OmegaOld);
				diacoeff    = vol * beta * OmegaOld * (2 * delta_t);
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
					printf("C-N scheme is currently not implemented for turbulence problems\n");
					printf("*** Choose BTCS scheme ... exiting the program *** \n");
					exit(10);
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

	free(deferCorrFlux);
	free(upwindFlux);
	free(diagConvect);

	return(res);
}




#include "APKCNterm.h"
/**
The function UpdateCNterm() computes and stores the convective and diffusive fluxes at current time level to
be used in the next timestep, when Crank Nicolson time stepping scheme is selected. The function is
split up into two for 2D and 3D domains and further into two cases for each depending on whether the Diffusion coefficient
is a constant or a variable.
Note: This function takes into consideration the multiplication factor "1/2" and also the necessary
signs for the fluxes. It sums up both the convective and diffusive fluxes.
*/

/*********************************************************************************************************
**********************************************************************************************************
########################## Call this before fillNewtoPreviousFlowsol #####################################
**********************************************************************************************************
*********************************************************************************************************/
int UpdateCNterm(Base *base,Zone *zone,int vartag)
{
	double *domainVariable;
	double *prevIterVar;
	double *Fc,*Fc_n,Fdiff[6],Fsum1,Fsum2;
	int NI,NJ,NK,index[3],indexP,scheme;
	int k,j,i,face,*ptrSize,flagdiffusionCoeff,flagrho,flagSpheat;
	double *deferCorrFlux,*upwindFlux,*diagConvect,convectivecoeff;
	double diffdia,rho;
	double alpha,delta_t;
	double *slabVariable[6],*slabVolume[6];
	double *vols,*beta;
	double coeffn[6];
	//double source,DiagSource;
	double *CNTerm;
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];

	ptrSize= zone->Size;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);

	alpha = base->BaseParameters->RatioCenUP;
	delta_t = base->BaseParameters->TimeStep;
	prevIterVar  =(double *)malloc(  (NK-1)*(NJ-1)*(NI-1)*sizeof(double));

	deferCorrFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	diagConvect  =(double *)malloc(  (NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	upwindFlux   =(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));

	domainVariable=MakeDomainArray(zone,vartag);//prevTstepVar stores prevFlowsol->variable 
	
	vols = zone->ZGeom->Volume;
	beta = zone->ZGeom->Beta;
	flagrho = zone->Flowsol->flagRho;
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	
	MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of u */
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/

	fillnewtoprev(domainVariable,prevIterVar,ptrSize);
	scheme=base->SolverScheme;
	
	CNTerm = MakeDomainArrayCNTerm(zone,vartag);
	if(base->BaseParameters->flag2d==1)
	{
		flagdiffusionCoeff = MakeDiffusionCoeffFlag(zone,vartag);
		if(flagdiffusionCoeff==0) /*2D case with Const Diff coeff*/
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index[0]=k,index[1]=j,index[2]=i;

				indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

				if(base->BaseParameters->ConvectionScheme==1)
					convfluxQ2D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				else if(base->BaseParameters->ConvectionScheme==0)
					convfluxC2D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
					
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;

				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				diffdia =(-1.0/rho)*Diffusion2dWithConstDiffCoeff(index,domainVariable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,beta,Fdiff,coeffn);

				Fc= (upwindFlux+indexP*6);
				Fc_n= (deferCorrFlux+indexP*6);

				Fsum1=0.0;
				Fsum2=0.0;

				for(face=0;face<4;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+(convectivecoeff)*(*(Fc+face));
					Fsum2=Fsum2+(convectivecoeff)*(*(Fc_n+face));
				}

				*(CNTerm+indexP)=(-1.0)*(Fsum2+Fsum1)/2.0;
				/*Calculated as a part of the Right Hand Side*/
			}
		}
		else
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
			MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index[0]=k,index[1]=j,index[2]=i;

				indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

				if(base->BaseParameters->ConvectionScheme==1)
					convfluxQ2D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				else if(base->BaseParameters->ConvectionScheme==0)
					convfluxC2D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
					
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;

				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				diffdia =(-1.0/rho)*Diffusion2dWithVaryingDiffCoeff(index,domainVariable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,Fdiff,coeffn);

				Fc= (upwindFlux+indexP*6);
				Fc_n= (deferCorrFlux+indexP*6);

				Fsum1=0.0;
				Fsum2=0.0;

				for(face=0;face<4;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+(convectivecoeff)*(*(Fc+face));
					Fsum2=Fsum2+(convectivecoeff)*(*(Fc_n+face));
				}

				*(CNTerm+indexP)=(-1.0)*(Fsum2+Fsum1)/2.0;
				/*Calculated as a part of the Right Hand Side*/
			}
		}
		
	}
	else
	{
		flagdiffusionCoeff = MakeDiffusionCoeffFlag(zone,vartag);
		if(flagdiffusionCoeff==0) /*3D case with Const Diff coeff*/
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index[0]=k,index[1]=j,index[2]=i;

				indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

				if(base->BaseParameters->ConvectionScheme==1)
					convfluxQ3D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				else if(base->BaseParameters->ConvectionScheme==0)
					convfluxC3D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				
					
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;

				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				diffdia =(-1.0/rho)*Diffusion3dWithConstDiffCoeff(index,domainVariable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,beta,Fdiff);

				Fc= (upwindFlux+indexP*6);
				Fc_n= (deferCorrFlux+indexP*6);

				Fsum1=0.0;
				Fsum2=0.0;

				for(face=0;face<6;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+(convectivecoeff)*(*(Fc+face));
					Fsum2=Fsum2+(convectivecoeff)*(*(Fc_n+face));
				}

				*(CNTerm+indexP)=(-1.0)*(Fsum2+Fsum1)/2.0;
				/*Calculated as a part of the Right Hand Side*/

			}
		}
		else
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
			MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index[0]=k,index[1]=j,index[2]=i;

				indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

				if(base->BaseParameters->ConvectionScheme==1)
					convfluxQ3D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				else if(base->BaseParameters->ConvectionScheme==0)
					convfluxC3D(prevIterVar,domainVariable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				
					
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*indexP);
				else
					convectivecoeff = 1.0;

				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				diffdia =(-1.0/rho)*Diffusion3dWithVaryingDiffCoeff(index,domainVariable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,Fdiff);

				Fc= (upwindFlux+indexP*6);
				Fc_n= (deferCorrFlux+indexP*6);

				Fsum1=0.0;
				Fsum2=0.0;

				for(face=0;face<6;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+(convectivecoeff)*(*(Fc+face));
					Fsum2=Fsum2+(convectivecoeff)*(*(Fc_n+face));
				}

				*(CNTerm+indexP)=(-1.0)*(Fsum2+Fsum1)/2.0;
				/*Calculated as a part of the Right Hand Side*/

			}
		}
	}
	
	free(deferCorrFlux);
	free(upwindFlux);
	free(diagConvect);
	free(prevIterVar);
	
	return 0;
}

int CalcPrevTStepFluxes(Base *base)
{
	Zone *zone;
	int cntZone,numZones,cntVar;
	numZones=base->Number_of_zones;
	for(cntZone=0;cntZone<numZones;cntZone++)
	{
		zone=base->CgnsZone+cntZone;
		

		if(zone->Flowsol->u!=NULL)
			UpdateCNterm(base,zone,1);

		if(zone->Flowsol->v!=NULL)
			UpdateCNterm(base,zone,2);

		if(zone->Flowsol->w!=NULL)
			UpdateCNterm(base,zone,3);
		
		if(zone->Flowsol->T!=NULL)
			UpdateCNterm(base,zone,6);
		
		for(cntVar=0;cntVar<10;cntVar++)
			if(zone->Flowsol->Scal[cntVar]!=NULL)
				UpdateCNterm(base,zone,600+cntVar);
	/*PrintCNFluxes(zone);*/
	}
	//printf("wrote down the CN fluxes\n");
	return 0;
}

double *MakeDomainArrayCNTerm(Zone *zone,int vartag)
{
	double *domain;
	
	switch(vartag)
	{
		case 1:/**/
			domain=zone->Flowsol->PrevTStepFluxes_u;
		break;
		case 2:/**/
			domain=zone->Flowsol->PrevTStepFluxes_v;
		break;
		case 3:/**/
			domain=zone->Flowsol->PrevTStepFluxes_w;
		break;
		case 6:
			domain=zone->Flowsol->PrevTStepFluxes_T;
		break;
		
		case 600:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[0];
		break;
		case 601:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[1];
		break;
		case 602:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[2];
		break;
		case 603:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[3];
		break;
		case 604:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[4];
		break;
		case 605:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[5];
		break;
		case 606:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[6];
		break;
		case 607:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[7];
		break;
		case 608:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[8];
		break;
		case 609:
			domain=zone->Flowsol->PrevTStepFluxes_Scal[9];
		break;
	}
	return (domain);
}


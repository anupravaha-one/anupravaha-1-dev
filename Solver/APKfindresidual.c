#include "APKfindresidual.h"

/**
This function solves for the mass velocities only using Gauss-Seidal point method and returns the
residual of the mass velocities for that particular iteration. This is implemented for both 2D &
3D cases having either constant or variable diffusion coefficient.
*/

double findresidual(Base *base,Zone *zone,double *massVelDomain,double *prevTstepVel,int *ptrsize,int vartag,double *prevIterMassVel)
{
	double *Fc,*Fc_n,Fdiff[6],Fsum1,Fsum2;
	double rho,vol,vel,velP,AX,RS,diacoeff;
	int NI,NJ,NK,index[3],INDEX,numBC,scheme;
	int k,j,i,face,flagrho,flagdiffusionCoeff;
	double residue,res=0.0;
	double *deferCorrFlux,*upwindFlux,*diagConvect;
	double condia,diffdia,DiagSource,source;
	double alpha,delta_t,omega=0.6;
	double *slabVelocity[6],*slabVolume[6];
	double *vols,*beta,CoeffN[6];
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];
	
	//ZoneBoundary *zonebc;
	//int nBc;
	FILE *Velomega;
	
	Velomega = fopen("Settings/VelocityRelaxation.dat","r");
	fscanf(Velomega,"%lf", &omega);
	fclose(Velomega);
	
	NI=*(ptrsize+0);
	NJ=*(ptrsize+1);
	NK=*(ptrsize+2);
	
	alpha = base->BaseParameters->RatioCenUP;
	delta_t = base->BaseParameters->TimeStep;
	scheme=base->SolverScheme;
	
	deferCorrFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	diagConvect=(double *)malloc((NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	upwindFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
		
	numBC= zone->Number_of_BC_Patches;
	
	vols= zone->ZGeom->Volume;
	beta= zone->ZGeom->Beta;
	flagrho = zone->Flowsol->flagRho;
		
	MakeSlabPtrArray(zone,slabVelocity,vartag*100);/*Called for making slabs of u */
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
		
	if(base->BaseParameters->flag2d==1)
	{
		flagdiffusionCoeff = MakeDiffusionCoeffFlag(zone,vartag*100);
		
		if(flagdiffusionCoeff==0 && base->flagTurb == 0) /*2D case with Const Diff coeff*/
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag*100);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				
				if(base->BaseParameters->ConvectionScheme==1)
				{
				convfluxQ2D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				else if(base->BaseParameters->ConvectionScheme==0)
				{
				convfluxC2D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}

				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel= *(massVelDomain+INDEX);
				velP = *(prevTstepVel+INDEX);

				diffdia =(-1.0/rho)*Diffusion2dWithConstDiffCoeff(index,massVelDomain,vols,ptrsize,slabVelocity,slabVolume,diffusionCoeffdomain,beta,Fdiff,CoeffN);

				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);

				condia = *(diagConvect+INDEX);

				Fsum1=0.0;
				Fsum2=0.0;
				for(face=0;face<4;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+*(Fc+face);
					Fsum2=Fsum2+*(Fc_n+face);
				}

				DiagSource=CalcSource(&source,(vartag*100),zone,base,index);

				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX= vol*vel + delta_t*Fsum1;
						RS=delta_t*source + vol*velP - delta_t*Fsum2;
					break;
					case 1:
						AX= vol*vel + delta_t*(Fsum1/2.0);
						RS= delta_t*source + vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in findresidual\n");
						exit(1);
					break;
				}

				residue=RS-AX;

				diacoeff= vol + delta_t*(diffdia+condia+DiagSource);

				*(massVelDomain+INDEX)= *(massVelDomain+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
		else /*2D case with varying Diff coeff*/
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag*100);
			MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag*100);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				if(base->BaseParameters->ConvectionScheme==1)
				{
				convfluxQ2D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				else if(base->BaseParameters->ConvectionScheme==0)
				{
				convfluxC2D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				
				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(massVelDomain+INDEX);
				velP=*(prevTstepVel+INDEX);
				
				diffdia =(-1.0/rho)*Diffusion2dWithVaryingDiffCoeff(index,massVelDomain,vols,ptrsize,slabVelocity,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,Fdiff,CoeffN);
				
				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);
				
				condia = *(diagConvect+INDEX);
				
				Fsum1=0.0;
				Fsum2=0.0;
				for(face=0;face<4;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+*(Fc+face);
					Fsum2=Fsum2+*(Fc_n+face);
					//printf("%d %d %d Fdiff[%d]=%e\tFconUP2=%e\tFconvdiff=%lf\n",k,j,i,face,Fdiff[face],*(Fc+face),*(Fc_n+face));
				}
				
				DiagSource=CalcSource(&source,(vartag*100),zone,base,index);
				//printf("%d %d %d Fsum1=%e\tFsum2=%e\tDiagSource=%lf\n",k,j,i,Fsum1,Fsum2,DiagSource);
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX= vol*vel + delta_t*Fsum1;
						RS=delta_t*source + vol*velP - delta_t*Fsum2;
					break;
					case 1:
						AX= vol*vel + delta_t*(Fsum1/2.0);
						RS= delta_t*source + vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in findresidual\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				//printf("source value=%lf\tDiagSource = %lf\n",source,DiagSource);
				//printf("residue = %lf\n",residue);
				diacoeff= vol + delta_t*(diffdia+condia+DiagSource);
				//printf("%d %d %d source value=%lf\tDiagSource = %lf\tdiacoeff=%lf\tresidue=%lf\n",k,j,i,source,DiagSource,diacoeff,residue);
				
				*(massVelDomain+INDEX)= *(massVelDomain+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
	}
	else /*3D case*/
	{
		flagdiffusionCoeff = MakeDiffusionCoeffFlag(zone,vartag*100);
		
		if(flagdiffusionCoeff==0 && base->flagTurb==0) /*3D case with Constant Diff coeff*/
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag*100);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				if(base->BaseParameters->ConvectionScheme==1)
				{
				convfluxQ3D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				else if(base->BaseParameters->ConvectionScheme==0)
				{
				convfluxC3D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				
				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(massVelDomain+INDEX);
				velP=*(prevTstepVel+INDEX);

				diffdia =(-1.0/rho)*Diffusion3dWithConstDiffCoeff(index,massVelDomain,vols,ptrsize,slabVelocity,slabVolume,diffusionCoeffdomain,beta,Fdiff);

				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);

				condia = *(diagConvect+INDEX);

				Fsum1=0.0;
				Fsum2=0.0;
				for(face=0;face<6;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+*(Fc+face);
					Fsum2=Fsum2+*(Fc_n+face);
				}

				DiagSource=CalcSource(&source,(vartag*100),zone,base,index);

				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX= vol*vel + delta_t*Fsum1;
						RS=delta_t*source + vol*velP - delta_t*Fsum2;
					break;
					case 1:
						AX= vol*vel + delta_t*(Fsum1/2.0);
						RS= delta_t*source + vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in findresidual\n");
						exit(1);
					break;
				}

				residue=RS-AX;

				diacoeff= vol + delta_t*(diffdia+condia+DiagSource);

				*(massVelDomain+INDEX)= *(massVelDomain+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
		else /*3D case with varying Diff coeff*/
		{
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag*100);
			MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag*100);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				if(base->BaseParameters->ConvectionScheme==1)
				{
				convfluxQ3D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				else if(base->BaseParameters->ConvectionScheme==0)
				{
				convfluxC3D(prevIterMassVel,massVelDomain,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag*100,scheme,index);
				}
				
				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(massVelDomain+INDEX);
				velP=*(prevTstepVel+INDEX);

				diffdia =(-1.0/rho)*Diffusion3dWithVaryingDiffCoeff(index,massVelDomain,vols,ptrsize,slabVelocity,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,Fdiff);

				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);

				condia = *(diagConvect+INDEX);

				Fsum1=0.0;
				Fsum2=0.0;
				for(face=0;face<6;face++)
				{
					Fsum1=Fsum1-(1.0/rho)*Fdiff[face]+*(Fc+face);
					Fsum2=Fsum2+*(Fc_n+face);
				}

				DiagSource=CalcSource(&source,(vartag*100),zone,base,index);

				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX= vol*vel + delta_t*Fsum1;
						RS=delta_t*source + vol*velP - delta_t*Fsum2;
					break;
					case 1:
						AX= vol*vel + delta_t*(Fsum1/2.0);
						RS= delta_t*source + vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in findresidual\n");
						exit(1);
					break;
				}

				residue=RS-AX;

				diacoeff= vol + delta_t*(diffdia+condia+DiagSource);

				*(massVelDomain+INDEX)= *(massVelDomain+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
	
	}
	
// 	ExtractMassVelocities(vartag*100,base,zone);
// 	for(nBc=0;nBc<numBC;nBc++)
// 	{
// 		zonebc = (zone->ZoneBc+nBc);	
// 		velocitymassBc(base,vartag,zone,zonebc);/* Ensure that the right tag is supplied for */
// 	}
	
	free(deferCorrFlux);
	free(upwindFlux);
	free(diagConvect);
	
	return(res);
}

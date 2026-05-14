#include "APKUpdateVar.h"

/**
This is a general function unlike "findresidual" (which only solves for mass velocities).
This function solves for any variable using Gauss-Seidal point method and returns the
residual of that variable for that particular iteration. This is implemented for both 2D &
3D cases having either constant or variable diffusion coefficient.
*/


double UpdateVar(Base *base,Zone *zone,double *variable,double *prevTstepVar,int vartag,double *prevIterVar)
{
	double *Fc,*Fc_n,Fdiff[6],Fsum1,Fsum2;
	double rho,vol,vel,velP,AX,RS,diacoeff=0.0;
	int NI,NJ,NK,index[3],INDEX,numBC,flagSpheat,flagdiffusionCoeff,flagrho,scheme;
	int k,j,i,face,*ptrSize,flagturbmodel;
	double residue,res=0.0;
	double *deferCorrFlux,*upwindFlux,*diagConvect;
	
	double condia,diffdia;
	double alpha,delta_t,omega=0.6;
	double *slabVariable[6],*slabVolume[6];
	double *vols,*beta,convectivecoeff;
	double varCoefficient,varCoefficientPrev;
	double source=0,DiagSource,CoeffN[6];
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];
	
	FILE *Velomega;
	
	//for Solidification and Melting problems
	double *defer_dele,*upwind_dele,*diag_dele;
	double *dele,Fsum1_dele,Fsum2_dele;

	//for LMMHD problems naveen
	double B_x,B_y,B_z,sigma;
	
	flagturbmodel=(base->flagTurbModel);
	
	alpha = base->BaseParameters->RatioCenUP;
	delta_t = base->BaseParameters->TimeStep;
	scheme=base->SolverScheme;
	
	if(vartag==600 || vartag==601)
		delta_t = base->BaseParameters->Turb->DualTimeStep;

	if(base->flagTurb==1 && (vartag==600 || vartag==601))
	{
		if(vartag==600)
			omega = *(base->BaseParameters->OMEGA+4);
		if(vartag==601)
			omega = *(base->BaseParameters->OMEGA+5);
		
// 		if((flagturbmodel==2)||(flagturbmodel==6))
// 		{
// 			Velomega = fopen("Settings/OmegaRelaxation.dat","r");
// 			fscanf(Velomega,"%lf", &omega);
// 		}
// 		else
// 		{
// 			Velomega = fopen("Settings/TurbRelaxation.dat","r");
// 			fscanf(Velomega,"%lf", &omega);
// 		}
	}
	else
	{
		Velomega = fopen("Settings/VelocityRelaxation.dat","r");
		fscanf(Velomega,"%lf", &omega);	
	}
	fclose(Velomega);

	ptrSize= zone->Size;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);

	
	if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
	{
		defer_dele =(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
		diag_dele  =(double *)malloc((NK-1)*(NJ-1)*(NI-1)*sizeof(double));
		upwind_dele=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	}
	
	deferCorrFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	diagConvect=(double *)malloc((NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	upwindFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	
	
	numBC= zone->Number_of_BC_Patches;
	vols= zone->ZGeom->Volume;
	beta= zone->ZGeom->Beta;
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	flagrho = zone->Flowsol->flagRho;
	
	dele    = zone->Flowsol->Latent_Heat_Content;
	
	MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of u*/
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	
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
				
				INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
				
				if(base->BaseParameters->ConvectionScheme==1)
				{
				convfluxQ2D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				}
				else if(base->BaseParameters->ConvectionScheme==0)
				{
				convfluxC2D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				}
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*INDEX);
				else
					convectivecoeff = 1.0;
				
				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				
				diffdia =(-1.0/rho)*Diffusion2dWithConstDiffCoeff(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,beta,Fdiff,CoeffN);
				
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
				
				DiagSource=CalcSource(&source,vartag,zone,base,index);

				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = varCoefficient*vol*vel + delta_t*Fsum1;
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*Fsum2;
					break;
					case 1:
						AX = varCoefficient*vol*vel + delta_t*(Fsum1/2.0);
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				diacoeff = varCoefficient*vol + delta_t*(diffdia+(convectivecoeff*(condia+DiagSource)));

				if(base->flagLMMHD==1)//naveen
				{
					sigma = *(zone->Flowsol->Eleconductivity+INDEX);
					if(vartag==1)
					{
						B_y = *(zone->Flowsol->MfieldY+INDEX);
						B_z = *(zone->Flowsol->MfieldZ+INDEX);
						diacoeff = diacoeff+vol*delta_t*(sigma/rho) *(B_y*B_y+B_z*B_z) ;
						residue = residue -vol*delta_t* (sigma/rho) *(B_y*B_y+B_z*B_z)*vel ;
					}
					if(vartag==2)
					{
						B_x = *(zone->Flowsol->MfieldX+INDEX);
						B_z = *(zone->Flowsol->MfieldZ+INDEX);
						diacoeff = diacoeff+vol*delta_t*(sigma/rho) *(B_x*B_x+B_z*B_z) ;
						residue = residue - vol*delta_t*(sigma/rho) *(B_x*B_x+B_z*B_z)*vel ;
					}
				}
				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
		else /*2D case with Varying Diff Coeff*/
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
				
				if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
				convfluxQ2D(dele,dele,0.0,zone,defer_dele,upwind_dele,diag_dele,21,scheme,index);
				
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*INDEX);
				else
					convectivecoeff = 1.0;

				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				
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
				if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
				{
					Fsum1_dele=0.0;
					Fsum2_dele=0.0;
					for(face=0;face<4;face++)
					{
						Fsum1_dele = Fsum1_dele + *(upwind_dele+ INDEX *6 +face);
						Fsum2_dele = Fsum2_dele + *(defer_dele + INDEX *6 +face);
					}
				}
				varCoefficient = CalcvarCoefficient(base,zone,vartag,INDEX);
				varCoefficientPrev = CalcvarCoefficientPrev(base,zone,vartag,INDEX);
				
				DiagSource = CalcSource(&source,vartag,zone,base,index);
				
				
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = varCoefficient*vol*vel + delta_t*Fsum1;
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*Fsum2;
						if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
						RS = RS - delta_t * Fsum1_dele;
					break;
					case 1:
						AX = varCoefficient*vol*vel + delta_t*(Fsum1/2.0);
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				
				diacoeff = varCoefficient*vol + delta_t*(diffdia+(convectivecoeff*(condia+DiagSource)));
				
				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				if(INDEX==0 && vartag==601)
				{
				printf("Diagsource = %e\tsource=%e\n",DiagSource,source);
				printf("residues=%e\tRS=%e\tAX=%e\tdiacoeff=%e\tvalue=%e\n",residue,RS,AX,diacoeff,*(variable+INDEX));
				}
				//printf("RS = %e\tAX=%e\tresidue=%e\tvalue=%e\tdiacoeff=%e\n",RS,AX,residue,*(variable+INDEX),diacoeff);
				res=res+(residue*residue);
			}
		}
		
	}
	else /*3D case*/
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
				
				INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
				
				if(base->BaseParameters->ConvectionScheme==1)
				{
				convfluxQ3D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				}
				else if(base->BaseParameters->ConvectionScheme==0)
				{
				convfluxC3D(prevIterVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
				}

				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*INDEX);
				else
					convectivecoeff = 1.0;

				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				
				diffdia =(-1.0/rho)*Diffusion3dWithConstDiffCoeff(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,beta,Fdiff);
				
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
				
				DiagSource=CalcSource(&source,vartag,zone,base,index);
				
				

				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = varCoefficient*vol*vel + delta_t*Fsum1;
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*Fsum2;
					break;
					case 1:
						AX = varCoefficient*vol*vel + delta_t*(Fsum1/2.0);
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				diacoeff = varCoefficient*vol + delta_t*(diffdia+(convectivecoeff*(condia+DiagSource)));
				if(base->flagLMMHD==1)//naveen
				{
					sigma = *(zone->Flowsol->Eleconductivity+INDEX);
					if(vartag==1)
					{
						B_y = *(zone->Flowsol->MfieldY+INDEX);
						B_z = *(zone->Flowsol->MfieldZ+INDEX);
						diacoeff = diacoeff+vol*delta_t*(sigma/rho) *(B_y*B_y+B_z*B_z) ;
						residue = residue -vol*delta_t* (sigma/rho) *(B_y*B_y+B_z*B_z)*vel ;
					}
					if(vartag==2)
					{
						B_x = *(zone->Flowsol->MfieldX+INDEX);
						B_z = *(zone->Flowsol->MfieldZ+INDEX);
						diacoeff = diacoeff+vol*delta_t*(sigma/rho) *(B_x*B_x+B_z*B_z) ;
						residue = residue -vol*delta_t* (sigma/rho) *(B_x*B_x+B_z*B_z)*vel ;
					}
					if(vartag==3)
					{
						B_x = *(zone->Flowsol->MfieldX+INDEX);
						B_y = *(zone->Flowsol->MfieldY+INDEX);
						diacoeff = diacoeff+vol*delta_t*(sigma/rho) *(B_x*B_x+B_y*B_y) ;
						residue = residue -vol*delta_t* (sigma/rho) *(B_x*B_x+B_y*B_y)*vel ;
					}
				}

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
				
				if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
				convfluxQ3D(dele,dele,0.0,zone,defer_dele,upwind_dele,diag_dele,21,scheme,index);
				
				if(vartag==6)
					convectivecoeff = *((zone->Flowsol->SpecificHeat)+flagSpheat*INDEX);
				else
					convectivecoeff = 1.0;
				
				rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				
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
				if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
				{
					Fsum1_dele=0.0;
					Fsum2_dele=0.0;
					for(face=0;face<6;face++)
					{
						Fsum1_dele = Fsum1_dele + *(upwind_dele+ INDEX *6 +face);
						Fsum2_dele = Fsum2_dele + *(defer_dele + INDEX *6 +face);
					}
				}
				varCoefficient = CalcvarCoefficient(base,zone,vartag,INDEX);
				varCoefficientPrev = CalcvarCoefficientPrev(base,zone,vartag,INDEX);

				DiagSource = CalcSource(&source,vartag,zone,base,index);
				
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = varCoefficient*vol*vel + delta_t*Fsum1;
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*Fsum2;
						if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
						RS = RS - delta_t * Fsum1_dele;
					break;
					case 1:
						AX = varCoefficient*vol*vel + delta_t*(Fsum1/2.0);
						RS = delta_t*source + varCoefficientPrev*vol*velP - delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				diacoeff = varCoefficient*vol + delta_t*(diffdia+(convectivecoeff*(condia+DiagSource)));

				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
	}
	free(deferCorrFlux);
	free(upwindFlux);
	free(diagConvect);
	if((base->flagSolidification==1) && (base->BaseParameters->flag_mushy)==1 && (vartag==7))
	{
		free(defer_dele);
		free(upwind_dele);
		free(diag_dele);
	}
	return(res);
}



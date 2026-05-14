#include "APKUpdateScaler.h"

/**
This function solves for any scaler for groundwater using Gauss-Seidal point method and returns the
residual of that variable for that particular iteration. This is implemented for both 2D &
3D cases and for anisotropic cases. It can solve contaminant in both saturated and unsaturated medium. It uses new function Dispersionflux2d and 3d for calculation of dispersion flux.
*/

double UpdateScaler(Base *base,Zone *zone,double *variable,double *prevTstepVar,int vartag,double *prevIterVar)
{
	double *Fc,*Fc_n,Fdiff[6],Fsum1,Fsum2;
	double vol,vel,velP,AX,RS,diacoeff;
	int NI,NJ,NK,index[3],INDEX,numBC,scheme;
	int k,j,i,face,*ptrSize,number;
	double residue,res=0.0;
	double *deferCorrFlux,*upwindFlux,*diagConvect;
	double condia,diffdia;
	double alpha,delta_t,omega=0.6;
	double *slabVariable[6],*slabVolume[6];
	double *vols;
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];
	double Vartag1,Vartag2,Vartag3;
	double  *slabGradX[6],*slabGradY[6],*slabGradZ[6],*GradX,*GradY,*GradZ,*SurfArea,*DeltaInverse,*slabDeltaInverse[6];
	double Porosity;
	int tag;
	double *theta,*theta_Pre;
	FILE *Velomega;

	Velomega = fopen("Settings/VelocityRelaxation.dat","r");
	fscanf(Velomega,"%lf", &omega);
	fclose(Velomega);

	ptrSize= zone->Size;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);

	alpha = base->BaseParameters->RatioCenUP;
	delta_t = base->BaseParameters->TimeStep;
	scheme=base->SolverScheme;

	deferCorrFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	diagConvect=(double *)malloc((NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	upwindFlux=(double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
	
	
	numBC= zone->Number_of_BC_Patches;
	vols= zone->ZGeom->Volume;
	
	SurfArea= zone->ZGeom->SurfaceArea;
	
	
	number=vartag-600;
	Vartag1=60000+number;
	Vartag2=60100+number;
	Vartag3=60200+number;
		
	GradX=zone->Flowsol->GradScalX[number];
	GradY=zone->Flowsol->GradScalY[number];
	GradZ=zone->Flowsol->GradScalZ[number];
	DeltaInverse=zone->ZGeom->DeltaInverse;
	
	MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of u*/
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	
	MakeSlabPtrArray(zone,slabGradX,Vartag1);
	MakeSlabPtrArray(zone,slabGradY,Vartag2);
	MakeSlabPtrArray(zone,slabGradZ,Vartag3);
	MakeSlabPtrArrayDeltaInverse(zone,slabDeltaInverse);
	
	
	
		
		
		tag=(base->flagSaturation);
		
		switch(tag)
		{
		case 0:
	        
		Porosity= (zone->Porosity); 	
		
		
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
				
				
					
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);
				
								
				diffdia=-Disperisonflux2d(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,SurfArea,Fdiff,GradX,GradY,slabGradX,slabGradY,DeltaInverse,slabDeltaInverse);
				
				

				condia = (*(diagConvect+INDEX));

				Fsum1=0.0;
				Fsum2=0.0;
				
				for(face=0;face<4;face++)
				{
					Fsum1=Fsum1-(1.0)*Fdiff[face]+(*(Fc+face));
					Fsum2=Fsum2+*(Fc_n+face);
				}
				
				
								
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = vol*vel*Porosity+ delta_t*Fsum1;
						
						RS = vol*velP*Porosity- delta_t*Fsum2;
					break;
					case 1:
						AX = vol*vel*Porosity+delta_t*(Fsum1/2.0);
						
					        RS =vol*velP*Porosity-delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				diacoeff =Porosity*vol+delta_t*(diffdia+condia);
				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
			
	        }
		else /*3D case*/
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


				
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				
             			diffdia =-Disperisonflux3d(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,SurfArea,Fdiff,GradX,GradY,GradZ,slabGradX,slabGradY,slabGradZ,DeltaInverse,slabDeltaInverse);
				
				
				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);

				condia = *(diagConvect+INDEX);

				Fsum1=0.0;
				Fsum2=0.0;
				for(face=0;face<6;face++)
				{
					Fsum1=Fsum1-Fdiff[face]+*(Fc+face);
					Fsum2=Fsum2+*(Fc_n+face);
				}
					
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = vol*vel*Porosity+ delta_t*Fsum1;
						
						RS = vol*velP*Porosity- delta_t*Fsum2;
					break;
					case 1:
						AX = vol*vel*Porosity+delta_t*(Fsum1/2.0);
						
					        RS =vol*velP*Porosity-delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				diacoeff =Porosity*vol + delta_t*(diffdia+condia);

				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
		   break;
		  
		case 1:           //unsaturated or variably saturated one
		
		theta=zone->Flowsol->WaterContent;
	        theta_Pre=zone->PrevFlowsol->WaterContent;
		
		
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
				
				
					
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);
				
								
				diffdia=-Disperisonflux2d(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,SurfArea,Fdiff,GradX,GradY,slabGradX,slabGradY,DeltaInverse,slabDeltaInverse);
				
				

				condia = (*(diagConvect+INDEX));

				Fsum1=0.0;
				Fsum2=0.0;
				
				for(face=0;face<4;face++)
				{
					Fsum1=Fsum1-(1.0)*Fdiff[face]+(*(Fc+face));
					Fsum2=Fsum2+*(Fc_n+face);
				}
				
				
								
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = vol*vel*(*theta+INDEX )+ delta_t*Fsum1;
						
						RS = vol*velP*(*theta_Pre+INDEX)- delta_t*Fsum2;
					break;
					case 1:
						AX = vol*vel*(*theta+INDEX )+delta_t*(Fsum1/2.0);
						
					        RS =vol*velP*(*theta_Pre+INDEX)-delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}

				residue=RS-AX;
				diacoeff =(*theta+INDEX )*vol+delta_t*(diffdia+condia);
				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
			
	        }
		else /*3D case*/
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


				
				vol = *(zone->ZGeom->Volume+INDEX);
				vel=*(variable+INDEX);
				velP=*(prevTstepVar+INDEX);
				
            			 diffdia =-Disperisonflux3d(index,variable,vols,ptrSize,slabVariable,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,SurfArea,Fdiff,GradX,GradY,GradZ,slabGradX,slabGradY,slabGradZ,DeltaInverse,slabDeltaInverse);
				
				
				Fc= (upwindFlux+INDEX*6);
				Fc_n= (deferCorrFlux+INDEX*6);

				condia = *(diagConvect+INDEX);

				Fsum1=0.0;
				Fsum2=0.0;
				for(face=0;face<6;face++)
				{
					Fsum1=Fsum1-Fdiff[face]+*(Fc+face);
					Fsum2=Fsum2+*(Fc_n+face);
				}
					
				switch(base->BaseParameters->TimeStepingSchemeTag)
				{
					case 0:
						AX = vol*vel*(*theta+INDEX )+ delta_t*Fsum1;
						
						RS = vol*velP*(*theta_Pre+INDEX)- delta_t*Fsum2;
					break;
					case 1:
						AX = vol*vel*(*theta+INDEX )+delta_t*(Fsum1/2.0);
						
					        RS =vol*velP*(*theta_Pre+INDEX)-delta_t*(Fsum2/2.0);
					break;
					default:
						printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
						printf("exiting in UpdateVar\n");
						exit(1);
					break;
				}
				residue=RS-AX;
				diacoeff = vol + delta_t*(diffdia+condia);

				*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
				res=res+(residue*residue);
			}
		}
		   break;
		}
	free(deferCorrFlux);
	free(upwindFlux);
	free(diagConvect);
	return(res);
}

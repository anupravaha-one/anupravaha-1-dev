#include "APKUpdatePressureHead2d.h"

/**The function UpdatePressureHead2d calculate the  residual for different conditon in both saturated and unsaturated  porous media. It check for the unsaturated and unsaturated medium , whether along datum plane or other plane then based on deformability of medium calculate the residual. for 2d problems.


*/

double UpdatePressureHead2d(Base *base,Zone *zone)
{
	double residual=0.0,residZone=0.0,*PressureHeadDomain,pDiff=0,k1,sumDiff,diffdia;
	double pressureDiff[6];
	int NI,NJ,NK,i,j,k,iFace,index[3],indexH,*ptrSize;
	double *slabVolume[6],*slabPressure[6],*slabDiffusionCoeff[6];
	double *beta,*volumeDomain,*diffusionCoeffdomain,V_p,s1,k2=0;
	double delta_t,CoeffN[6],omega=1.0,*theta1,*theta2;
	
	double *datumflux,sumdatumflux=0,diffTheta,St,Por,S1,s2,s3,s4,*Grad;
	double diffdia1;
	int tag;
	FILE *Pomega;

	Pomega = fopen("Settings/PressureRelaxation.dat","r");// chk whether file opening or not
	fscanf(Pomega,"%lf", &omega);
	fclose(Pomega);
        
	delta_t = base->BaseParameters->TimeStep;
	PressureHeadDomain = zone->Flowsol->H;
	volumeDomain = zone->ZGeom->Volume;
	ptrSize = zone->Size ;
	beta=zone->ZGeom->Beta;
	
	St = (zone->Flowsol->Storativity);
	
	
	k1 =  *(zone->Flowsol->SatHydConductivity);
	
	Por= (zone->Porosity);
    
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	datumflux=(double *)malloc(6*(NJ-1)*(NI-1)*(NK-1)*sizeof(double));
	
	MakeSlabPtrArray(zone,slabPressure,19);
	MakeSlabPtrArray(zone,slabVolume,5);//Called for making slabs of Pressure
	  
	tag=(base->flagSaturation);
	  
	switch(tag)
	{               
		case 0:
	                
			 diffusionCoeffdomain =MakeDiffusionCoeffDomain(zone,19);
						 
			 if(base->BaseParameters->flagdat==1)//check whethwe in vertical plane or horizontal plane 
	                 datum2d(zone,datumflux,k1);
	                 else
	                 datum2d(zone,datumflux,k2);
			
			
			if(base->flagdeformable ==1 && base->BaseParameters->flagdat ==1)
			{
			
			
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
				index[0]=k,index[1]=j,index[2]=i;
				
				
				V_p =*(volumeDomain+indexH);
				s1=k1*(delta_t/(V_p*St)); 
								
			diffdia=-s1*Diffusion2dWithConstDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff,CoeffN); 
			
				diffdia=diffdia+1; //storativiy	
				
				sumDiff=0.0;
				sumdatumflux=0.0;
				pDiff =0.0;
				
				
				pDiff= *(PressureHeadDomain+indexH)-*(zone->PrevFlowsol->H+indexH);
					
				for(iFace=0;iFace<4;iFace++)
				{
					sumDiff=sumDiff+pressureDiff[iFace];
					sumdatumflux = sumdatumflux+*(datumflux+6*indexH+iFace);
				}
				
				sumdatumflux=(delta_t/(V_p*St))*sumdatumflux; 
				
				residual= -pDiff+s1*sumDiff+sumdatumflux;
				
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH) + omega*residual/diffdia;
		
				residZone=residZone+residual*residual;
			
			}
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
			}
			else if(base->flagdeformable ==1 && base->BaseParameters->flagdat ==0)
			{	
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				V_p =*(volumeDomain+indexH);
				
				s1=k1*(delta_t/(V_p*St));  
				
				
		diffdia=-s1*Diffusion2dWithConstDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff,CoeffN);
			
				
				diffdia=diffdia+1;
				sumDiff=0.0;
				pDiff =0.0;
				
				pDiff=*(PressureHeadDomain+indexH)-*(zone->PrevFlowsol->H+indexH);
					
				for(iFace=0;iFace<4;iFace++)
				sumDiff=sumDiff+pressureDiff[iFace];
				residual= -pDiff+s1*sumDiff;
		
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH)+ omega*residual/diffdia;
		
				residZone=residZone+residual*residual;
			
			}
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
			
			}
			
			else if(base->flagdeformable ==0 && base->BaseParameters->flagdat ==1)
			{
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				
				V_p =*(volumeDomain+indexH);
				s1=k1 ;
						
				diffdia=-s1*Diffusion2dWithConstDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff,CoeffN); 
			
				sumDiff=0.0;
				sumdatumflux=0.0;

				for(iFace=0;iFace<4;iFace++)
				{
					sumDiff=sumDiff+pressureDiff[iFace]; 
					sumdatumflux = sumdatumflux+*(datumflux+6*indexH+iFace);

				}
				
				residual= s1*sumDiff+sumdatumflux;
				
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH) + omega*residual/diffdia;
			
				residZone=residZone+residual*residual;
			
				}
				
				residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));

			}
			else if(base->flagdeformable ==0 && base->BaseParameters->flagdat ==0)
			{       
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				V_p =*(volumeDomain+indexH);
				s1=k1;
						
			diffdia=-s1*Diffusion2dWithConstDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff,CoeffN);
			
				sumDiff=0.0;
	
				
				for(iFace=0;iFace<4;iFace++)
				sumDiff=sumDiff+pressureDiff[iFace];
				
				
				residual= s1*sumDiff;
		
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH) + omega*residual/diffdia;
		
				residZone=residZone+residual*residual;
			
			        }
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
		         }
			break;
			
	case 1:// for h- based form
	         
		 
		 theta1=zone->Flowsol->WaterContent;
	         theta2=zone->PrevFlowsol->WaterContent;
		 diffusionCoeffdomain= MakeDiffusionCoeffDomain(zone,19);
		 MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,719);
	         datumUnsat2d(base,zone,volumeDomain,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,datumflux);
		 Grad=zone->Flowsol->GradTheta;
		
		 if(base->flagdeformable ==0 && base->BaseParameters->flagdat ==1)
		 {
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
			indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
				
			V_p =*(volumeDomain+indexH);
			
			s4=*(Grad+indexH);
			
			s2=*(theta1+indexH); 
			
			s3=*(theta2+indexH); 
		       
	
	diffdia=-Diffusion2dWithVaryingDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,pressureDiff,CoeffN); 
	  		
		
		sumDiff=0.0;
		
		sumdatumflux=0.0;
		
		for(iFace=0;iFace<4;iFace++)
		{
			sumDiff=sumDiff+pressureDiff[iFace];
			
			sumdatumflux=sumdatumflux+*(datumflux+6*indexH+iFace);
		}
	 
			
			s1=k1*(delta_t/V_p); 
			diffTheta=(s2-s3);
			diffdia1=s4+s1*diffdia;
			residual= s1*(sumDiff+sumdatumflux)-diffTheta;
		        
		
		*(PressureHeadDomain+indexH)=*(PressureHeadDomain+indexH)+residual/diffdia1;
	      	
		residZone=residZone+residual*residual;
			
		}
		
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
			
					
		}	       
	       
		  
	         else if(base->flagdeformable ==1 && base->BaseParameters->flagdat ==1)
			{
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
				index[0]=k,index[1]=j,index[2]=i;
		 		V_p =*(volumeDomain+indexH);
				s2=*(theta1+indexH);//watercontent at present 
				s3=*(theta2+indexH);//watercontent at previous time
				S1=s2*St/Por;
				s1=k1*(delta_t/(V_p*S1)); 
				
				
				diffTheta=(s2-s3)*1.0/S1;
				
		diffdia=-s1*Diffusion2dWithVaryingDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,pressureDiff,CoeffN); 
			
				diffdia=diffdia+1.0; //storativiy	
				
				sumDiff=0.0;
				sumdatumflux=0.0;
				pDiff =0.0;
				
				
				pDiff= *(PressureHeadDomain+indexH)-*(zone->PrevFlowsol->H+indexH);
					
				for(iFace=0;iFace<4;iFace++)
				{
					sumDiff=sumDiff+pressureDiff[iFace];
					sumdatumflux = sumdatumflux+*(datumflux+6*indexH+iFace);
				}
				
				sumdatumflux=(delta_t/(V_p*S1))*sumdatumflux; 
				
				residual=-pDiff+s1*sumDiff+sumdatumflux-diffTheta;
		
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH) +omega*residual/diffdia;
		
				residZone=residZone+residual*residual;
			
			}
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
			
			}
		 else if(base->flagdeformable ==1 && base->BaseParameters->flagdat ==0)
			{	
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				V_p =*(volumeDomain+indexH);
				s2=*(theta1+indexH); //watercontent at present 
				s3=*(theta2+indexH); //watercontent at previous time
				S1=s2*St/Por;
			
				s1=k1*(delta_t/(V_p*S1)); 
				diffTheta=(s2-s3)*1.0/S1;
				
				
                           diffdia=-s1*Diffusion2dWithVaryingDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,pressureDiff,CoeffN); 
			
				diffdia=diffdia+1.0;
								
					
				sumDiff=0.0;
				
				pDiff =0.0;
				
				pDiff=*(PressureHeadDomain+indexH)-*(zone->PrevFlowsol->H+indexH);
					
				for(iFace=0;iFace<4;iFace++)
				sumDiff=sumDiff+pressureDiff[iFace];
				
				
				residual= -pDiff+s1*sumDiff-diffTheta;
		
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH) + omega*residual/diffdia;
		
				residZone=residZone+residual*residual;
			
			        }
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
			
			}
		else if(base->flagdeformable ==0 && base->BaseParameters->flagdat ==0)
			{
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				index[0]=k,index[1]=j,index[2]=i;
				
				V_p =*(volumeDomain+indexH);
				
				
				s2=*(theta1+indexH); //watercontent at present 
				s3=*(theta2+indexH); //
				s1=k1*(delta_t/(V_p)); 
				diffTheta=s2-s3;
			
						
		diffdia=-s1*Diffusion2dWithVaryingDiffCoeff(index,PressureHeadDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,slabDiffusionCoeff,beta,pressureDiff,CoeffN); 
			
				sumDiff=0.0;
				
				
				for(iFace=0;iFace<4;iFace++)
				sumDiff=sumDiff+pressureDiff[iFace];
				
				
				residual= s1*sumDiff-diffTheta;
		
				*(PressureHeadDomain +indexH) = *(PressureHeadDomain+indexH) + omega*residual/diffdia;
		
				residZone=residZone+residual*residual;
			
			        }
			residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
		         }
			 
			break;
		
			
	}	
		
	free(datumflux);
	return(residZone);
}


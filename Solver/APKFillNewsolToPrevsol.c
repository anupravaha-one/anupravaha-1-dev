#include "APKFillNewsolToPrevsol.h"

/**
This function updates the old (previous time step) values to new (current time step) values
of all variables to be solved. This function loops over all the cell centers and updates the
old values with the recent values of the corresponding variable.
*/

int FillNewsolToPrevsol(Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index,count;
	int i,j,k,nophase,nop;
	Zone *zone;
	nop=base->nphase;
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
	
		if(base->BaseParameters->TimeStepingSchemeTag==1)
			CalcPrevTStepFluxes(base);
		
		for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;

			if(base->flagMPFV==0)
			{
				*((zone->PrevFlowsol)->u+index) =*((zone->Flowsol)->u+index);
				*((zone->PrevFlowsol)->v+index) =*((zone->Flowsol)->v+index);
				if((zone->Flowsol)->w!=NULL)
					*((zone->PrevFlowsol)->w+index) =*((zone->Flowsol)->w+index);
			}
			
			if(base->flaginduced==1)//rajendra
			{
			*((zone->PrevFlowsol)->TotalMfieldX+index) =*((zone->Flowsol)->TotalMfieldX+index);//rajendra
			*((zone->PrevFlowsol)->TotalMfieldY+index) =*((zone->Flowsol)->TotalMfieldY+index);
			*((zone->PrevFlowsol)->TotalMfieldZ+index) =*((zone->Flowsol)->TotalMfieldZ+index);
			}

			//added by adhi+shash for disperse phase 
			if((base->flagDroplet==1)||(base->flagParticle==1))
			{		
			/*	*((zone->PrevFlowsol)->MomentumDispersex+index) =*((zone->Flowsol)->MomentumDispersex+index);

				*((zone->PrevFlowsol)->MomentumDispersey+index) =*((zone->Flowsol)->MomentumDispersey+index);
			
				if((zone->Flowsol)->MomentumDispersez!=NULL)
				{
					*((zone->PrevFlowsol)->MomentumDispersez+index) =*((zone->Flowsol)->MomentumDispersez+index);
			
					*((zone->PrevFlowsol)->W_Disperse+index) =*((zone->Flowsol)->W_Disperse+index);
				}
			*/
				*((zone->PrevFlowsol)->VarthetaDisperse+index) =*((zone->Flowsol)->VarthetaDisperse+index);
							
				*((zone->PrevFlowsol)->U_Disperse+index) =*((zone->Flowsol)->U_Disperse+index);
				*((zone->PrevFlowsol)->V_Disperse+index) =*((zone->Flowsol)->V_Disperse+index);

				if((zone->Flowsol)->W_Disperse!=NULL)
				{	
					*((zone->PrevFlowsol)->W_Disperse+index) =*((zone->Flowsol)->W_Disperse+index);
				}
			
				
				if(((base->BaseParameters->flagVariableDia) == 1)&&(base->flagDroplet==1))
				{
					*((zone->PrevFlowsol)->ParticleDia+index*(base->BaseParameters->flagVariableDia))=*((zone->Flowsol)->ParticleDia+index*(base->BaseParameters->flagVariableDia));	
			
					*((zone->PrevFlowsol)->Tdrop+index) =*((zone->Flowsol)->Tdrop+index);
			
					if(base->BaseParameters->flagDropFuelMassFrac==1)
					{
						*((zone->PrevFlowsol)->DropFuelMassFrac+index) =*((zone->Flowsol)->DropFuelMassFrac+index);
					}
				}
			
			/*	if((base->flagDroplet==1))
				{
					*((zone->PrevFlowsol)->NumDensity+index) =*((zone->Flowsol)->NumDensity+index);
				}
			*/
			}
			
			
			else if(base->flagMPFV==1)//Added by sushil for multphase flow
			{		
				*((zone->PrevFlowsol)->Momentumx+index) =*((zone->Flowsol)->Momentumx+index);
				*((zone->PrevFlowsol)->Momentumy+index) =*((zone->Flowsol)->Momentumy+index);
				*((zone->PrevFlowsol)->Momentumz+index) =*((zone->Flowsol)->Momentumz+index);
				
				for(nophase=0;nophase<nop;nophase++)
				{
					*((zone->PrevFlowsol)->Momentumofphasex[nophase]+index) =*((zone->Flowsol)->Momentumofphasex[nophase]+index);
					
					*((zone->PrevFlowsol)->Momentumofphasey[nophase]+index) =*((zone->Flowsol)->Momentumofphasey[nophase]+index);
					
					*((zone->PrevFlowsol)->Momentumofphasez[nophase]+index) =*((zone->Flowsol)->Momentumofphasez[nophase]+index);
					
					*((zone->PrevFlowsol)->VarthetaPhase[nophase]+index) =*((zone->Flowsol)->VarthetaPhase[nophase]+index);
					//printf("\nindex=%d\tvarthea=%e",index,*((zone->PrevFlowsol)->VarthetaPhase[nophase]+index));
					
				}
				
			}
			
			if(base->flagGroundwater ==0)
			{
				if((zone->Flowsol->T)!=NULL)
					*((zone->PrevFlowsol)->T+index) =*((zone->Flowsol)->T+index);
				if((zone->Flowsol->enthalpy)!=NULL)
					*((zone->PrevFlowsol)->enthalpy+index) =*((zone->Flowsol)->enthalpy+index);
				if((zone->Flowsol->Latent_Heat_Content)!=NULL)
					*((zone->PrevFlowsol)->Latent_Heat_Content+index) =*((zone->Flowsol)->Latent_Heat_Content+index);
				
				if((zone->Flowsol->P)!=NULL)
					*((zone->PrevFlowsol)->P+index) =*((zone->Flowsol)->P+index);
			}
			else if(base->flagGroundwater==1 && base->flagVelocity==0)
			{	
				if((zone->Flowsol->H)!=NULL) // for groundwater
					*((zone->PrevFlowsol)->H+index) =*((zone->Flowsol)->H+index);	
				
				if((zone->Flowsol->WaterContent)!=NULL)// for groundwater
					*((zone->PrevFlowsol)->WaterContent+index) =*((zone->Flowsol)->WaterContent+index);
			}
				
			for(count=0;count<10;count++)
			{
				if((zone->Flowsol->Scal[count])!=NULL)
					*((zone->PrevFlowsol)->Scal[count]+index) =*((zone->Flowsol)->Scal[count]+index);
			}

			if(base->flagGroundwater ==0)
			{
                        //Shyam, START
				if((zone->Flowsol->MixFrac)!=NULL)
					*((zone->PrevFlowsol)->MixFrac+index) =*((zone->Flowsol)->MixFrac+index);
	
				if((zone->Flowsol->FuelMassFrac)!=NULL)
					*((zone->PrevFlowsol)->FuelMassFrac+index) =*((zone->Flowsol)->FuelMassFrac+index);
	
				if((zone->Flowsol->OxiMassFrac)!=NULL)
					*((zone->PrevFlowsol)->OxiMassFrac+index) =*((zone->Flowsol)->OxiMassFrac+index);
			}

			for(count=0;count<base->BaseParameters->NIngases;count++)
                        {
				if(zone->PrevFlowsol->InMassFrac[count] != NULL) 
	                        	*((zone->PrevFlowsol)->InMassFrac[count]+index) = *((zone->Flowsol)->InMassFrac[count]+index);                               
                        }
                        //Shyam, END
		}
		
		if((zone->Flowsol->flagRho)==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->Rho+index) =*((zone->Flowsol)->Rho+index);
				//printf("\nindex=%d\tdenval=%e",index,*((zone->PrevFlowsol)->Rho+index));
			}
		}
		
		if((zone->Flowsol->flagViscoeff)==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->Viscoeff+index) =*((zone->Flowsol)->Viscoeff+index);
				
				if(base->flagMPFV==1)
				{
					for(nophase=0;nophase<nop;nophase++)
					{
						*((zone->PrevFlowsol)->ViscoeffDummy[nophase]+index) =*((zone->Flowsol)->ViscoeffDummy[nophase]+index);
						//printf("\nindex=%d\tdenval=%e",index,*((zone->PrevFlowsol)->ViscoeffDummy[nophase]+index));
					}
				}
			}
		}
		
		if((zone->Flowsol->flagDiffcoeff)==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->Diffcoeff+index) =*((zone->Flowsol)->Diffcoeff+index);
			}
		}
		
		if((zone->Flowsol->flagSpecificHeat)==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->SpecificHeat+index) =*((zone->Flowsol)->SpecificHeat+index);
			}
		}
		if((zone->Flowsol->flagCPdrop)==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->CPdrop+index) =*((zone->Flowsol)->CPdrop+index);
			}
		}
		for(count=0;count<10;count++)
		{
			if(zone->Flowsol->flagTestSaturation==4)
			{
				if((zone->Flowsol->flagScalcoeff[count])==1)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
						*((zone->PrevFlowsol)->Scalcoeff[count]+index) =*((zone->Flowsol)->Scalcoeff[count]+index);
					}
				}
			}	
		}
		
		//Electrical Parameters
		if((zone->Flowsol->Potential)!=NULL)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->Potential+index) =*((zone->Flowsol)->Potential+index);
			}
		}
		if((zone->Flowsol->RhoIon)!=NULL)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*((zone->PrevFlowsol)->RhoIon+index) =*((zone->Flowsol)->RhoIon+index);
			}
		}
		//End of Electrical Parameters
		
		//LMMHD parameters(naveen)
		if(base->flagLMMHD==1)
		{
			if((zone->Flowsol->PotentialMHD)!=NULL)
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*((zone->PrevFlowsol)->PotentialMHD+index) =*((zone->Flowsol)->PotentialMHD+index);
                                }       
			}
                
		}
	
		if(base->flaginduced==1)//rajendra
		{
			if((zone->Flowsol->Psi)!=NULL)
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*((zone->PrevFlowsol)->Psi+index) =*((zone->Flowsol)->Psi+index);
				}
			}
		}	
		//End of LMMHD parameters
		if(base->flagGroundwater ==0)
		{
			UpdateBoundarySlabs(zone);
			UpdateOldZSlabs(zone);
		}	
	}
	return 0;
}

/**
This function updates the old (previous time step) values to new (current time step) values
of K and Epsilon for the turbulent flow problem.
*/

int FillNewsolToPrevsolKepsilon(Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index,count;
	int i,j,k;
	Zone *zone;
	
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
	
		for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
			for(count=0;count<2;count++)
			{
				if((zone->Flowsol->Scal[count])!=NULL)
					*((zone->PrevFlowsol)->Scal[count]+index) =*((zone->Flowsol)->Scal[count]+index);
			}
		}
		UpdateOldZSlabsKEpsilon(zone);
	}
	return 0;
}



#include "APKMakeDomainArray.h"
int FillNewsolToPrevsol_T(Base *base,int vartag)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k;
	Zone *zone;
	double *old,*new;
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		old=MakeDomainArray(zone,vartag*10);
		new=MakeDomainArray(zone,vartag);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
	
		for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			*(old+index) =*(new+index);
		}
	}
	return 0;
}



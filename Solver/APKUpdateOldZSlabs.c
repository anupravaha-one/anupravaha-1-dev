#include "APKUpdateOldZSlabs.h"
/**
This function UpdateOldZSlabs() is called to fill the previous timestep slab values (uOld ,umassOld, vOld, vmassOld etc)
with current timestep slab values(u,umass,v,vmass etc). This function also does the above operation for Interface slabs
also 
*/
int UpdateOldZSlabs(Zone *zone)
{
	int k,j,i,NK,NJ,NI,iface,*ptrSize,iScal,nophase,nop;
	char tag,name[33];
	void *ptr;
	ZoneBoundary *bcPtr;
	Interface *IfacePtr;
	nop = zone->Nphase;		
	ptrSize = zone->Size;
		
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	for(iface=0;iface<6;iface++)
	{
		ptr = FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		
		if(tag=='b')
		{
			bcPtr=(ZoneBoundary *)ptr;
			
			switch(iface+1)
			{
				case 1:

					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							*(bcPtr->ZSlab->uOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->u+k*(NJ-1)+j);

						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							*(bcPtr->ZSlab->vOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->v+k*(NJ-1)+j);

						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							*(bcPtr->ZSlab->wOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->w+k*(NJ-1)+j);
						
						if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							*(bcPtr->ZSlab->TOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->T+k*(NJ-1)+j);
							
						if((zone->Flowsol->enthalpy)!=NULL && (bcPtr->ZSlab->enthalpyOld)!=NULL)
							*(bcPtr->ZSlab->enthalpyOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->enthalpy+k*(NJ-1)+j);


						if((zone->Flowsol->umass)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->umass+k*(NJ-1)+j);

						if((zone->Flowsol->vmass)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->vmass+k*(NJ-1)+j);
							
						if((zone->Flowsol->wmass)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
						*(bcPtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->wmass+k*(NJ-1)+j);

					/*	if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->uOld+k*(NJ-1)+j);

						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->vOld+k*(NJ-1)+j);
							
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->wOld+k*(NJ-1)+j);
					*/


						if((zone->Flowsol->TotalMfieldX)!=NULL && (bcPtr->ZSlab->TotalMfieldXold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldXold+k*(NJ-1)+j)=*(bcPtr->ZSlab->TotalMfieldX+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (bcPtr->ZSlab->TotalMfieldYold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldYold+k*(NJ-1)+j)=*(bcPtr->ZSlab->TotalMfieldY+k*(NJ-1)+j);//raj


						if((zone->Flowsol->TotalMfieldZ)!=NULL && (bcPtr->ZSlab->TotalMfieldZold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldZold+k*(NJ-1)+j)=*(bcPtr->ZSlab->TotalMfieldZ+k*(NJ-1)+j);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (bcPtr->ZSlab->PotentialMHDold)!=NULL)
							*(bcPtr->ZSlab->PotentialMHDold+k*(NJ-1)+j)=*(bcPtr->ZSlab->PotentialMHD+k*(NJ-1)+j);//raj
						if((zone->Flowsol->Psi)!=NULL && (bcPtr->ZSlab->Psiold)!=NULL)
							*(bcPtr->ZSlab->Psiold+k*(NJ-1)+j)=*(bcPtr->ZSlab->Psi+k*(NJ-1)+j);//raj
						
						//Added by sushil for multiphase flow	
						if((zone->Flowsol->Momentumx)!=NULL && (bcPtr->ZSlab->MomentumxOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumx+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumy)!=NULL && (bcPtr->ZSlab->MomentumyOld)!=NULL)
							*(bcPtr->ZSlab->MomentumyOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumy+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (bcPtr->ZSlab->MomentumzOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumz+k*(NJ-1)+j);
							
						
								
						if((zone->Flowsol->Momentumxstar)!=NULL && (bcPtr->ZSlab->MomentumxstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxstarOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumxstar+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (bcPtr->ZSlab->MomentumystarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumystarOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumystar+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (bcPtr->ZSlab->MomentumzstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzstarOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumzstar+k*(NJ-1)+j);
							
						for(nophase=0;nophase<nop;nophase++)
						{	
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasexOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumofphasex[nophase]+k*(NJ-1)+j);
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
								
								*(bcPtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumofphasey[nophase]+k*(NJ-1)+j);
									
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasezOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumofphasez[nophase]+k*(NJ-1)+j);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (bcPtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->VarthetaPhase[nophase]+k*(NJ-1)+j);
						}	
						//end
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->U_Disperse+k*(NJ-1)+j);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->V_Disperse+k*(NJ-1)+j);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->W_Disperse+k*(NJ-1)+j);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NJ-1)+j);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NJ-1)+j);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Tdrop+k*(NJ-1)+j);
			

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 2:
					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							*(bcPtr->ZSlab->uOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->u+k*(NJ-1)+j);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							*(bcPtr->ZSlab->vOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->v+k*(NJ-1)+j);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							*(bcPtr->ZSlab->wOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->w+k*(NJ-1)+j);
						
						if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							*(bcPtr->ZSlab->TOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->T+k*(NJ-1)+j);
						
						if((zone->Flowsol->enthalpy)!=NULL && (bcPtr->ZSlab->enthalpyOld)!=NULL)
							*(bcPtr->ZSlab->enthalpyOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->enthalpy+k*(NJ-1)+j);
							
						if((zone->Flowsol->umass)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->umass+k*(NJ-1)+j);
						if((zone->Flowsol->vmass)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->vmass+k*(NJ-1)+j);
							
						if((zone->Flowsol->wmass)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
						*(bcPtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->wmass+k*(NJ-1)+j);

					/*	if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->uOld+k*(NJ-1)+j);

						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->vOld+k*(NJ-1)+j);
							
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->wOld+k*(NJ-1)+j);
					*/

						if((zone->Flowsol->TotalMfieldX)!=NULL && (bcPtr->ZSlab->TotalMfieldXold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldXold+k*(NJ-1)+j)=*(bcPtr->ZSlab->TotalMfieldX+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (bcPtr->ZSlab->TotalMfieldYold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldYold+k*(NJ-1)+j)=*(bcPtr->ZSlab->TotalMfieldY+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (bcPtr->ZSlab->TotalMfieldZold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldZold+k*(NJ-1)+j)=*(bcPtr->ZSlab->TotalMfieldZ+k*(NJ-1)+j);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (bcPtr->ZSlab->PotentialMHDold)!=NULL)
							*(bcPtr->ZSlab->PotentialMHDold+k*(NJ-1)+j)=*(bcPtr->ZSlab->PotentialMHD+k*(NJ-1)+j);//raj
						if((zone->Flowsol->Psi)!=NULL && (bcPtr->ZSlab->Psiold)!=NULL)
							*(bcPtr->ZSlab->Psiold+k*(NJ-1)+j)=*(bcPtr->ZSlab->Psi+k*(NJ-1)+j);//raj

						
						//Added by sushil for multiphase flow
						if((zone->Flowsol->Momentumx)!=NULL && (bcPtr->ZSlab->MomentumxOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumx+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumy)!=NULL && (bcPtr->ZSlab->MomentumyOld)!=NULL)
							*(bcPtr->ZSlab->MomentumyOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumy+k*(NJ-1)+j);	
							
						if((zone->Flowsol->Momentumz)!=NULL && (bcPtr->ZSlab->MomentumzOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumz+k*(NJ-1)+j);	
							
							
						if((zone->Flowsol->Momentumxstar)!=NULL && (bcPtr->ZSlab->MomentumxstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxstarOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumxstar+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (bcPtr->ZSlab->MomentumystarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumystarOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumystar+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (bcPtr->ZSlab->MomentumzstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzstarOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumzstar+k*(NJ-1)+j);
							
							
						for(nophase=0;nophase<nop;nophase++)
						{	
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasexOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumofphasex[nophase]+k*(NJ-1)+j);
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
								
								*(bcPtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumofphasey[nophase]+k*(NJ-1)+j);
									
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasezOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Momentumofphasez[nophase]+k*(NJ-1)+j);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (bcPtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NJ-1)+j)=*(bcPtr->ZSlab->VarthetaPhase[nophase]+k*(NJ-1)+j);
						}		
						//end
						//Added by adhi+shash for disperse phase flow//mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->U_Disperse+k*(NJ-1)+j);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->V_Disperse+k*(NJ-1)+j);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->W_Disperse+k*(NJ-1)+j);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NJ-1)+j);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NJ-1)+j);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Tdrop+k*(NJ-1)+j);

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 3:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							*(bcPtr->ZSlab->uOld+k*(NI-1)+i)=*(bcPtr->ZSlab->u+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							*(bcPtr->ZSlab->vOld+k*(NI-1)+i)=*(bcPtr->ZSlab->v+k*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							*(bcPtr->ZSlab->wOld+k*(NI-1)+i)=*(bcPtr->ZSlab->w+k*(NI-1)+i);
						if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							*(bcPtr->ZSlab->TOld+k*(NI-1)+i)=*(bcPtr->ZSlab->T+k*(NI-1)+i);
						if((zone->Flowsol->enthalpy)!=NULL && (bcPtr->ZSlab->enthalpyOld)!=NULL)
							*(bcPtr->ZSlab->enthalpyOld+k*(NI-1)+i)=*(bcPtr->ZSlab->enthalpy+k*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->umass+k*(NI-1)+i);
						if((zone->Flowsol->vmass)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->vmass+k*(NI-1)+i);	
						if((zone->Flowsol->wmass)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
						*(bcPtr->ZSlab->wmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->wmass+k*(NI-1)+i);	


					/*	if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->uOld+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->vOld+k*(NI-1)+i);	
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->wOld+k*(NI-1)+i);	
					*/
						if((zone->Flowsol->TotalMfieldX)!=NULL && (bcPtr->ZSlab->TotalMfieldXold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldXold+k*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldX+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (bcPtr->ZSlab->TotalMfieldYold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldYold+k*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldY+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (bcPtr->ZSlab->TotalMfieldZold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldZold+k*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldZ+k*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (bcPtr->ZSlab->PotentialMHDold)!=NULL)
							*(bcPtr->ZSlab->PotentialMHDold+k*(NI-1)+i)=*(bcPtr->ZSlab->PotentialMHD+k*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (bcPtr->ZSlab->Psiold)!=NULL)
							*(bcPtr->ZSlab->Psiold+k*(NI-1)+i)=*(bcPtr->ZSlab->Psi+k*(NI-1)+i);//raj

						
						//Added by sushil for multiphase flow	
						if((zone->Flowsol->Momentumx)!=NULL && (bcPtr->ZSlab->MomentumxOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumx+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumy)!=NULL && (bcPtr->ZSlab->MomentumyOld)!=NULL)
							*(bcPtr->ZSlab->MomentumyOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumy+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (bcPtr->ZSlab->MomentumzOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumz+k*(NI-1)+i);	
							
							
						if((zone->Flowsol->Momentumxstar)!=NULL && (bcPtr->ZSlab->MomentumxstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxstarOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumxstar+k*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (bcPtr->ZSlab->MomentumystarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumystarOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumystar+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (bcPtr->ZSlab->MomentumzstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzstarOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumzstar+k*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{	
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasexOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasex[nophase]+k*(NI-1)+i);
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
								
								*(bcPtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasey[nophase]+k*(NI-1)+i);
									
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasezOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasez[nophase]+k*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (bcPtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaPhase[nophase]+k*(NI-1)+i);
						}			
							//end
						//Added by adhi+shash for disperse phase flow	//mrunalini

						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+k*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+k*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+k*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+k*(NI-1)+i);

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 4:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							*(bcPtr->ZSlab->uOld+k*(NI-1)+i)=*(bcPtr->ZSlab->u+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							*(bcPtr->ZSlab->vOld+k*(NI-1)+i)=*(bcPtr->ZSlab->v+k*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							*(bcPtr->ZSlab->wOld+k*(NI-1)+i)=*(bcPtr->ZSlab->w+k*(NI-1)+i);
						if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							*(bcPtr->ZSlab->TOld+k*(NI-1)+i)=*(bcPtr->ZSlab->T+k*(NI-1)+i);
						if((zone->Flowsol->enthalpy)!=NULL && (bcPtr->ZSlab->enthalpyOld)!=NULL)
							*(bcPtr->ZSlab->enthalpyOld+k*(NI-1)+i)=*(bcPtr->ZSlab->enthalpy+k*(NI-1)+i);

						if((zone->Flowsol->umass)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->umass+k*(NI-1)+i);
						if((zone->Flowsol->vmass)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->vmass+k*(NI-1)+i);
						if((zone->Flowsol->wmass)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
						*(bcPtr->ZSlab->wmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->wmass+k*(NI-1)+i);

					/*	if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->uOld+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->vOld+k*(NI-1)+i);	
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->wOld+k*(NI-1)+i);	
					*/
						if((zone->Flowsol->TotalMfieldX)!=NULL && (bcPtr->ZSlab->TotalMfieldXold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldXold+k*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldX+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (bcPtr->ZSlab->TotalMfieldYold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldYold+k*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldY+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (bcPtr->ZSlab->TotalMfieldZold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldZold+k*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldZ+k*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (bcPtr->ZSlab->PotentialMHDold)!=NULL)
							*(bcPtr->ZSlab->PotentialMHDold+k*(NI-1)+i)=*(bcPtr->ZSlab->PotentialMHD+k*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (bcPtr->ZSlab->Psiold)!=NULL)
							*(bcPtr->ZSlab->Psiold+k*(NI-1)+i)=*(bcPtr->ZSlab->Psi+k*(NI-1)+i);//raj


						
						//Added by sushil for multiphase flow
						
						
						if((zone->Flowsol->Momentumx)!=NULL && (bcPtr->ZSlab->MomentumxOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumx+k*(NI-1)+i);
							
						
						if((zone->Flowsol->Momentumy)!=NULL && (bcPtr->ZSlab->MomentumyOld)!=NULL)
							*(bcPtr->ZSlab->MomentumyOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumy+k*(NI-1)+i);
							
						
						if((zone->Flowsol->Momentumz)!=NULL && (bcPtr->ZSlab->MomentumzOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumz+k*(NI-1)+i);	
						
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (bcPtr->ZSlab->MomentumxstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxstarOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumxstar+k*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (bcPtr->ZSlab->MomentumystarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumystarOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumystar+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (bcPtr->ZSlab->MomentumzstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzstarOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumzstar+k*(NI-1)+i);
							
							
						for(nophase=0;nophase<nop;nophase++)
						{	
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasexOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasex[nophase]+k*(NI-1)+i);
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
								
								*(bcPtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasey[nophase]+k*(NI-1)+i);
									
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasezOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasez[nophase]+k*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (bcPtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaPhase[nophase]+k*(NI-1)+i);
						}				
						
						//end	
						//Added by adhi+shash for disperse phase flow	//mrunalini

						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+k*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+k*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+k*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+k*(NI-1)+i);

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 5:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							*(bcPtr->ZSlab->uOld+j*(NI-1)+i)=*(bcPtr->ZSlab->u+j*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							*(bcPtr->ZSlab->vOld+j*(NI-1)+i)=*(bcPtr->ZSlab->v+j*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							*(bcPtr->ZSlab->wOld+j*(NI-1)+i)=*(bcPtr->ZSlab->w+j*(NI-1)+i);
						
						if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							*(bcPtr->ZSlab->TOld+j*(NI-1)+i)=*(bcPtr->ZSlab->T+j*(NI-1)+i);
						if((zone->Flowsol->enthalpy)!=NULL && (bcPtr->ZSlab->enthalpyOld)!=NULL)
							*(bcPtr->ZSlab->enthalpyOld+j*(NI-1)+i)=*(bcPtr->ZSlab->enthalpy+j*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->umass+j*(NI-1)+i);
							
						if((zone->Flowsol->vmass)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->vmass+j*(NI-1)+i);
						if((zone->Flowsol->wmass)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->wmass+j*(NI-1)+i);


					/*	if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->uOld+j*(NI-1)+i);
							
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->vOld+j*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->wOld+j*(NI-1)+i);

					*/

						if((zone->Flowsol->TotalMfieldX)!=NULL && (bcPtr->ZSlab->TotalMfieldXold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldXold+j*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldX+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (bcPtr->ZSlab->TotalMfieldYold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldYold+j*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldY+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (bcPtr->ZSlab->TotalMfieldZold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldZold+j*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldZ+j*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (bcPtr->ZSlab->PotentialMHDold)!=NULL)
							*(bcPtr->ZSlab->PotentialMHDold+j*(NI-1)+i)=*(bcPtr->ZSlab->PotentialMHD+j*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (bcPtr->ZSlab->Psiold)!=NULL)
							*(bcPtr->ZSlab->Psiold+j*(NI-1)+i)=*(bcPtr->ZSlab->Psi+j*(NI-1)+i);//raj

						
						//Added by sushil for multiphase flow
						if((zone->Flowsol->Momentumx)!=NULL && (bcPtr->ZSlab->MomentumxOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumx+j*(NI-1)+i);
						
						if((zone->Flowsol->Momentumy)!=NULL && (bcPtr->ZSlab->MomentumyOld)!=NULL)
							*(bcPtr->ZSlab->MomentumyOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumy+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumz)!=NULL && (bcPtr->ZSlab->MomentumzOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumz+j*(NI-1)+i);		
						
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (bcPtr->ZSlab->MomentumxstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxstarOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumxstar+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (bcPtr->ZSlab->MomentumystarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumystarOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumystar+j*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (bcPtr->ZSlab->MomentumzstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzstarOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumzstar+j*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{	
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasexOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasex[nophase]+j*(NI-1)+i);
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
								
								*(bcPtr->ZSlab->MomentumofphaseyOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasey[nophase]+j*(NI-1)+i);
									
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasezOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasez[nophase]+j*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (bcPtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->VarthetaPhaseOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaPhase[nophase]+j*(NI-1)+i);
						}				
						
						//end
						//Added by adhi+shash for disperse phase flow	//mrunalini

						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+j*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+j*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+j*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+j*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+j*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+j*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+j*(NI-1)+i);


						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
				case 6:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							*(bcPtr->ZSlab->uOld+j*(NI-1)+i)=*(bcPtr->ZSlab->u+j*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							*(bcPtr->ZSlab->vOld+j*(NI-1)+i)=*(bcPtr->ZSlab->v+j*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							*(bcPtr->ZSlab->wOld+j*(NI-1)+i)=*(bcPtr->ZSlab->w+j*(NI-1)+i);
						
						if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							*(bcPtr->ZSlab->TOld+j*(NI-1)+i)=*(bcPtr->ZSlab->T+j*(NI-1)+i);
							
						if((zone->Flowsol->enthalpy)!=NULL && (bcPtr->ZSlab->enthalpyOld)!=NULL)
							*(bcPtr->ZSlab->enthalpyOld+j*(NI-1)+i)=*(bcPtr->ZSlab->enthalpy+j*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->umass+j*(NI-1)+i);
							
						if((zone->Flowsol->vmass)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->vmass+j*(NI-1)+i);
							
						if((zone->Flowsol->wmass)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
						*(bcPtr->ZSlab->wmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->wmass+j*(NI-1)+i);

					/*	if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							*(bcPtr->ZSlab->umassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->uOld+j*(NI-1)+i);
							
						if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							*(bcPtr->ZSlab->vmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->vOld+j*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							*(bcPtr->ZSlab->wmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->wOld+j*(NI-1)+i);
					*/


						if((zone->Flowsol->TotalMfieldX)!=NULL && (bcPtr->ZSlab->TotalMfieldXold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldXold+j*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldX+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (bcPtr->ZSlab->TotalMfieldYold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldYold+j*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldY+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (bcPtr->ZSlab->TotalMfieldZold)!=NULL)
							*(bcPtr->ZSlab->TotalMfieldZold+j*(NI-1)+i)=*(bcPtr->ZSlab->TotalMfieldZ+j*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (bcPtr->ZSlab->PotentialMHDold)!=NULL)
							*(bcPtr->ZSlab->PotentialMHDold+j*(NI-1)+i)=*(bcPtr->ZSlab->PotentialMHD+j*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (bcPtr->ZSlab->Psiold)!=NULL)
							*(bcPtr->ZSlab->Psiold+j*(NI-1)+i)=*(bcPtr->ZSlab->Psi+j*(NI-1)+i);//raj

						
						//Added by sushil for multipahse flow
						if((zone->Flowsol->Momentumx)!=NULL && (bcPtr->ZSlab->MomentumxOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumx+j*(NI-1)+i);
						
						if((zone->Flowsol->Momentumy)!=NULL && (bcPtr->ZSlab->MomentumyOld)!=NULL)
							*(bcPtr->ZSlab->MomentumyOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumy+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumz)!=NULL && (bcPtr->ZSlab->MomentumzOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumz+j*(NI-1)+i);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (bcPtr->ZSlab->MomentumxstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumxstarOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumxstar+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (bcPtr->ZSlab->MomentumystarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumystarOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumystar+j*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (bcPtr->ZSlab->MomentumzstarOld)!=NULL)
							*(bcPtr->ZSlab->MomentumzstarOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumzstar+j*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{	
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasexOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasex[nophase]+j*(NI-1)+i);
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
								
								*(bcPtr->ZSlab->MomentumofphaseyOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasey[nophase]+j*(NI-1)+i);
									
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (bcPtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->MomentumofphasezOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->Momentumofphasez[nophase]+j*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (bcPtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(bcPtr->ZSlab->VarthetaPhaseOld[nophase]+j*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaPhase[nophase]+j*(NI-1)+i);
						}					
						
						//end
						//Added by adhi+shash for disperse phase flow	//mrunalini

						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+j*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+j*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+j*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+j*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+j*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+j*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+j*(NI-1)+i);

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
			}
		}
		else
		if(tag=='i')
		{
			IfacePtr=(Interface *)ptr;
			
			switch(iface+1)
			{
				case 1:

					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->uOld)!=NULL)
							*(IfacePtr->ZSlab->uOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->u+k*(NJ-1)+j);
						
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vOld)!=NULL)
							*(IfacePtr->ZSlab->vOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->v+k*(NJ-1)+j);

						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wOld)!=NULL)
							*(IfacePtr->ZSlab->wOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->w+k*(NJ-1)+j);
						
						if((zone->Flowsol->T)!=NULL && (IfacePtr->ZSlab->TOld)!=NULL)
							*(IfacePtr->ZSlab->TOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->T+k*(NJ-1)+j);
							
						if((zone->Flowsol->enthalpy)!=NULL && (IfacePtr->ZSlab->enthalpyOld)!=NULL)
							*(IfacePtr->ZSlab->enthalpyOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->enthalpy+k*(NJ-1)+j);


						if((zone->Flowsol->umass)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->umass+k*(NJ-1)+j);
							
						if((zone->Flowsol->vmass)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->vmass+k*(NJ-1)+j);
							
						if((zone->Flowsol->wmass)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
						*(IfacePtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->wmass+k*(NJ-1)+j);


					/*	if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->uOld+k*(NJ-1)+j);
							
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->vOld+k*(NJ-1)+j);
							
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
							*(IfacePtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->wOld+k*(NJ-1)+j);
					*/
						if((zone->Flowsol->TotalMfieldX)!=NULL && (IfacePtr->ZSlab->TotalMfieldXold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldXold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->TotalMfieldX+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (IfacePtr->ZSlab->TotalMfieldYold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldYold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->TotalMfieldY+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (IfacePtr->ZSlab->TotalMfieldZold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldZold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->TotalMfieldZ+k*(NJ-1)+j);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (IfacePtr->ZSlab->PotentialMHDold)!=NULL)
							*(IfacePtr->ZSlab->PotentialMHDold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->PotentialMHD+k*(NJ-1)+j);//raj
						if((zone->Flowsol->Psi)!=NULL && (IfacePtr->ZSlab->Psiold)!=NULL)
							*(IfacePtr->ZSlab->Psiold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Psi+k*(NJ-1)+j);//raj

						
						//Added by sushil for multiphase flow	
						if((zone->Flowsol->Momentumx)!=NULL && (IfacePtr->ZSlab->MomentumxOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumx+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumy)!=NULL && (IfacePtr->ZSlab->MomentumyOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumyOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumy+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (IfacePtr->ZSlab->MomentumzOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumz+k*(NJ-1)+j);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (IfacePtr->ZSlab->MomentumxstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxstarOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumxstar+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (IfacePtr->ZSlab->MomentumystarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumystarOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumystar+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (IfacePtr->ZSlab->MomentumzstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzstarOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumzstar+k*(NJ-1)+j);
							
						for(nophase=0;nophase<nop;nophase++)
						{
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasexOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumofphasex[nophase]+k*(NJ-1)+j);
							
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumofphasey[nophase]+k*(NJ-1)+j);
							
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasezOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumofphasez[nophase]+k*(NJ-1)+j);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (IfacePtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->VarthetaPhase[nophase]+k*(NJ-1)+j);
						}	
						
						//end
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (IfacePtr->ZSlab->U_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->U_DisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->U_Disperse+k*(NJ-1)+j);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (IfacePtr->ZSlab->V_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->V_DisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->V_Disperse+k*(NJ-1)+j);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (IfacePtr->ZSlab->W_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->W_DisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->W_Disperse+k*(NJ-1)+j);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (IfacePtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(IfacePtr->ZSlab->VarthetaDisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->VarthetaDisperse+k*(NJ-1)+j);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (IfacePtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(IfacePtr->ZSlab->DropFuelMassFracOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->DropFuelMassFrac+k*(NJ-1)+j);

						if((zone->Flowsol->Tdrop)!=NULL && (IfacePtr->ZSlab->TdropOld)!=NULL)
							*(IfacePtr->ZSlab->TdropOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Tdrop+k*(NJ-1)+j);		
								
						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 2:
					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->uOld)!=NULL)
							*(IfacePtr->ZSlab->uOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->u+k*(NJ-1)+j);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vOld)!=NULL)
							*(IfacePtr->ZSlab->vOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->v+k*(NJ-1)+j);
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wOld)!=NULL)
							*(IfacePtr->ZSlab->wOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->w+k*(NJ-1)+j);
						
						if((zone->Flowsol->T)!=NULL && (IfacePtr->ZSlab->TOld)!=NULL)
							*(IfacePtr->ZSlab->TOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->T+k*(NJ-1)+j);
						
						if((zone->Flowsol->enthalpy)!=NULL && (IfacePtr->ZSlab->enthalpyOld)!=NULL)
							*(IfacePtr->ZSlab->enthalpyOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->enthalpy+k*(NJ-1)+j);
							
						if((zone->Flowsol->umass)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->umass+k*(NJ-1)+j);
							
						if((zone->Flowsol->vmass)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->vmass+k*(NJ-1)+j);
							
						if((zone->Flowsol->wmass)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
						*(IfacePtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->wmass+k*(NJ-1)+j);
					/*

						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->uOld+k*(NJ-1)+j);
							
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->vOld+k*(NJ-1)+j);
							
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
							*(IfacePtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->wOld+k*(NJ-1)+j);

					*/
						if((zone->Flowsol->TotalMfieldX)!=NULL && (IfacePtr->ZSlab->TotalMfieldXold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldXold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->TotalMfieldX+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (IfacePtr->ZSlab->TotalMfieldYold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldYold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->TotalMfieldY+k*(NJ-1)+j);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (IfacePtr->ZSlab->TotalMfieldZold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldZold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->TotalMfieldZ+k*(NJ-1)+j);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (IfacePtr->ZSlab->PotentialMHDold)!=NULL)
							*(IfacePtr->ZSlab->PotentialMHDold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->PotentialMHD+k*(NJ-1)+j);//raj
						if((zone->Flowsol->Psi)!=NULL && (IfacePtr->ZSlab->Psiold)!=NULL)
							*(IfacePtr->ZSlab->Psiold+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Psi+k*(NJ-1)+j);//raj

						
						//Added by sushil for multiphase flow
						
						if((zone->Flowsol->Momentumx)!=NULL && (IfacePtr->ZSlab->MomentumxOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumx+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumy)!=NULL && (IfacePtr->ZSlab->MomentumyOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumyOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumy+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (IfacePtr->ZSlab->MomentumzOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumz+k*(NJ-1)+j);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (IfacePtr->ZSlab->MomentumxstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxstarOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumxstar+k*(NJ-1)+j);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (IfacePtr->ZSlab->MomentumystarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumystarOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumystar+k*(NJ-1)+j);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (IfacePtr->ZSlab->MomentumzstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzstarOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumzstar+k*(NJ-1)+j);
							
						for(nophase=0;nophase<nop;nophase++)
						{
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasexOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumofphasex[nophase]+k*(NJ-1)+j);
							
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumofphasey[nophase]+k*(NJ-1)+j);
							
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasezOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Momentumofphasez[nophase]+k*(NJ-1)+j);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (IfacePtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->VarthetaPhase[nophase]+k*(NJ-1)+j);
						}	
						
						
						//end
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (IfacePtr->ZSlab->U_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->U_DisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->U_Disperse+k*(NJ-1)+j);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (IfacePtr->ZSlab->V_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->V_DisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->V_Disperse+k*(NJ-1)+j);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (IfacePtr->ZSlab->W_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->W_DisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->W_Disperse+k*(NJ-1)+j);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (IfacePtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(IfacePtr->ZSlab->VarthetaDisperseOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->VarthetaDisperse+k*(NJ-1)+j);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (IfacePtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(IfacePtr->ZSlab->DropFuelMassFracOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->DropFuelMassFrac+k*(NJ-1)+j);

						if((zone->Flowsol->Tdrop)!=NULL && (IfacePtr->ZSlab->TdropOld)!=NULL)
							*(IfacePtr->ZSlab->TdropOld+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Tdrop+k*(NJ-1)+j);
			
						

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 3:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->uOld)!=NULL)
							*(IfacePtr->ZSlab->uOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->u+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vOld)!=NULL)
							*(IfacePtr->ZSlab->vOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->v+k*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wOld)!=NULL)
							*(IfacePtr->ZSlab->wOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->w+k*(NI-1)+i);
						if((zone->Flowsol->T)!=NULL && (IfacePtr->ZSlab->TOld)!=NULL)
							*(IfacePtr->ZSlab->TOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->T+k*(NI-1)+i);
							
						if((zone->Flowsol->enthalpy)!=NULL && (IfacePtr->ZSlab->enthalpyOld)!=NULL)
							*(IfacePtr->ZSlab->enthalpyOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->enthalpy+k*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->umass+k*(NI-1)+i);
						if((zone->Flowsol->vmass)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->vmass+k*(NI-1)+i);
							
						if((zone->Flowsol->wmass)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
						*(IfacePtr->ZSlab->wmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->wmass+k*(NI-1)+i);


					/*	if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->uOld+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->vOld+k*(NI-1)+i);
							
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
							*(IfacePtr->ZSlab->wmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->wOld+k*(NI-1)+i);
					*/

						if((zone->Flowsol->TotalMfieldX)!=NULL && (IfacePtr->ZSlab->TotalMfieldXold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldXold+k*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldX+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (IfacePtr->ZSlab->TotalMfieldYold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldYold+k*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldY+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (IfacePtr->ZSlab->TotalMfieldZold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldZold+k*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldZ+k*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (IfacePtr->ZSlab->PotentialMHDold)!=NULL)
							*(IfacePtr->ZSlab->PotentialMHDold+k*(NI-1)+i)=*(IfacePtr->ZSlab->PotentialMHD+k*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (IfacePtr->ZSlab->Psiold)!=NULL)
							*(IfacePtr->ZSlab->Psiold+k*(NI-1)+i)=*(IfacePtr->ZSlab->Psi+k*(NI-1)+i);//raj

						
						
						//Added by sushil for multiphase flow	
							
						if((zone->Flowsol->Momentumx)!=NULL && (IfacePtr->ZSlab->MomentumxOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumx+k*(NI-1)+i);
							
						if((zone->Flowsol->Momentumy)!=NULL && (IfacePtr->ZSlab->MomentumyOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumyOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumy+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (IfacePtr->ZSlab->MomentumzOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumz+k*(NI-1)+i);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (IfacePtr->ZSlab->MomentumxstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxstarOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumxstar+k*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (IfacePtr->ZSlab->MomentumystarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumystarOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumystar+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (IfacePtr->ZSlab->MomentumzstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzstarOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumzstar+k*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasexOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasex[nophase]+k*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasey[nophase]+k*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasezOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasez[nophase]+k*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (IfacePtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaPhase[nophase]+k*(NI-1)+i);
						}	
						
							//end
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (IfacePtr->ZSlab->U_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->U_DisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->U_Disperse+k*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (IfacePtr->ZSlab->V_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->V_DisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->V_Disperse+k*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (IfacePtr->ZSlab->W_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->W_DisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->W_Disperse+k*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (IfacePtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(IfacePtr->ZSlab->VarthetaDisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaDisperse+k*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (IfacePtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(IfacePtr->ZSlab->DropFuelMassFracOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->DropFuelMassFrac+k*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (IfacePtr->ZSlab->TdropOld)!=NULL)
							*(IfacePtr->ZSlab->TdropOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Tdrop+k*(NI-1)+i);
			

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 4:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->uOld)!=NULL)
							*(IfacePtr->ZSlab->uOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->u+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vOld)!=NULL)
							*(IfacePtr->ZSlab->vOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->v+k*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wOld)!=NULL)
							*(IfacePtr->ZSlab->wOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->w+k*(NI-1)+i);
						
						if((zone->Flowsol->T)!=NULL && (IfacePtr->ZSlab->TOld)!=NULL)
							*(IfacePtr->ZSlab->TOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->T+k*(NI-1)+i);
							
						if((zone->Flowsol->enthalpy)!=NULL && (IfacePtr->ZSlab->enthalpyOld)!=NULL)
							*(IfacePtr->ZSlab->enthalpyOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->enthalpy+k*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->umass+k*(NI-1)+i);
							
						if((zone->Flowsol->vmass)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->vmass+k*(NI-1)+i);	
							
						if((zone->Flowsol->wmass)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
						*(IfacePtr->ZSlab->wmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->wmass+k*(NI-1)+i);	

					/*	
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->uOld+k*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->vOld+k*(NI-1)+i);
							
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
							*(IfacePtr->ZSlab->wmassOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->wOld+k*(NI-1)+i);
					*/

						if((zone->Flowsol->TotalMfieldX)!=NULL && (IfacePtr->ZSlab->TotalMfieldXold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldXold+k*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldX+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (IfacePtr->ZSlab->TotalMfieldYold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldYold+k*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldY+k*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (IfacePtr->ZSlab->TotalMfieldZold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldZold+k*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldZ+k*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (IfacePtr->ZSlab->PotentialMHDold)!=NULL)
							*(IfacePtr->ZSlab->PotentialMHDold+k*(NI-1)+i)=*(IfacePtr->ZSlab->PotentialMHD+k*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (IfacePtr->ZSlab->Psiold)!=NULL)
							*(IfacePtr->ZSlab->Psiold+k*(NI-1)+i)=*(IfacePtr->ZSlab->Psi+k*(NI-1)+i);//raj
						

						
						//Added by sushil for multiphase flow
						
						
						if((zone->Flowsol->Momentumx)!=NULL && (IfacePtr->ZSlab->MomentumxOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumx+k*(NI-1)+i);
							
						if((zone->Flowsol->Momentumy)!=NULL && (IfacePtr->ZSlab->MomentumyOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumyOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumy+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (IfacePtr->ZSlab->MomentumzOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumz+k*(NI-1)+i);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (IfacePtr->ZSlab->MomentumxstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxstarOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumxstar+k*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (IfacePtr->ZSlab->MomentumystarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumystarOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumystar+k*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (IfacePtr->ZSlab->MomentumzstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzstarOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumzstar+k*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasexOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasex[nophase]+k*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphaseyOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasey[nophase]+k*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasezOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasez[nophase]+k*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (IfacePtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->VarthetaPhaseOld[nophase]+k*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaPhase[nophase]+k*(NI-1)+i);
						}	
						
						//end	
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (IfacePtr->ZSlab->U_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->U_DisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->U_Disperse+k*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (IfacePtr->ZSlab->V_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->V_DisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->V_Disperse+k*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (IfacePtr->ZSlab->W_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->W_DisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->W_Disperse+k*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (IfacePtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(IfacePtr->ZSlab->VarthetaDisperseOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaDisperse+k*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (IfacePtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(IfacePtr->ZSlab->DropFuelMassFracOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->DropFuelMassFrac+k*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (IfacePtr->ZSlab->TdropOld)!=NULL)
							*(IfacePtr->ZSlab->TdropOld+k*(NI-1)+i)=*(IfacePtr->ZSlab->Tdrop+k*(NI-1)+i);

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 5:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->uOld)!=NULL)
							*(IfacePtr->ZSlab->uOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->u+j*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vOld)!=NULL)
							*(IfacePtr->ZSlab->vOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->v+j*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wOld)!=NULL)
							*(IfacePtr->ZSlab->wOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->w+j*(NI-1)+i);
						
						if((zone->Flowsol->T)!=NULL && (IfacePtr->ZSlab->TOld)!=NULL)
							*(IfacePtr->ZSlab->TOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->T+j*(NI-1)+i);
							
						if((zone->Flowsol->enthalpy)!=NULL && (IfacePtr->ZSlab->enthalpyOld)!=NULL)
							*(IfacePtr->ZSlab->enthalpyOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->enthalpy+j*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->umass+j*(NI-1)+i);
							
						if((zone->Flowsol->vmass)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->vmass+j*(NI-1)+i);	
							
						if((zone->Flowsol->wmass)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
						*(IfacePtr->ZSlab->wmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->wmass+j*(NI-1)+i);

					/*
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->uOld+j*(NI-1)+i);
							
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->vOld+j*(NI-1)+i);	
							
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
							*(IfacePtr->ZSlab->wmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->wOld+j*(NI-1)+i);
					*/

						if((zone->Flowsol->TotalMfieldX)!=NULL && (IfacePtr->ZSlab->TotalMfieldXold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldXold+j*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldX+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (IfacePtr->ZSlab->TotalMfieldYold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldYold+j*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldY+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (IfacePtr->ZSlab->TotalMfieldZold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldZold+j*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldZ+j*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (IfacePtr->ZSlab->PotentialMHDold)!=NULL)
							*(IfacePtr->ZSlab->PotentialMHDold+j*(NI-1)+i)=*(IfacePtr->ZSlab->PotentialMHD+j*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (IfacePtr->ZSlab->Psiold)!=NULL)
							*(IfacePtr->ZSlab->Psiold+j*(NI-1)+i)=*(IfacePtr->ZSlab->Psi+j*(NI-1)+i);//raj

						
						//Added by sushil for multiphase flow
						if((zone->Flowsol->Momentumx)!=NULL && (IfacePtr->ZSlab->MomentumxOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumx+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumy)!=NULL && (IfacePtr->ZSlab->MomentumyOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumyOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumy+j*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (IfacePtr->ZSlab->MomentumzOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumz+j*(NI-1)+i);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (IfacePtr->ZSlab->MomentumxstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxstarOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumxstar+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (IfacePtr->ZSlab->MomentumystarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumystarOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumystar+j*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (IfacePtr->ZSlab->MomentumzstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzstarOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumzstar+j*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasexOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasex[nophase]+j*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphaseyOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasey[nophase]+j*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasezOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasez[nophase]+j*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (IfacePtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->VarthetaPhaseOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaPhase[nophase]+j*(NI-1)+i);
						}	
						//end
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (IfacePtr->ZSlab->U_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->U_DisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->U_Disperse+j*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (IfacePtr->ZSlab->V_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->V_DisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->V_Disperse+j*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (IfacePtr->ZSlab->W_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->W_DisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->W_Disperse+j*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (IfacePtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(IfacePtr->ZSlab->VarthetaDisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaDisperse+j*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (IfacePtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(IfacePtr->ZSlab->DropFuelMassFracOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->DropFuelMassFrac+j*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (IfacePtr->ZSlab->TdropOld)!=NULL)
							*(IfacePtr->ZSlab->TdropOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Tdrop+j*(NI-1)+i);


						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
				case 6:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->uOld)!=NULL)
							*(IfacePtr->ZSlab->uOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->u+j*(NI-1)+i);
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vOld)!=NULL)
							*(IfacePtr->ZSlab->vOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->v+j*(NI-1)+i);
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wOld)!=NULL)
							*(IfacePtr->ZSlab->wOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->w+j*(NI-1)+i);
						
						if((zone->Flowsol->T)!=NULL && (IfacePtr->ZSlab->TOld)!=NULL)
							*(IfacePtr->ZSlab->TOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->T+j*(NI-1)+i);
							
						if((zone->Flowsol->enthalpy)!=NULL && (IfacePtr->ZSlab->enthalpyOld)!=NULL)
							*(IfacePtr->ZSlab->enthalpyOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->enthalpy+j*(NI-1)+i);
						
						if((zone->Flowsol->umass)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->umass+j*(NI-1)+i);
							
						if((zone->Flowsol->vmass)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->vmass+j*(NI-1)+i);	
							
						if((zone->Flowsol->wmass)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
						*(IfacePtr->ZSlab->wmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->wmass+j*(NI-1)+i);		

					/*
						if((zone->Flowsol->u)!=NULL && (IfacePtr->ZSlab->umassOld)!=NULL)
							*(IfacePtr->ZSlab->umassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->uOld+j*(NI-1)+i);
							
						if((zone->Flowsol->v)!=NULL && (IfacePtr->ZSlab->vmassOld)!=NULL)
							*(IfacePtr->ZSlab->vmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->vOld+j*(NI-1)+i);	
							
						if((zone->Flowsol->w)!=NULL && (IfacePtr->ZSlab->wmassOld)!=NULL)
							*(IfacePtr->ZSlab->wmassOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->wOld+j*(NI-1)+i);
					*/

						if((zone->Flowsol->TotalMfieldX)!=NULL && (IfacePtr->ZSlab->TotalMfieldXold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldXold+j*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldX+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldY)!=NULL && (IfacePtr->ZSlab->TotalMfieldYold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldYold+j*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldY+j*(NI-1)+i);//raj

						if((zone->Flowsol->TotalMfieldZ)!=NULL && (IfacePtr->ZSlab->TotalMfieldZold)!=NULL)
							*(IfacePtr->ZSlab->TotalMfieldZold+j*(NI-1)+i)=*(IfacePtr->ZSlab->TotalMfieldZ+j*(NI-1)+i);//raj
						if((zone->Flowsol->PotentialMHD)!=NULL && (IfacePtr->ZSlab->PotentialMHDold)!=NULL)
							*(IfacePtr->ZSlab->PotentialMHDold+j*(NI-1)+i)=*(IfacePtr->ZSlab->PotentialMHD+j*(NI-1)+i);//raj
						if((zone->Flowsol->Psi)!=NULL && (IfacePtr->ZSlab->Psiold)!=NULL)
							*(IfacePtr->ZSlab->Psiold+j*(NI-1)+i)=*(IfacePtr->ZSlab->Psi+j*(NI-1)+i);//raj

						//Added by sushil for multipahse flow
						if((zone->Flowsol->Momentumx)!=NULL && (IfacePtr->ZSlab->MomentumxOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumx+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumy)!=NULL && (IfacePtr->ZSlab->MomentumyOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumyOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumy+j*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumz)!=NULL && (IfacePtr->ZSlab->MomentumzOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumz+j*(NI-1)+i);
						
						if((zone->Flowsol->Momentumxstar)!=NULL && (IfacePtr->ZSlab->MomentumxstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumxstarOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumxstar+j*(NI-1)+i);
							
						if((zone->Flowsol->Momentumystar)!=NULL && (IfacePtr->ZSlab->MomentumystarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumystarOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumystar+j*(NI-1)+i);	
						
						if((zone->Flowsol->Momentumzstar)!=NULL && (IfacePtr->ZSlab->MomentumzstarOld)!=NULL)
							*(IfacePtr->ZSlab->MomentumzstarOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumzstar+j*(NI-1)+i);
							
						for(nophase=0;nophase<nop;nophase++)
						{
							if((zone->Flowsol->Momentumofphasex[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasexOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasexOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasex[nophase]+j*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasey[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphaseyOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphaseyOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasey[nophase]+j*(NI-1)+i);
							
							
							if((zone->Flowsol->Momentumofphasez[nophase])!=NULL && (IfacePtr->ZSlab->MomentumofphasezOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->MomentumofphasezOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Momentumofphasez[nophase]+j*(NI-1)+i);
							
							if((zone->Flowsol->VarthetaPhase[nophase])!=NULL && (IfacePtr->ZSlab->VarthetaPhaseOld[nophase])!=NULL)
							
								*(IfacePtr->ZSlab->VarthetaPhaseOld[nophase]+j*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaPhase[nophase]+j*(NI-1)+i);
						}	
						
						//end
						//Added by adhi+shash for disperse phase flow //changed by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (IfacePtr->ZSlab->U_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->U_DisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->U_Disperse+j*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (IfacePtr->ZSlab->V_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->V_DisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->V_Disperse+j*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (IfacePtr->ZSlab->W_DisperseOld)!=NULL)
							*(IfacePtr->ZSlab->W_DisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->W_Disperse+j*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (IfacePtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(IfacePtr->ZSlab->VarthetaDisperseOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->VarthetaDisperse+j*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (IfacePtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(IfacePtr->ZSlab->DropFuelMassFracOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->DropFuelMassFrac+j*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (IfacePtr->ZSlab->TdropOld)!=NULL)
							*(IfacePtr->ZSlab->TdropOld+j*(NI-1)+i)=*(IfacePtr->ZSlab->Tdrop+j*(NI-1)+i);

						for(iScal=0;iScal<10;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
			}
		}
	}
	
	return 0;
}




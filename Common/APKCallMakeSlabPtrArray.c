#include "APKCallMakeSlabPtrArray.h"

/**
This function is called in Main() just before calling the solvers to fix the Interface and boundary slabs i.e., to assign
respective pointers to variables SlabUmass, SlabVmass, SlabWmass, SlabU, SlabV, SlabW, SlabPress, SlabT and SlabScalar
which will be further directly used during various Gauss-Seidel iterations like UpdateVar, findresidual and UpdatePressure.
1. This function has to be taken care while adding any new variable to the solver
*/

int CallMakeSlabPtrArray(Base *base)
{
	Zone *zone;
	int nZones,numZones,nop=0,iScal,nophase;
	
	nZones=(base->Number_of_zones);
	
	nop=base->nphase;
	
	
	
	
	if(base->BaseParameters->flag2d==1) //2D
	{
		for(numZones=0;numZones<nZones;numZones++)
		{	zone = (base->CgnsZone+numZones);
			if(base->flagMPFV==0)
			{
				
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabUmass,100);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVmass,200);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabPress,4);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabU,1);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabV,2);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVol,5);
					if(base->CgnsZone->Flowsol->T!=NULL)///This statement needs to be checked for combustion module. This looks bit odd.  Bhanesh.
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabT,6);
				if((zone->Flowsol)->flagRho==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDensity,710);
				
				if((zone->Flowsol)->flagViscoeff==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Visc,711);
				
				if((zone->Flowsol)->flagDiffcoeff==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_K,712);
				if((zone->Flowsol)->flagSpecificHeat==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_CP,713);
				
				if(base->flagTurb==1)
				{
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Eff_Visc,714);
					
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Tvisc,715);
				}
				//double *Slab_EffScal[6];
				//double *Slab_TurbCond[6];
				
				
				//############################################################################
				//adhi+shash
				if((base->flagDroplet==1)||(base->flagParticle==1))//added by Adhiraj for
				{                                                  //Disperse phase 
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVxDisperse,62);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVyDisperse,63);
				//MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVzDisperse,64);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVarthetaDisperse,65);
				MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVarthetaDisperse,650);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabRhoDisperse,66);

				MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVxDisperse,620);	
				MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVyDisperse,630);

				//Added by chandrasekhar
				if(base->flagTurb==1)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDummyViscosity_disperse,726);
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDiffVolfrac_disperse,727);
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVxDisperse,620);	
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVyDisperse,630);
					}
				//End	
				//the following lines are commented as the number density equation is no longer solved for in gas-droplet flows module
				/*if(base->flagDroplet==1)
				{
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabNumDensity,67);
				}
				*/

				if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
				{
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabUDmass,6200);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVDmass,6300);
						
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabTdrop,69);
					MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabTdrop,690);
				}
				if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1)&&(base->BaseParameters->flagDropFuelMassFrac==1))
	                        {
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDropFuelMassFrac,68);
					MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabDropFuelMassFrac,680);

				}
				}
				//############################################################################
				//	if(base->CgnsZone->Flowsol->T!=NULL && base->flagTurbComb==0)///This statement needs to be checked for combustion module. This looks bit odd.  Bhanesh.
				//	MakeSlabPtrArray_com(zone,zone->Flowsol->SlabT,6);
				
				if(base->CgnsZone->Flowsol->enthalpy!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Enthalpy,7);
				
				if(base->flagSolidification==1 && base->CgnsZone->Flowsol->Latent_Heat_Content!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Latent_Heat_Content,21);
	
				if(base->CgnsZone->Flowsol->Potential!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Potential,8);
				if(base->CgnsZone->Flowsol->RhoIon!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_RhoIon,9);
	
				if(base->CgnsZone->Flowsol->MixFrac!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabMixFrac,14);
	
				if(base->CgnsZone->Flowsol->FuelMassFrac!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabFuelFrac,15);
				
				if(base->CgnsZone->Flowsol->PotentialMHD!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Potential_MHD,17);
				if(base->CgnsZone->Flowsol->TotalMfieldX!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Total_MfieldX,56);
				if(base->CgnsZone->Flowsol->TotalMfieldY!=NULL)//rajendra
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Total_MfieldY,57);
				if(base->CgnsZone->Flowsol->TotalMfieldZ!=NULL)//rajendra
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Total_MfieldZ,58);
				if(base->CgnsZone->Flowsol->Psi!=NULL)//rajendra
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Psi_MHD,59);
				if(base->CgnsZone->Flowsol->Flow_Potential!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Flow_Potential,990);
				if(base->flagGroundwater==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_H,19);
				
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabU_old,10);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabV_old,20);
			}
				//Added by sushil....For multiphase flow variables......
			
			else
			{
				
				
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVol,5);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDensity,710);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabPress,4);
				if(base->CgnsZone->Flowsol->Momentumx!=NULL)
					{	
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVx,31);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabVxstar),3100);
					
					
					}
				if(base->CgnsZone->Flowsol->Momentumy!=NULL)
					{	
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabVy),32);
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabVystar),3200);
					}
				
				//Added by chandrasekhar --- To make slabs of Turbulent viscosity and effective viscosity in Flowsol datastructure
				if(base->flagTurb==1)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Eff_Visc,714);
					
						MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Tvisc,715);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[0]),730);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[1]),731);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[2]),732);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[3]),733);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[4]),734);
						
					}
				//End --- chandrasekhar
				
				for(nophase=0;nophase<nop;nophase++)
				{
						if(base->CgnsZone->Flowsol->VarthetaPhase[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVartheta[nophase],34+nophase);
					
						if(base->CgnsZone->Flowsol->Momentumofphasex[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabVxPhase[nophase]),41+3*nophase);
							
						if(base->CgnsZone->Flowsol->Momentumofphasey[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabVyPhase[nophase]),41+3*nophase+1);
				}
					
				//Added by chandrasekhar--To compute the gradients in case of turbulent MPFV	
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVx_Old,310);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVy_Old,320);
				
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVcx_Old,410);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVcy_Old,420);
				//End --- chandrasekhar
			}	//End of multiphase flow varibles.........
			for(iScal=0;iScal<10;iScal++)
			{
				if(base->CgnsZone->Flowsol->Scal[iScal]!=NULL)
				{
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabScalar[iScal],600+iScal);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabScalar_old[iScal],(600+iScal)*10);
				}
			}
		}
	}
	else if(base->BaseParameters->flag2d==0) //3D
	/**This condition can be removed by inserting another 'if' condition in upper 'if' loop**/
	{
		for(numZones=0;numZones<nZones;numZones++)
		{	if(base->flagMPFV==0)
			{
				zone = (base->CgnsZone+numZones);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabUmass,100);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVmass,200);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabWmass,300);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabPress,4);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabU,1);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabV,2);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabW,3);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVol,5);
				
				if((zone->Flowsol)->flagRho==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDensity,710);
				
				if((zone->Flowsol)->flagViscoeff==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Visc,711);
				
				if((zone->Flowsol)->flagDiffcoeff==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_K,712);
				if((zone->Flowsol)->flagSpecificHeat==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_CP,713);
				
				if(base->flagTurb==1)
				{
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Eff_Visc,714);
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Tvisc,715);
				}
				
				//############################################################################
				//adhi+shash
				if((base->flagDroplet==1)||(base->flagParticle==1))//added by Adhiraj for
				{                                                  //Disperse phase 
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVxDisperse,62);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVyDisperse,63);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVzDisperse,64);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVarthetaDisperse,65);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabRhoDisperse,66);

					MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVxDisperse,620);
					MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVyDisperse,630);
					MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVzDisperse,640);	

					//chandrasekhar
					if(base->flagTurb==1)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDummyViscosity_disperse,726);
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDiffVolfrac_disperse,727);
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVxDisperse,620);	
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVyDisperse,630);
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVzDisperse,640);					
					}
					//End
					MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVarthetaDisperse,650);
					//the following lines are commented as the number density equation is no longer solved for in gas-droplet flows module
					/*if(base->flagDroplet==1)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabNumDensity,67);
					}*/
					if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
	                                {
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabUDmass,6200);	
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabVDmass,6300);
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabWDmass,6400);	
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabTdrop,69);
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabTdrop,690);
					}
					if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1)&&(base->BaseParameters->flagDropFuelMassFrac==1))
	                                {
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDropFuelMassFrac,68);
						MakeSlabPtrArray_com(zone,zone->PrevFlowsol->SlabDropFuelMassFrac,680);
					}
					
				}
				//############################################################################
				if(base->CgnsZone->Flowsol->T!=NULL)///IF condition is not identical for temperature for 2d and 3d solvers.  This has to be checked.  Bhanesh.
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabT,6);
				if(base->CgnsZone->Flowsol->enthalpy!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Enthalpy,7);
				
				if(base->flagSolidification==1 && base->CgnsZone->Flowsol->Latent_Heat_Content!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Latent_Heat_Content,21);
	
	
				if(base->CgnsZone->Flowsol->Potential!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Potential,8);
				if(base->CgnsZone->Flowsol->RhoIon!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_RhoIon,9);
				if(base->CgnsZone->Flowsol->MixFrac!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabMixFrac,14);
				
				if(base->CgnsZone->Flowsol->FuelMassFrac!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabFuelFrac,15);
				
				
				if(base->CgnsZone->Flowsol->PotentialMHD!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Potential_MHD,17);
				
				if(base->CgnsZone->Flowsol->TotalMfieldX!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Total_MfieldX,56);
				if(base->CgnsZone->Flowsol->TotalMfieldY!=NULL)//rajendra
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Total_MfieldY,57);
				if(base->CgnsZone->Flowsol->TotalMfieldZ!=NULL)//rajendra
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Total_MfieldZ,58);
				if(base->CgnsZone->Flowsol->Psi!=NULL)//rajendra
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Psi_MHD,59);
				
				if(base->CgnsZone->Flowsol->Flow_Potential!=NULL)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Flow_Potential,990);
				
				if(base->flagGroundwater==1)
					MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_H,19);
				
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabU_old,10);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabV_old,20);
					MakeSlabPtrArray_com(zone,zone->Flowsol->SlabW_old,30);

				//Added by sushil....For multiphase flow variables......
			}
			else
			{	
				zone = (base->CgnsZone+numZones);//added --chandrasekhar
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVol,5);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabDensity,710);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabPress,4);
				if(base->CgnsZone->Flowsol->Momentumx!=NULL)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVx,31);
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVxstar,3100);
					}
				if(base->CgnsZone->Flowsol->Momentumy!=NULL)
					{	
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVy,32);
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVystar,3200);
					}
				if(base->CgnsZone->Flowsol->Momentumz!=NULL)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVz,33);
						MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVzstar,3300);
					}
				//Added by chandrasekhar --- To make slabs of Turbulent viscosity and effective viscosity in Flowsol datastructure
				if(base->flagTurb==1)
					{
						MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Eff_Visc,714);
					
						MakeSlabPtrArray_com(zone,zone->Flowsol->Slab_Tvisc,715);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[0]),730);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[1]),731);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[2]),732);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[3]),733);
						
						MakeSlabPtrArray_com(zone,(zone->Flowsol->SlabDummyeffViscoeffofphase[4]),734);
						
					}
				//End --- chandrasekhar
				for(nophase=0;nophase<nop;nophase++)
				{
						
						if(base->CgnsZone->Flowsol->VarthetaPhase[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVartheta[nophase],34+nophase);
						
						if(base->CgnsZone->Flowsol->Momentumofphasex[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVxPhase[nophase],41+3*nophase);
							
						if(base->CgnsZone->Flowsol->Momentumofphasey[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVyPhase[nophase],41+3*nophase+1);
						
						if(base->CgnsZone->Flowsol->Momentumofphasez[nophase]!=NULL)
							MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVzPhase[nophase],41+3*nophase+2);
				}
				//Added by chandrasekhar--To compute the gradients in case of turbulent MPFV	
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVx_Old,310);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVy_Old,320);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVz_Old,330);
				
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVcx_Old,410);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVcy_Old,420);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabVcz_Old,430);
				//End --- chandrasekhar	
			}	//End of multiphase flow varibles.........
			
			
			
			for(iScal=0;iScal<10;iScal++)
			{
				if(base->CgnsZone->Flowsol->Scal[iScal]!=NULL)
				{
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabScalar[iScal],600+iScal);
				MakeSlabPtrArray_com(zone,zone->Flowsol->SlabScalar_old[iScal],(600+iScal)*10);
				}
			}
		}
	}
	

	return(0);
}

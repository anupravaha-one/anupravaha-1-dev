#include "APKRefreshVariableCoefficients.h"

/**
This function is similar to "RefresSurfaceCells". but the only difference is this is used for 
coefficients which are varying like scalar coefficients, variable specific heats, densities
viscocities etc.
*/

int RefreshVariableCoefficients(int variableTag,Base *base)
{
	int cntZone,numZones;
	Zone *zone;
	
	numZones=base->Number_of_zones;
	
	//loop for variable Scalar Diffusion coefficients
	if(variableTag >= 700 && variableTag <= 709)
	{
		if((variableTag-700) >= base->BaseParameters->Turb->numTurbScalars)
		{
			if(base->flagGroundwater==0)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
				{
					zone = base->CgnsZone+cntZone;
					
					if(zone->Flowsol->flagScalcoeff[variableTag-700]==1)
					{
						ExtractVariableCoefficients(variableTag,base,zone);
					}
				}
			}
			else
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
				{
					zone = base->CgnsZone+cntZone;
					
					//printf("zone->Difftag=%d\n", zone->Difftag);
					if(zone->Difftag!=0)
					{
						//if(zone->Flowsol->flagScalcoeff[variableTag-700]==1) <<< This should be activated, Check before doing so <<< Narsimha 18th April, 2008.
							Extract_Var_Coeff_PorousMedia(variableTag,base,zone);
					}
				}
			}
		}
	}
	else if(variableTag ==712)  /*It T is not NULL then refresh Thermal conductivity*/
	{
		if(base->CgnsZone->Flowsol->Diffcoeff!=NULL)
		{
			if(base->CgnsZone->Flowsol->flagDiffcoeff==1)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if(variableTag ==713)  /*It T is not NULL then refresh Specific Heat*/
	{
		if(base->CgnsZone->Flowsol->SpecificHeat!=NULL)
		{
			if(base->CgnsZone->Flowsol->flagSpecificHeat==1)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if(variableTag == 710)           /*This case is for Density*/
	{
		if(((base->CgnsZone->Flowsol)->flagRho)==1)
		{
			for(cntZone=0;cntZone < numZones;cntZone++)
				ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
		}
	}
	else if(variableTag == 711)           /*This case is for Viscosity*/
	{
		if(((base->CgnsZone->Flowsol)->flagViscoeff)==1)
		{
			for(cntZone=0;cntZone < numZones;cntZone++)
				ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
		}
	}
	else if(variableTag ==717)  /*It T is not NULL then refresh Specific Heat*/
	{
		if(base->CgnsZone->Flowsol->ChargeDiffcoeff!=NULL)
		{
			if(base->CgnsZone->Flowsol->flagChargeDiffCoeff==1)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if((variableTag == 719)||(variableTag>=1&&variableTag<=3))
	{
		if(base->CgnsZone->Flowsol->UnSatPermeability!=NULL)//shash this condition has to be imposed
		{
		for(cntZone=0;cntZone < numZones;cntZone++)
			  Extract_Var_Coeff_PorousMedia(variableTag,base,base->CgnsZone+cntZone);
		}
	}
	
	//for multiphase flow added by sushil
	else if((variableTag == 720)&&(base->flagMPFV == 1))           /*This case is for Viscosity*/
	{	
		if(((base->CgnsZone->Flowsol)->flagViscoeff)==1)
		{
			if(base->CgnsZone->Flowsol->ViscoeffDummy[0]!=NULL)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if((variableTag == 721)&&(base->flagMPFV == 1))           /*This case is for Viscosity*/
	{
		if(((base->CgnsZone->Flowsol)->flagViscoeff)==1)
		{
			if(base->CgnsZone->Flowsol->ViscoeffDummy[1]!=NULL)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if((variableTag == 722)&&(base->flagMPFV == 1))           /*This case is for Viscosity*/
	{
		if(((base->CgnsZone->Flowsol)->flagViscoeff)==1)
		{
			if(base->CgnsZone->Flowsol->ViscoeffDummy[2]!=NULL)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if((variableTag == 723)&&(base->flagMPFV == 1))           /*This case is for Viscosity*/
	{
		if(((base->CgnsZone->Flowsol)->flagViscoeff)==1)
		{
			if(base->CgnsZone->Flowsol->ViscoeffDummy[3]!=NULL)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if((variableTag == 724)&&(base->flagMPFV == 1))           /*This case is for Viscosity*/
	{
		if(((base->CgnsZone->Flowsol)->flagViscoeff)==1)
		{
			if(base->CgnsZone->Flowsol->ViscoeffDummy[4]!=NULL)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	else if(variableTag == 725)  /*This case is for Cp of droplet temperature equation & 726 may be used for mass diffusivity of mass fraction equation*/
	{
		if(base->CgnsZone->Flowsol->CPdrop!=NULL)
		{
			if(base->CgnsZone->Flowsol->flagCPdrop==1)
			{
				for(cntZone=0;cntZone < numZones;cntZone++)
					ExtractVariableCoefficients(variableTag,base,base->CgnsZone+cntZone);
			}
		}
	}
	
	return 0;
}

		

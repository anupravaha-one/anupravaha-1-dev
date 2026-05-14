/***************************************************************************************************

                                     WRITTEN BY BHANESH edited by anji
				     
****************************************************************************************************/				     
#include "APKinitialize_TwoPhase.h"

/**
This functions is called for both solidification-melting and Flow boiling-HEM problems, which deal with Two phase and 
phase change. Although it can be observed that initialization for Flow boiling problem with HEM (flag twophase) and
Solidification is different.
*/
int initialize_TwoPhase(Base *base,int *count)
{
/*                              vairable 'count' is not incremented in this routine
                                as no other solvable variable is initiated.........*/

	double density,k,constant=1.0,viscosity,E_fm,E_sm,C_p;
	int j,num_zone,numiface,numBc,iface,zBc,dir,i,NI,NJ,NK,total,*size;
	Zone *z,*zone;
	num_zone= base->Number_of_zones;
	z= base->CgnsZone;
	if ((base->flagSolidification==1)||(base->flagTwoPhase==1))
	{
		for(j=0;j<num_zone;j++)
		{       
			/*Memory allocation in the Flow Solution Node X Direction Velocity*/
			if (((z+j)->Flowsol)->flagRho!=1)
			{
				density = *(((z+j)->Flowsol)->Rho);
				((z+j)->Flowsol)->Rho= (double *)initializeVarCoefficient(1,(z+j),density,1);
				((z+j)->PrevFlowsol)->Rho= (double *)initializeVarCoefficient(1,(z+j),density,1);
				if(base->flagSolidification==1)
				((z+j)->PrevFlowsol)->Rho_prev= (double *)initializeVarCoefficient(1,(z+j),density,1);
				((z+j)->Flowsol)->flagRho=1;
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
						if(((z+j)->Iface+iface)->Direction==dir)
						(((z+j)->Iface+iface)->ZSlab)->Rho=(double *)initializeSlabVarCoeff(1,(z+j),density,1,dir,((z+j)->Iface+iface)->Point_Range);
					for(zBc=0;zBc<numBc;zBc++)
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						(((z+j)->ZoneBc+zBc)->ZSlab)->Rho = (double *)initializeSlabVarCoeff(1,(z+j),density,1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
		for(j=0;j<num_zone;j++)
		{       
			//Memory allocation in the Flow Solution Node X Direction Velocity
			if (((z+j)->Flowsol)->flagViscoeff!=1)
			{
				viscosity = *(((z+j)->Flowsol)->Viscoeff);
				base->BaseParameters->Mu_f = viscosity;
				base->BaseParameters->Mu_s = viscosity;
				((z+j)->Flowsol)->Viscoeff= (double *)initializeVarCoefficient(1,(z+j),viscosity,1);
				((z+j)->PrevFlowsol)->Viscoeff= (double *)initializeVarCoefficient(1,(z+j),viscosity,1);
				((z+j)->Flowsol)->flagViscoeff=1;
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
						if(((z+j)->Iface+iface)->Direction==dir)
						(((z+j)->Iface+iface)->ZSlab)->Viscoeff=(double *)initializeSlabVarCoeff(1,(z+j),viscosity,1,dir,((z+j)->Iface+iface)->Point_Range);
					for(zBc=0;zBc<numBc;zBc++)
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						(((z+j)->ZoneBc+zBc)->ZSlab)->Viscoeff =(double *)initializeSlabVarCoeff(1,(z+j),viscosity,1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
			printf("Viscosity entered =%lf\n",viscosity);
		}
		//print_enthalpy(base);
		//exit(0);
		for(j=0;j<num_zone;j++)
		{       
			if (((z+j)->Flowsol)->flagDiffcoeff!=1)
			{
				k = *(((z+j)->Flowsol)->Diffcoeff);
				//printf("k address=%ld\n",((z+j)->Flowsol)->Diffcoeff);
				((z+j)->Flowsol)->Diffcoeff= (double *)initializeVarCoefficient(1,(z+j),k,1);
				((z+j)->PrevFlowsol)->Diffcoeff= (double *)initializeVarCoefficient(1,(z+j),k,1);
				((z+j)->Flowsol)->flagDiffcoeff=1;
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
						if(((z+j)->Iface+iface)->Direction==dir)
						(((z+j)->Iface+iface)->ZSlab)->Diffcoeff=(double *)initializeSlabVarCoeff(1,(z+j),k,1,dir,((z+j)->Iface+iface)->Point_Range);
					for(zBc=0;zBc<numBc;zBc++)
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						(((z+j)->ZoneBc+zBc)->ZSlab)->Diffcoeff = (double *)initializeSlabVarCoeff(1,(z+j),k,1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
		for(j=0;j<num_zone;j++)
		{       
			if (((z+j)->Flowsol)->flagSpecificHeat!=1)
			{
				C_p = *(((z+j)->Flowsol)->SpecificHeat);
				
				((z+j)->Flowsol)->SpecificHeat= (double *)initializeVarCoefficient(1,(z+j),C_p,1);
				((z+j)->PrevFlowsol)->SpecificHeat= (double *)initializeVarCoefficient(1,(z+j),C_p,1);
				((z+j)->Flowsol)->flagSpecificHeat=1;
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
						if(((z+j)->Iface+iface)->Direction==dir)
						(((z+j)->Iface+iface)->ZSlab)->SpecificHeat=(double *)initializeSlabVarCoeff(1,(z+j),C_p,1,dir,((z+j)->Iface+iface)->Point_Range);
					for(zBc=0;zBc<numBc;zBc++)
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						(((z+j)->ZoneBc+zBc)->ZSlab)->SpecificHeat = (double *)initializeSlabVarCoeff(1,(z+j),C_p,1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
		for(j=0;j<num_zone;j++)
		{       
			((z+j)->Flowsol)->f= (double *)initialize_variable(1,(z+j),constant);
			//printf("f address=%ld\n",((z+j)->Flowsol)->f);
			//((z+j)->PrevFlowsol)->f= (double *)initializeVarCoefficient(1,(z+j),1.0,1);
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
					if(((z+j)->Iface+iface)->Direction==dir)
					(((z+j)->Iface+iface)->ZSlab)->f=(double *)initializeslab(1,dir,((z+j)->Iface+iface)->Point_Range);
				for(zBc=0;zBc<numBc;zBc++)
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					(((z+j)->ZoneBc+zBc)->ZSlab)->f = (double *)initializeslab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
			}
		}
		
		if (base->flagSolidification==1)
		{
			printf("\nIs the problem is mushy region problem? (1. yes 0. No):  ");
			scanf ("%d",&(base->BaseParameters->flag_mushy));
			if ((base->BaseParameters->flag_mushy)==1)
			{
				printf("\nEnter Solidus Temperature: ");
				scanf ("%lf",&(base->BaseParameters->Solidus_Temp));
				printf("\nEnter Liquidus Temperature: ");
				scanf ("%lf",&(base->BaseParameters->Liquidus_Temp));
			}
			else
			{
			
				printf("Enter the Melting Temperature of the solid: ");
				scanf ("%lf",&(base->BaseParameters->Melt_Temp));
			}
			printf("\nVariation of C_p, K are assumed to be varying linearly in the melting region\n");
			printf("\nEnter the Specific heat value in the solid  region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->CP_s));
			printf("\nEnter the Specific heat value in the liquid region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->CP_f));
			printf("\nEnter the Density value in the solid  region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->Rho_s));
			printf("\nEnter the Density value in the liquid region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->Rho_f));
			printf("\nEnter the Thermal Conductivity value in the solid  region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->K_s));
			printf("\nEnter the Thermal Conductivity value in the liquid region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->K_f));
			printf("\nEnter the Enthalpy value in the solid  region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->E_sm));
			printf("\nEnter the Enthalpy value in the liquid region at Melting Temp: ");
			scanf ("%lf",&(base->BaseParameters->E_fm));
			
			printf("\nEnter the value of 'C' in Carman Kozney equation for flow in porous media");
			printf("\nLower value indicated higher permeability of mush in case of mushy region phase change problems\n");
			scanf ("%lf",&(base->BaseParameters->C));
			
			printf("\nEnthalpy values scanned are %e %e\n",(base->BaseParameters->E_sm),(base->BaseParameters->E_fm));
			
		
		//Created by bhans Algo 2
		for(j=0;j<num_zone;j++)
		{       
			((z+j)->Flowsol)->Latent_Heat_Content    = (double *)initialize_variable(1,(z+j),constant);
			((z+j)->PrevFlowsol)->Latent_Heat_Content= (double *)initialize_variable(1,(z+j),constant);
			zone = (base->CgnsZone+j);
			size=(base->CgnsZone+j)->Size;
			NI = *(size);
			NJ = *(size+1);
			NK = *(size+2);
			total = (NI-1)*(NJ-1)*(NK-1);
			E_sm=(base->BaseParameters->E_sm);
			E_fm=(base->BaseParameters->E_fm);
			for (i=0;i<total;i++)
			{
				if (*((zone->Flowsol)->enthalpy+i)<=E_sm)
				*((zone->Flowsol)->Latent_Heat_Content+i) = 0.0;
				else
				*((zone->Flowsol)->Latent_Heat_Content+i) = E_fm - E_sm;
				//printf("%e\t",*((zone->Flowsol)->Latent_Heat_Content+i));
			}
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
					if(((z+j)->Iface+iface)->Direction==dir)
					(((z+j)->Iface+iface)->ZSlab)->Latent_Heat_Content=(double *)initializeslab(1,dir,((z+j)->Iface+iface)->Point_Range);
				for(zBc=0;zBc<numBc;zBc++)
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					(((z+j)->ZoneBc+zBc)->ZSlab)->Latent_Heat_Content = (double *)initializeslab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
			}
		}
		initialize_Latent_Heat(base);
		//exit(0);
		//Solidification_updatevariable(base); This is hashed because of folder
		//priorities that has to be provided.  This is called directly in main.
		}
		else if(base->flagTwoPhase==1)
		{		
				printf("\nEnter the  reference pressure for the Saturation  Temparature: ");
				scanf ("%lf",&(base->BaseParameters->P_o));
				printf("\nEnter the Density value in the vapor  region at Saturation  Temparature: ");
				scanf ("%lf",&(base->BaseParameters->Rho_s));
				printf("\nEnter the Density value in the liquid region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->Rho_f));
				printf("\nEnter the Viscosity value in the vapor  region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->Mu_s));
				printf("\nEnter the Viscosity value in the liquid region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->Mu_f));
				printf("\nEnter the Thermal Conductivity value in the vapor  region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->K_s));
				printf("\nEnter the Thermal Conductivity value in the liquid region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->K_f));
				printf("\nEnter the Specific Heat value in the vapor  region at Saturation  Temparature: ");
				scanf ("%lf",&(base->BaseParameters->CP_s));
				printf("\nEnter the Specific Heat value in the Liquid region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->CP_f));
				printf("\nEnter the Internal Energy value in the vapor  region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->E_sm));
				printf("\nEnter the Internal Energy value in the liquid region at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->E_fm));
				printf("\nEnter the del(Rho)/del(P) value for liquid at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->Pc_l));
				printf("\nEnter the del(Rho)/del(P) value for vapor at Saturation Temparature: ");
				scanf ("%lf",&(base->BaseParameters->Pc_v));
				//printf("\n Latent Enthalpy value in the vapor  region at Saturation Temparature: ");
				
				//scanf ("%lf",&(base->BaseParameters->e_L));
				for(j=0;j<num_zone;j++)
				{
					((z+j)->Flowsol)->Latent_Heat_Content= (double *)initialize_variable(0,(z+j),constant);
					//((z+j)->PrevFlowsol)->Temp= (double *)initializeVarCoefficient(1,(z+j),1.0,1);
					numiface=(z+j)->Number_of_Interfaces;
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
							if(((z+j)->Iface+iface)->Direction==dir)
							(((z+j)->Iface+iface)->ZSlab)->Latent_Heat_Content=(double *)initializeslab(0,dir,((z+j)->Iface+iface)->Point_Range);
						for(zBc=0;zBc<numBc;zBc++)
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							(((z+j)->ZoneBc+zBc)->ZSlab)->Latent_Heat_Content = (double *)initializeslab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}

		}
	}
	else
	{	
		for(j=0;j<num_zone;j++)
		{       
			((z+j)->Flowsol)->f= (double *)initialize_variable(0,(z+j),constant);
			((z+j)->PrevFlowsol)->Rho_prev= NULL;
			//((z+j)->PrevFlowsol)->f= (double *)initializeVarCoefficient(1,(z+j),1.0,1);
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
					if(((z+j)->Iface+iface)->Direction==dir)
					(((z+j)->Iface+iface)->ZSlab)->f=(double *)initializeslab(0,dir,((z+j)->Iface+iface)->Point_Range);
				for(zBc=0;zBc<numBc;zBc++)
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					(((z+j)->ZoneBc+zBc)->ZSlab)->f = (double *)initializeslab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range);
			}
		}
		for(j=0;j<num_zone;j++)
		{       
			((z+j)->Flowsol)->Latent_Heat_Content= (double *)initialize_variable(0,(z+j),constant);
			//((z+j)->PrevFlowsol)->Temp= (double *)initializeVarCoefficient(1,(z+j),1.0,1);
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
					if(((z+j)->Iface+iface)->Direction==dir)
					(((z+j)->Iface+iface)->ZSlab)->Latent_Heat_Content=(double *)initializeslab(0,dir,((z+j)->Iface+iface)->Point_Range);
				for(zBc=0;zBc<numBc;zBc++)
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					(((z+j)->ZoneBc+zBc)->ZSlab)->Latent_Heat_Content = (double *)initializeslab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range);
			}
		}
	}
	return 0;
}



int initialize_Latent_Heat(Base *base)
{
	int numzones,i,total,NI,NJ,NK,no,numBc,numiface,index;
	int number,direction;
	Zone *zone;
	int *size,*Point_Range;
	double E_sm,E_fm;
	char tag,name[33];
	void *ptr;
	Interface *facePtr;
	Zone *ptrZoneNeighbor;
	E_sm = base->BaseParameters->E_sm;
	E_fm = base->BaseParameters->E_fm;
	numzones = base->Number_of_zones;
	for (no=0;no<numzones;no++)
	{
		zone = (base->CgnsZone+no);
		size=(base->CgnsZone+no)->Size;
		NI = *(size);
		NJ = *(size+1);
		NK = *(size+2);
		total = (NI-1)*(NJ-1)*(NK-1);
		for (i=0;i<total;i++)
		{
			if(*((zone->Flowsol)->enthalpy+i)>E_sm)
			*((zone->Flowsol)->Latent_Heat_Content+i) = E_fm -E_sm;
			else
			*((zone->Flowsol)->Latent_Heat_Content+i) = 0.0;
			
		}
		numBc=zone->Number_of_BC_Patches;
		for (index=0;index<numBc;index++)
		{
			direction= (zone->ZoneBc+index)->Direction;
			Point_Range=((zone)->ZoneBc+index)->Point_Range;
			if(direction==1||direction==2)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
			if(direction==3||direction==4)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
			if(direction==5||direction==6)
			number=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

			
			for (i=0;i<number;i++)
			{
				if(*(((zone->ZoneBc+index)->ZSlab)->enthalpy+i)>E_sm)
				*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i) = E_fm - E_sm;
				else
				*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i) = 0.0;
			}

		}
		
		numiface=zone->Number_of_Interfaces;
		for (index=0;index<numiface;index++)
		{
			
			direction= (zone->Iface+index)->Direction;
			ptr = (void *)FindInterfaceBCPtrMod(zone,direction,name,&tag);
			if (tag!='i')
			{
				printf("Some thing wrong in update_dele_iteration\n");
				exit(1);
			}
			facePtr=(Interface *)ptr;
			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);
			FillerInterfaceVariables(21,facePtr,ptrZoneNeighbor,direction-1);
		}
	}
	
	return 0;
	
}



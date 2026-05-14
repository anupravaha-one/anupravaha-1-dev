#include "APKInitializeNonSolveCombVariables.h"


int InitializeNonSolveCombVariables(Base *CgnsBase)
{
int i,j,k,num_zone,numiface,numBc,dir,iface,zBc;
double constant = 0.0;
Zone *z;

num_zone = CgnsBase->Number_of_zones;
z = CgnsBase->CgnsZone;

	//Temperature is initialized for all types of reactions
	i=1;
	for(j=0;j<num_zone;j++)
	{
		constant = *(CgnsBase->BaseParameters->InitialTemp);
		((z+j)->Flowsol)->T= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->T= (double *)initialize_variable(i,(z+j),constant);

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
                        	if(((z+j)->Iface+iface)->Direction==dir)
				{
					(((z+j)->Iface+iface)->ZSlab)->T = (double *)initializeSlabVarCoeff(1,(z+j),constant,4,dir,((z+j)->Iface+iface)->Point_Range);
                                        (((z+j)->Iface+iface)->ZSlab)->TOld = (double *)initializeSlabVarCoeff(1,(z+j),constant,4,dir,((z+j)->Iface+iface)->Point_Range);

                                        if((((z+j)->Iface+iface)->ZSlab)->T == NULL)
                                        {
                                        printf("Memory is not alloted to T, exiting\n");
                                        exit(1);
                                        }
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
                        	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{

					(((z+j)->ZoneBc+zBc)->ZSlab)->T = (double *)initializeSlabVarCoeff(1,(z+j),constant,4,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        (((z+j)->ZoneBc+zBc)->ZSlab)->TOld = (double *)initializeSlabVarCoeff(1,(z+j),constant,4,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}

	if(CgnsBase->flagDiffKinetCont==1) //Diffusion controlled reaction
	{

		for(j=0;j<num_zone;j++)
		{
                constant = 1.0;

                ((z+j)->Flowsol)->enthalpy = (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->FuelMassFrac = (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->OxiMassFrac = (double *)initialize_variable(i,(z+j),constant);
                for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                {
                        ((z+j)->Flowsol)->InMassFrac[k] = (double *)initialize_variable(i,(z+j),constant);
                }
                ((z+j)->Flowsol)->ProdMassFrac = (double *)initialize_variable(i,(z+j),constant);

		((z+j)->PrevFlowsol)->enthalpy= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);
                ((z+j)->PrevFlowsol)->OxiMassFrac = (double *)initialize_variable(i,(z+j),constant);
                for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                {
                	((z+j)->PrevFlowsol)->InMassFrac[k] = (double *)initialize_variable(i,(z+j),constant);
                }
                ((z+j)->PrevFlowsol)->ProdMassFrac = (double *)initialize_variable(i,(z+j),constant);

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;

                	for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
                        		if(((z+j)->Iface+iface)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab for Mixture Fraction(Z)*/
					(((z+j)->Iface+iface)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        (((z+j)->Iface+iface)->ZSlab)->FuelMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        (((z+j)->Iface+iface)->ZSlab)->OxiMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        (((z+j)->Iface+iface)->ZSlab)->OxiMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
						(((z+j)->Iface+iface)->ZSlab)->InMassFrac[k] = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
						(((z+j)->Iface+iface)->ZSlab)->InMassFracOld[k] = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
					(((z+j)->Iface+iface)->ZSlab)->ProdMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->ProdMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->enthalpy = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->enthalpyOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
                                	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Boundary Slab for Mixture Farction(Z)*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
						(((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac[k] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        }
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
					(((z+j)->ZoneBc+zBc)->ZSlab)->InMassFracOld[k] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        }
					(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpy = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpyOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
	}
        else if(CgnsBase->flagDiffKinetCont==2) //Kinetically Controlled
        {
        	if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                {
			for(j=0;j<num_zone;j++)
			{
                	constant = 1.0;

			((z+j)->Flowsol)->OxiMassFrac = (double *)initialize_variable(i,(z+j),constant);
                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                        {
                		((z+j)->Flowsol)->InMassFrac[k] = (double *)initialize_variable(i,(z+j),constant);
                        }
                	((z+j)->Flowsol)->ProdMassFrac = (double *)initialize_variable(i,(z+j),constant);

                	((z+j)->PrevFlowsol)->OxiMassFrac = (double *)initialize_variable(i,(z+j),constant);
                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                        {
                                ((z+j)->PrevFlowsol)->InMassFrac[k] = (double *)initialize_variable(i,(z+j),constant);
                        }
                	((z+j)->PrevFlowsol)->ProdMassFrac = (double *)initialize_variable(i,(z+j),constant);

			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;

                		for(dir=1;dir<=6;dir++)
				{
				for(iface=0;iface<numiface;iface++)
				{
                        		if(((z+j)->Iface+iface)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab for Mixture Fraction(Z)*/
                                        (((z+j)->Iface+iface)->ZSlab)->OxiMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        (((z+j)->Iface+iface)->ZSlab)->OxiMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                        		for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
					(((z+j)->Iface+iface)->ZSlab)->InMassFrac[k] = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
					(((z+j)->Iface+iface)->ZSlab)->InMassFracOld[k] = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
					(((z+j)->Iface+iface)->ZSlab)->ProdMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->ProdMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
                                	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Boundary Slab for Mixture Farction(Z)*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
                                                (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac[k] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        }
                                        for(k=0;k<CgnsBase->BaseParameters->NIngases;k++)
                                        {
					(((z+j)->ZoneBc+zBc)->ZSlab)->InMassFracOld[k] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                        }
					(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
				}
			}
                }
               /* else //Premixed
                {
			for(j=0;j<num_zone;j++)
			{
                	constant = 1.0;

			((z+j)->Flowsol)->OxiMassFrac = (double *)initialize_variable(i,(z+j),constant);
                	((z+j)->Flowsol)->ProdMassFrac = (double *)initialize_variable(i,(z+j),constant);

                	((z+j)->PrevFlowsol)->OxiMassFrac = (double *)initialize_variable(i,(z+j),constant);
                	((z+j)->PrevFlowsol)->ProdMassFrac = (double *)initialize_variable(i,(z+j),constant);

			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;

                		for(dir=1;dir<=6;dir++)
				{
				for(iface=0;iface<numiface;iface++)
				{
                        		if(((z+j)->Iface+iface)->Direction==dir)
					{
					//Memory allocation in the Interface Slab for Mixture Fraction(Z)
                                        (((z+j)->Iface+iface)->ZSlab)->OxiMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        (((z+j)->Iface+iface)->ZSlab)->OxiMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->ProdMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->ProdMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                        }
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
                                	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					//Memory allocation in the Boundary Slab for Mixture Farction(Z)
					(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
				}
			}
                }*/	
                for(j=0;j<num_zone;j++) //Making unused InMassFrac array to NULL
                {

                        for(k=CgnsBase->BaseParameters->NIngases;k<10;k++)
                        {

                                ((z+j)->Flowsol)->InMassFrac[k] = NULL;

                                ((z+j)->PrevFlowsol)->InMassFrac[k] = NULL;
                        }

                        numiface=(z+j)->Number_of_Interfaces;
                        numBc=(z+j)->Number_of_BC_Patches;

                                for(dir=1;dir<=6;dir++)
                                {
                                        for(iface=0;iface<numiface;iface++)
                                        {
                                                if(((z+j)->Iface+iface)->Direction==dir)
                                                {
                                                        for(k=CgnsBase->BaseParameters->NIngases;k<10;k++)
                                                        {
                                                                (((z+j)->Iface+iface)->ZSlab)->InMassFrac[k] = NULL;

                                                                (((z+j)->Iface+iface)->ZSlab)->InMassFracOld[k] = NULL;
                                                        }
                                                }
                                        }

                                        for(zBc=0;zBc<numBc;zBc++)
                                        {
                                                if(((z+j)->ZoneBc+zBc)->Direction==dir)
                                                {
                                                        for(k=CgnsBase->BaseParameters->NIngases;k<10;k++)
                                                        {
                                                                (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac[k] = NULL;

                                                                (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFracOld[k] = NULL;
                                                        }
                                                }
                                        }
                                }
                }

        }


return(0);

}



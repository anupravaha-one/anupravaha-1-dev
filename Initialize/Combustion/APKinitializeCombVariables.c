#include "APKinitializeCombVariables.h"
/**
This function is for initializing solvable combustion variables like mixture fraction/passive variable, fuel mass fraction and inert mass fraction.
*/

int initializeCombVariables(Base *CgnsBase, FILE *var2, int *count)
{

int i,j,dir,iface,numBc,numiface,zBc,num_zone;
double constant;
Zone *z;

 z=CgnsBase->CgnsZone;
 num_zone=CgnsBase->Number_of_zones;

 		if(CgnsBase->flagDiffKinetCont==1) //==1 for diffusion controlled reactions
		{

			/*The memory allocation for mixture fraction(Z) for Flowsol node and Slabs at boundary and Interface */
			printf("\nCREATING MEMORY FOR MIXTURE FRACTION(Z)");
			i=1;
			//fprintf(var2,"%d\n",i);

			CgnsBase->VARTAG[*count]=14;//Variable Tag-14 is for mixture fraction(Z) by shyam
			*count=*count+1;
			for(j=0;j<num_zone;j++)
			{
				printf("\nEnter the Initialization value of mixture fraction(Z):");
				scanf("%lf",&constant);
				//fprintf(var2,"%lf\n",constant);
				/*Memory allocation in the Flow Solution Node For Mixtur Fraction(Z)*/
				((z+j)->Flowsol)->MixFrac= (double *)initialize_variable(i,(z+j),constant);
				((z+j)->PrevFlowsol)->MixFrac= (double *)initialize_variable(i,(z+j),constant);

				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{	if(((z+j)->Iface+iface)->Direction==dir)
						{
						((z+j)->Iface+iface)->ConjugateRms = 0.0;
						/*Memory allocation in the Interface Slab for Mixture Fraction(Z)*/
						(((z+j)->Iface+iface)->ZSlab)->MixFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                (((z+j)->Iface+iface)->ZSlab)->MixFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Boundary Slab for Mixture Farction(Z)*/
						(((z+j)->ZoneBc+zBc)->ZSlab)->MixFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                (((z+j)->ZoneBc+zBc)->ZSlab)->MixFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						}
					}
				}
			}
			//Making fuel mass fraction pointer to NULL in all zones,because it is calculated from correlations for this type
			i=0;
			//fprintf(var2,"%d\n",i);
			for(j=0;j<num_zone;j++)
			{
				constant=0.0;
				//fprintf(var2,"%lf\n",constant);
				/*Memory allocation in the Flow Solution Node For Fuel mass Fraction(m_fu)*/
				((z+j)->Flowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);
				((z+j)->PrevFlowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);
                                //((z+j)->Flowsol)->InMassFrac= (double *)initialize_variable(i,(z+j),constant);
                                //((z+j)->PrevFlowsol)->InMassFrac= (double *)initialize_variable(i,(z+j),constant);


				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{	if(((z+j)->Iface+iface)->Direction==dir)
						{
						/*Memory allocation in the Interface Slab for Mixture Fraction(Z)*/
						//((z+j)->Iface+iface)->ConjugateRms = 0.0;
						(((z+j)->Iface+iface)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                //(((z+j)->Iface+iface)->ZSlab)->InMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Boundary Slab for Mixture Farction(Z)*/
						(((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                               // (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						}
					}
				}
			}


                }
                else if(CgnsBase->flagDiffKinetCont==2)
                {
                        if(CgnsBase->flagPreNonPre==2)
                        {
                                printf("\n CREATING MEMORY FOR MIXTURE FRACTION(Z)");
                                i=1;
                                //fprintf(var2,"%d\n",i);

                                CgnsBase->VARTAG[*count]=14;//Variable Tag-14 is for mixture fraction(Z) by shyam
                                *count=*count+1;
                                for(j=0;j<num_zone;j++)
                                {
                                        printf("\nEnter the Initialization value of mixture fraction(Z):");
                                        scanf("%lf",&constant);
                                        //fprintf(var2,"%lf\n",constant);
                                        /*Memory allocation in the Flow Solution Node For Mixtur Fraction(Z)*/
                                        ((z+j)->Flowsol)->MixFrac= (double *)initialize_variable(i,(z+j),constant);
                                        ((z+j)->PrevFlowsol)->MixFrac= (double *)initialize_variable(i,(z+j),constant);

                                        numiface=(z+j)->Number_of_Interfaces;
                                        numBc=(z+j)->Number_of_BC_Patches;
                                        for(dir=1;dir<=6;dir++)
                                        {
                                                for(iface=0;iface<numiface;iface++)
                                                {	if(((z+j)->Iface+iface)->Direction==dir)
                                                        {
                                                        //((z+j)->Iface+iface)->ConjugateRms = 0.0;
                                                        /*Memory allocation in the Interface Slab for Mixture Fraction(Z)*/
                                                        (((z+j)->Iface+iface)->ZSlab)->MixFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        (((z+j)->Iface+iface)->ZSlab)->MixFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        }
                                                }

                                                for(zBc=0;zBc<numBc;zBc++)
                                                {	if(((z+j)->ZoneBc+zBc)->Direction==dir)
                                                        {
                                                        /*Memory allocation in the Boundary Slab for Mixture Farction(Z)*/
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->MixFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->MixFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        }
                                                }
                                        }
                                }

                                /*The memory allocation for fuel mass fraction for Flowsol node and Slabs at boundary and Interface */
                                printf("\n CREATING MEMORY FOR FUEL MASS FRACTION");
                                i=1;
                                //fprintf(var2,"%d\n",i);
                                CgnsBase->VARTAG[*count]=15;//Variable Tag-15 is for fuel mass fraction by shyam
                                *count=*count+1;
                                for(j=0;j<num_zone;j++)
                                {
                                        printf("\nEnter the Initialization value of fuel mass fraction:");
                                        scanf("%lf",&constant);
                                        //constant = *(CgnsBase->BaseParameters->InitialMassFracs+0);
                                        //fprintf(var2,"%lf\n",constant);
                                        /*Memory allocation in the Flow Solution Node For Fuel mass Fraction(m_fu)*/
                                        ((z+j)->Flowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);
                                        ((z+j)->PrevFlowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);


                                        numiface=(z+j)->Number_of_Interfaces;
                                        numBc=(z+j)->Number_of_BC_Patches;
                                        for(dir=1;dir<=6;dir++)
                                        {
                                                for(iface=0;iface<numiface;iface++)
                                                {	if(((z+j)->Iface+iface)->Direction==dir)
                                                        {
                                                        /*Memory allocation in the Interface Slab for Mixture Fraction(Z)*/
                                                        ((z+j)->Iface+iface)->ConjugateRms = 0.0;
                                                        (((z+j)->Iface+iface)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        (((z+j)->Iface+iface)->ZSlab)->FuelMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        }
                                                }

                                                for(zBc=0;zBc<numBc;zBc++)
                                                {	if(((z+j)->ZoneBc+zBc)->Direction==dir)
                                                        {
                                                        /*Memory allocation in the Boundary Slab for Mixture Farction(Z)*/
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        }
                                                }
                                        }
                                }

                               /* i=0;
                                //fprintf(var2,"%d\n",i);
                                for(j=0;j<num_zone;j++)
                                {
                                        constant=0.0;
                                        //fprintf(var2,"%lf\n",constant);
                                        //Memory allocation in the Flow Solution Node For Fuel mass Fraction(m_fu)
                                        ((z+j)->Flowsol)->InMassFrac= (double *)initialize_variable(i,(z+j),constant);
                                        ((z+j)->PrevFlowsol)->InMassFrac= (double *)initialize_variable(i,(z+j),constant);


                                        numiface=(z+j)->Number_of_Interfaces;
                                        numBc=(z+j)->Number_of_BC_Patches;
                                        for(dir=1;dir<=6;dir++)
                                        {
                                                for(iface=0;iface<numiface;iface++)
                                                {	if(((z+j)->Iface+iface)->Direction==dir)
                                                        {
                                                        //Memory allocation in the Interface Slab for Mixture Fraction(Z)
                                                        //((z+j)->Iface+iface)->ConjugateRms = 0.0;
                                                        (((z+j)->Iface+iface)->ZSlab)->InMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        }
                                                }

                                                for(zBc=0;zBc<numBc;zBc++)
                                                {	if(((z+j)->ZoneBc+zBc)->Direction==dir)
                                                        {
                                                        //Memory allocation in the Boundary Slab for Mixture Farction(Z)
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        }
                                                }
                                        }
                                }*/


                        }

                   /*else
                   {
			printf("\n CREATING MEMORY FOR MIXTURE FRACTION(Z)");
			i=1;
			//fprintf(var2,"%d\n",i);

			CgnsBase->VARTAG[*count]=14;//Variable Tag-14 is for mixture fraction(Z) by shyam
			*count=*count+1;
			for(j=0;j<num_zone;j++)
			{
				printf("\nEnter the Initialization value of mixture fraction(Z):");
				scanf("%lf",&constant);
				//fprintf(var2,"%lf\n",constant);
				//Memory allocation in the Flow Solution Node For Mixtur Fraction(Z)
				((z+j)->Flowsol)->MixFrac= (double *)initialize_variable(i,(z+j),constant);
				((z+j)->PrevFlowsol)->MixFrac= (double *)initialize_variable(i,(z+j),constant);


				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{	if(((z+j)->Iface+iface)->Direction==dir)
						{
						((z+j)->Iface+iface)->ConjugateRms = 0.0;
						//Memory allocation in the Interface Slab for Mixture Fraction(Z)
						(((z+j)->Iface+iface)->ZSlab)->MixFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                (((z+j)->Iface+iface)->ZSlab)->MixFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						//Memory allocation in the Boundary Slab for Mixture Farction(Z)
						(((z+j)->ZoneBc+zBc)->ZSlab)->MixFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                (((z+j)->ZoneBc+zBc)->ZSlab)->MixFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						}
					}
				}
			}

			printf("\n CREATING MEMORY FOR FUEL MASS FRACTION");
			i=1;
			//fprintf(var2,"%d\n",i);
			CgnsBase->VARTAG[*count]=15;//Variable Tag-15 is for fuel mass fraction by shyam
			*count=*count+1;
			for(j=0;j<num_zone;j++)
			{
				//printf("\nEnter the Initialization value of fuel mass fraction:");
				//scanf("%lf",&constant);
                                constant = *(CgnsBase->BaseParameters->InitialMassFracs+0);
				//fprintf(var2,"%lf\n",constant);
				//Memory allocation in the Flow Solution Node For Fuel mass Fraction(m_fu)
				((z+j)->Flowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);
				((z+j)->PrevFlowsol)->FuelMassFrac= (double *)initialize_variable(i,(z+j),constant);

				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{	if(((z+j)->Iface+iface)->Direction==dir)
						{
						//Memory allocation in the Interface Slab for Mixture Fraction(Z)
                                                ((z+j)->Iface+iface)->ConjugateRms = 0.0;
						(((z+j)->Iface+iface)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                (((z+j)->Iface+iface)->ZSlab)->FuelMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						//Memory allocation in the Boundary Slab for Mixture Farction(Z)
						(((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                (((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						}
					}
				}
			}

                        printf("\n CREATING MEMORY FOR INERT GAS MASS FRACTION");
			i=1;
			//fprintf(var2,"%d\n",i);
			CgnsBase->VARTAG[*count]=16;//Variable Tag-16 is for inert gas mass fraction by shyam
			*count=*count+1;

                                for(j=0;j<num_zone;j++)
                                {
                                        //printf("\nEnter the Initialization value of fuel mass fraction:");
                                        //scanf("%lf",&constant);
                                        constant = *(CgnsBase->BaseParameters->InitialMassFracs+2);
                                        //fprintf(var2,"%lf\n",constant);
                                        //Memory allocation in the Flow Solution Node For Fuel mass Fraction(m_fu)
                                        ((z+j)->Flowsol)->InMassFrac= (double *)initialize_variable(i,(z+j),constant);
                                        ((z+j)->PrevFlowsol)->InMassFrac= (double *)initialize_variable(i,(z+j),constant);


                                        numiface=(z+j)->Number_of_Interfaces;
                                        numBc=(z+j)->Number_of_BC_Patches;
                                        for(dir=1;dir<=6;dir++)
                                        {
                                                for(iface=0;iface<numiface;iface++)
                                                {	if(((z+j)->Iface+iface)->Direction==dir)
                                                        {
                                                        //Memory allocation in the Interface Slab for Mixture Fraction(Z)
                                                        ((z+j)->Iface+iface)->ConjugateRms = 0.0;
                                                        (((z+j)->Iface+iface)->ZSlab)->InMassFrac = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        (((z+j)->Iface+iface)->ZSlab)->InMassFracOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
                                                        }
                                                }

                                                for(zBc=0;zBc<numBc;zBc++)
                                                {	if(((z+j)->ZoneBc+zBc)->Direction==dir)
                                                        {
                                                        //Memory allocation in the Boundary Slab for Mixture Farction(Z)
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        (((z+j)->ZoneBc+zBc)->ZSlab)->InMassFracOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
                                                        }
                                                }
                                        }
                                }


                   }*/
                }

return(0);

}

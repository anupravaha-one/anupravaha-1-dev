#include "APKinitialize_enthalpy.h"
int initialize_enthalpy(Base *CgnsBase,FILE *var2)
{
		int j,num_zone,dir,iface,numiface,numBc,zBc;
		double constant;
		Zone *z;
		z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
		num_zone=(CgnsBase)->Number_of_zones;
	
		for(j=0;j<num_zone;j++)
		{
		   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
			printf("\nEnter the Initialization value of Entahlpy:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		   } else { constant = 0.0;}
			/*Memory allocation in the Flow Solution Node For Temperature*/
			((z+j)->Flowsol)->enthalpy= (double *)initialize_variable(1,(z+j),constant);
			((z+j)->PrevFlowsol)->enthalpy= (double *)initialize_variable(1,(z+j),constant);
			((z+j)->Flowsol)->T= (double *)initialize_variable(0,(z+j),constant);
			((z+j)->PrevFlowsol)->T= (double *)initialize_variable(0,(z+j),constant);


			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab for Temperature*/
					(((z+j)->Iface+iface)->ZSlab)->enthalpy = (double *)initializeslab(1,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->enthalpyOld = (double *)initializeslab(1,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->T = (double *)initializeslab(0,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Boundary Slab for Temperature*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->T = (double *)initializeslab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->GradT = (double *)initializeslab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->TOld = (double *)initializeslab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpy = (double *)initializeslab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpyOld = (double *)initializeslab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
	return 0;
}



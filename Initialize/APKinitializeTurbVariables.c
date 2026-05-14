#include "APKinitializeTurbVariables.h"
/**
This function provides initial guess for the turbulence variables viz., turbulent kinetic energy (k) and dissipation rate of turbulent kinetic energy (epsilon) using the assumption of equillibrium between turbulent production and dissipation. Note that memory for the variables is not allocated here as it is already done in the function initialize().

This function is called in case the restart option is not choosen and only for the first time step. This function is called after calculating the turbulent production term without multiplying with the turbulent visocity. Here, the turbulent viscosity is replaced by the molecular viscosity to obtain the initial guess for epsilon.
*/
int initializeTurbVariables(Base *base)
{
	int nZones,numZones,i,j,k,NI,NJ,NK,*ptrSize,index;
	int flagRho,flagViscoeff;
	double mu,mu_t,*KE,*Epsilon,rho,Production,volume;
	Zone *zone;
	//Interface *IfacePtr;
	//ZoneBoundary *bcPtr;
	
	if(base->flagTurbModel == 1 || base->flagTurbModel == 3 || base->flagTurbModel == 4 || base->flagTurbModel == 5)
	{
		nZones = base->Number_of_zones;
		
		for(numZones=0;numZones<nZones;numZones++)
		{
			zone = base->CgnsZone+numZones;
			
			if(zone->flagSolidFluid==0)
			{
				ptrSize = zone->Size;
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
				
				KE = zone->Flowsol->Scal[0];
				Epsilon = zone->Flowsol->Scal[1];
				
				flagRho = zone->Flowsol->flagRho;
				flagViscoeff = zone->Flowsol->flagViscoeff;
				rho = *(zone->Flowsol->Rho);
				mu  = *(zone->Flowsol->Viscoeff);
				
				index = 0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					if(flagRho == 1)
						rho = *(zone->Flowsol->Rho+index);
					
					if(flagViscoeff == 1)
						mu  = *(zone->Flowsol->Viscoeff+index);
					
					mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+index);
					volume = *(zone->ZGeom->Volume+index);
					Production = *(zone->Flowsol->Tcoeff->TProduction+index);
					
					//mu_t = 10.0*mu;
					*(Epsilon+index) = mu_t * Production / volume;
					*(KE+index) = sqrt(mu_t * *(Epsilon+index) / (0.09*rho));
					
					index++;
				}
			}
		}
	}
	
	RefreshSurfaceCells(600,base);
	///Calculating friction velocity which is used in the boundary condition for epsilon
	if(!(base->flagTurbModel==4 || base->flagTurbModel==5))
	{
		CalcUtau(base);
	}
	RefreshSurfaceCells(601,base);
	
	return 0;
}

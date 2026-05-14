#include"APKTURBCalcGamma.h"
/**
This function calculates the value of (EpsilonOld/Kold) outside the dual time stepping loop of k-epsilon equations
and stores it in zone->Flowsol->Gamma for use in the dual time stepping loop.  This method is supposed to increase the 
stability as the value (EpsilonOld/Kold) is kept constant throught the dual time stepping loop.

NOTE: This method is currently deactivated due to the convergence problemes encountered as this method essentially
decouples the k-epsilon equations inside the dual time stepping loop which caused the solution (for both momentum and k-epsilon equations)
not to converge beyond the order of 1.0e-3 -> 1.0e-4.  This method also observed a wide range of fluctuations in both the
moementum and k-epsilon equations.
To activate this method remove the comments ('//') that have been put before CalcGamma() [LINE: 337] and freeGamma() [LINE: 390]
in solverSemiCoupled() and also use 'Gamma = *(zone->Flowsol->Gamma+INDEX);' [LINE: 113,257] instead of 'Gamma = EpsilonOld/Kold'
in the function kturb() [Do the similar thing in the function Omega() for k-Omega model];
However, on experience, it is not suggested to use this method.
*/
int CalcGamma(Base *base)
{
	int NI,NJ,NK,nz,cz;
	int i,j,k,*ptrSize,INDEX;
	Zone *zone;
	double Kold,Epsold,GammaMin,GammaMax,mu,rho,vol,TProduction;

	nz=(base->Number_of_zones);

	if((base->flagLogKEps == 0 || base->flagTurbModel == 2) && base->BaseParameters->Turb->numTurbScalars >= 2)	//Standard k-epsilon model
	{
		base->BaseParameters->Turb->TimeScaleMinLimitCount = 0;
		base->BaseParameters->Turb->TimeScaleMaxLimitCount = 0;
		for(cz=0;cz<nz;cz++)
		{
			zone = base->CgnsZone+cz;
			if(zone->flagSolidFluid==0)
			{
				ptrSize = zone->Size;
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);

				zone->Flowsol->Gamma = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				if(zone->Flowsol->Gamma == NULL)
				{
					printf("Memory Allocation Failed in CalcGamma.c\n");
					exit(100);
				}
				
				INDEX = 0;
				if(base->BaseParameters->Turb->flagLimTimeScale == 0)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						//Kold   = *(zone->PrevFlowsol->Scal[0]+INDEX);
						//Epsold = *(zone->PrevFlowsol->Scal[1]+INDEX);
						Kold   = *(zone->Flowsol->Scal[0]+INDEX);
						Epsold = *(zone->Flowsol->Scal[1]+INDEX);
						
						*(zone->Flowsol->Gamma+INDEX) = Epsold/Kold;
						INDEX++;
					}
				}
				else
				{
					rho = *(zone->PrevFlowsol->Rho);
					mu  = *(zone->Flowsol->Viscoeff);
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						Kold   = *(zone->Flowsol->Scal[0]+INDEX);
						Epsold = *(zone->Flowsol->Scal[1]+INDEX);

						*(zone->Flowsol->Gamma+INDEX) = Epsold/Kold;
						
						if(zone->Flowsol->flagRho == 1)
						{
							rho = *(zone->PrevFlowsol->Rho+INDEX);
						}
						if(zone->Flowsol->flagViscoeff == 1)
						{
							mu = *(zone->Flowsol->Viscoeff+INDEX);
						}
						vol = *(zone->ZGeom->Volume+INDEX);
						TProduction = *(zone->Flowsol->Tcoeff->TProduction+INDEX);
						
						GammaMin = 1.5*0.09*sqrt(8.0/3.0*(0.5*TProduction*rho/vol));
						GammaMax = sqrt((rho*Epsold)/mu);
						
						if(*(zone->Flowsol->Gamma+INDEX) < GammaMin)
						{
	// 						if(base->BaseParameters->Time/base->BaseParameters->TimeStep >=20)
							{
								*(zone->Flowsol->Gamma+INDEX) = GammaMin;
								base->BaseParameters->Turb->TimeScaleMaxLimitCount++;
	// 							printf("%s k=%d j=%d i=%d GammaMin=%e\n", zone->zonename,k,j,i,GammaMin);
							}
						}
						if(*(zone->Flowsol->Gamma+INDEX) > GammaMax)
						{
							*(zone->Flowsol->Gamma+INDEX) = GammaMax;
							base->BaseParameters->Turb->TimeScaleMinLimitCount++;
						}
						INDEX++;
					}
				}
			}
		}
	}

	return(100);
}

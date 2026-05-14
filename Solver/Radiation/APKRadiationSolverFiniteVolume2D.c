#include "APKRadiationSolverFiniteVolume2D.h"

/**After the preprocessing stage where reading of the grid CGNS file, asking for spherical discretization, and calculating the cell area, volume
and control angle, this function is called for 3D calculation of radiative problem in participating medium. This function again calls the functions
FillNewRadsolToPrevRadsol() for making the present calculation of Intensity to previous value for calculation of rms value of error between two iteration
of Gauss Seidel, CalculateIntensity(), where whole RTE is solved, FillOutgoingZoneBCIntensity(), to fill the boundary patch with outgoing ray Intensity, 
incoming ray intensity is calculated by boundary condition and finally ExtractVariables() to apply the boundary condition and handling the Interface patches.*/

int RadiationSolverFiniteVolume2D(Root CgnsRoot)
	{
	FILE *Conver;
	int i,j,k,l,m,nt,np,ni,nj,nk,ind_base, ind_zone, GSiterCount, nBase, nZone,cellindex;
	double ConCrit, zoneResidual, baseResidual, *IntNew, *IntOld, sumN, sumO;
	Zone *z;
	
	Conver=fopen("Settings/Convergence.dat","r");
	fscanf(Conver,"%lf",&ConCrit);
	fclose(Conver);
	
	nBase=CgnsRoot.Number_of_bases;
		
	GSiterCount=0;
	
	
	do
	{
		baseResidual=0.0;
		for(ind_base=0;ind_base<nBase;ind_base++)
		{
		
		nZone=(CgnsRoot.CgnsBase+ind_base)->Number_of_zones;
		
		nt=(CgnsRoot.CgnsBase+ind_base)->BaseParameters->ntheta;
		np=(CgnsRoot.CgnsBase+ind_base)->BaseParameters->nphi;
		
	       	zoneResidual=0.0;
		
			for(ind_zone=0;ind_zone<nZone;ind_zone++)
				{
				
				 z=(CgnsRoot.CgnsBase+ind_base)->CgnsZone+ind_zone;
				 
				 ni=*(z->Size+0);
				 nj=*(z->Size+1);
				 nk=*(z->Size+2);
					
				IntNew=z->Flowsol->RadIntensity;
				IntOld=z->PrevFlowsol->RadIntensity;
				  
				 
				 FillNewRadsolToPrevRadsol((CgnsRoot.CgnsBase+ind_base),z);
				
				 CalculateIntensity((CgnsRoot.CgnsBase+ind_base),z,16);
				 
				 FillOutgoingZoneBCIntensity((CgnsRoot.CgnsBase+ind_base), z);
				 
				 ExtractVariables(16, CgnsRoot.CgnsBase, z);
				 
			sumN=0.0;
				 sumO=0.0;
				 for(l=0;l<nt;l++)
				 for(m=0;m<np;m++)
				 	for(k=0;k<nk-1;k++)
					for(j=0;j<nj-1;j++)
					for(i=0;i<ni-1;i++)
						{
						cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						sumN=sumN+*(IntNew+cellindex*nt*np+np*l+m);
						sumO=sumO+*(IntOld+cellindex*nt*np+np*l+m);
						}
				 
			zoneResidual=zoneResidual+fabs(sumN-sumO)/fabs(sumN); 
				}
		baseResidual=baseResidual+zoneResidual;
		}
		
		GSiterCount++;
		if(GSiterCount%10==0)
		printf("\tIntensityResidual=%e\n",baseResidual);
	}while(baseResidual>ConCrit);
	
	//printDirections(CgnsRoot.CgnsBase);
	//printRadiationSlabs(CgnsRoot.CgnsBase);
		 
	return 0;
	}

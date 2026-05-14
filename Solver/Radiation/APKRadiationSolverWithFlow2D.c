#include "APKRadiationSolverWithFlow2D.h"

/** This function is only called when RTE is solved with the flow and does the same thing as the function RadiationSolverFiniteVolume2D()*/

int RadiationSolverWithFlow2D(Base *CgnsBase)
	{
	FILE *Conver;
	int i,j,k,l,m,nt,np,ni,nj,nk,ind_zone, GSiterCount, nZone,cellindex;
	double ConCrit, baseResidual, *IntNew, *IntOld, sumN, sumO;
	Zone *z;
	
	Conver=fopen("Settings/Convergence.dat","r");
	fscanf(Conver,"%lf",&ConCrit);
	fclose(Conver);
	
	
	GSiterCount=0;
	printf("\n\tSolving For Radiative Heat Transfer\n");
	
	do
	{
		baseResidual=0.0;
		
		
		nZone=CgnsBase->Number_of_zones;
		
		nt=CgnsBase->BaseParameters->ntheta;
		np=CgnsBase->BaseParameters->nphi;
		
	        	for(ind_zone=0;ind_zone<nZone;ind_zone++)
				{
				
				 z=CgnsBase->CgnsZone+ind_zone;
				 
				 ni=*(z->Size+0);
				 nj=*(z->Size+1);
				 nk=*(z->Size+2);
					
				IntNew=z->Flowsol->RadIntensity;
				IntOld=z->PrevFlowsol->RadIntensity;
				  
				 
				 FillNewRadsolToPrevRadsol(CgnsBase,z);
				
				 CalculateIntensity(CgnsBase,z,16);
				 
				 FillOutgoingZoneBCIntensity(CgnsBase, z);
				 
				 ExtractVariables(16, CgnsBase, z);
				 
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
				 
				baseResidual=baseResidual+fabs(sumN-sumO)/fabs(sumN);
			}
		
		GSiterCount++;
		if(GSiterCount%50==0)
		printf("\tIntensityResidual=%e\n",baseResidual);
	}while(baseResidual>ConCrit);
	
	//printDirections(CgnsRoot.CgnsBase);
	//printRadiationSlabs(CgnsBase);
	 
	return 0;
	}

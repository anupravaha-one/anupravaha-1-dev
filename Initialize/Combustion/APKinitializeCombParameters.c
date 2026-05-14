#include "APKinitializeCombParameters.h"

/**
This function, asks for no. of species that participate in the reaction, Mol. wts. of species, no. of moles of species that participate in the stoichometric single step reaction, Cp-coefficients of individual species and initial concentrations of species. Calculates stoichometric constant, net initial Cp.
*/

int initializeCombParameters(Base *CgnsBase, FILE *var2)
{
  int nzones,zonecount,i,j,NI,NJ,NK,*size,totalcells,numiface,numBc,dir,iface,zBc;
  int *Point_Range;
  double temp,*Cp,netCp,*MolWteq,*Cpi,A,B;
  //double *Cp1;
  Zone *z;
	  double a,b,c,d,e,f,g,h,T;			//nitesh

	printf("Enter the initial temperature(K)\n");
        scanf("%lf",&temp);
	printf("Enter no. of species participate in combustion(including inert gases)\n");
        scanf("%d",&(CgnsBase->BaseParameters->Nspecies));
        printf("Enter no. of inert gases present(Both in fuel stream + air stream)\n");
        scanf("%d",&(CgnsBase->BaseParameters->NIngases));
        fprintf(var2,"%d\n",CgnsBase->BaseParameters->NIngases);
        printf("Enter whether the fuel is a single one or mixture of fuels(0->single; 1->mixture of fuels)\n");
        scanf("%d",&(CgnsBase->BaseParameters->flagfueltype));

        CgnsBase->BaseParameters->MolWt = (double *)malloc(CgnsBase->BaseParameters->Nspecies*sizeof(double));
	CgnsBase->BaseParameters->Cp_ref = (double *)malloc(CgnsBase->BaseParameters->Nspecies*sizeof(double));  //nitesh
        CgnsBase->BaseParameters->Moles = (double *)malloc(CgnsBase->BaseParameters->Nspecies*sizeof(double));
	//CgnsBase->BaseParameters->CpCoeff = (double *)malloc(CgnsBase->BaseParameters->Nspecies*3*sizeof(double));
	CgnsBase->BaseParameters->CpCoeff = (double *)malloc(CgnsBase->BaseParameters->Nspecies*8*sizeof(double));    //nitesh

        CgnsBase->BaseParameters->rhoinfinity = (double *)malloc(sizeof(double));
        for(i=0;i<CgnsBase->BaseParameters->NIngases;i++)
        {
                CgnsBase->BaseParameters->InertgasFuelInlet[i] = (double *)malloc(sizeof(double));
                CgnsBase->BaseParameters->InertgasOxiInlet[i] = (double *)malloc(sizeof(double));
        }
        if(CgnsBase->flagPreNonPre == 2)
        {
	CgnsBase->BaseParameters->FuelMassFracInlet = (double *)malloc(sizeof(double));
        CgnsBase->BaseParameters->OxiMassFracInlet = (double *)malloc(sizeof(double));
        CgnsBase->BaseParameters->Zst = (double *)malloc(sizeof(double));
        }
        else
        {
	CgnsBase->BaseParameters->FuelMassFracInlet = NULL;
        CgnsBase->BaseParameters->OxiMassFracInlet  = NULL;
        CgnsBase->BaseParameters->Zst = NULL;
        }
        /*if(CgnsBase->flagDiffKinetCont==1)
        {
		CgnsBase->BaseParameters->enthalpyFuelInlet = (double *)malloc(sizeof(double));
                CgnsBase->BaseParameters->enthalpyOxiInlet = (double *)malloc(sizeof(double));
        }
        else
        {
		CgnsBase->BaseParameters->enthalpyFuelInlet = NULL;
                CgnsBase->BaseParameters->enthalpyOxiInlet = NULL;
        }*/

	CgnsBase->BaseParameters->enthalpyFuelInlet = (double *)malloc(sizeof(double));                            //nitesh
        CgnsBase->BaseParameters->enthalpyOxiInlet = (double *)malloc(sizeof(double));                             //nitesh
	
      
        CgnsBase->BaseParameters->delH = (double *)malloc(sizeof(double));

	z = CgnsBase->CgnsZone;

        //Calculation of ambient density
	*(CgnsBase->BaseParameters->rhoinfinity) = (101325.0*37.9968064)/(8314.0*temp);
	

        printf("Enter Heat of the Reaction\n");
        scanf("%lf",CgnsBase->BaseParameters->delH);
        fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->delH));

        printf("Enter parameters(Mol.Wt,Moles,CoCoeff) of fuel,oxidant,inertgas, and products.\n");
        printf("Maintain the above order while entering the parameters related to species.\n");

        for(j=0;j<CgnsBase->BaseParameters->Nspecies;j++)
        {
        	printf("Enter Mol.Wt of species[%d]\n",j+1);
                scanf("%lf",CgnsBase->BaseParameters->MolWt+j);
                fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->MolWt+j));

		/*printf("Enter Specific heat Cp of species [%d] at reference(initial) temperature\n",j+1);
                scanf("%lf",CgnsBase->BaseParameters->Cp_ref+j);
                fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->Cp_ref+j));                                          //nitesh
		*/
		
		printf("Enter No. of Moles of species[%d]\n",j+1);
               	scanf("%lf",CgnsBase->BaseParameters->Moles+j);
                fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->Moles+j));


//                 printf("Enter three CpCoeffs. of species[%d]\n",j+1);
//         	for(i=0;i<3;i++)
//                 {
//                         scanf("%lf",CgnsBase->BaseParameters->CpCoeff+3*j+i);
//                         fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->CpCoeff+3*j+i));
//                 }
		
		printf("Enter eight CpCoeffs. of species[%d]\n",j+1);                                                   //nitesh
        	for(i=0;i<8;i++)
                {
                        scanf("%lf",CgnsBase->BaseParameters->CpCoeff+8*j+i);
                        fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->CpCoeff+8*j+i));
                }
        }
        if(CgnsBase->BaseParameters->FuelMassFracInlet != NULL)
        {
	printf("Enter  mass fraction of fuel at fuel stream\n");
        scanf("%lf",CgnsBase->BaseParameters->FuelMassFracInlet);
	fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->FuelMassFracInlet));
        }
        if(CgnsBase->BaseParameters->OxiMassFracInlet!= NULL)
        {
        printf("Enter  mass fraction of oxidant at oxidant stream\n");
        scanf("%lf",CgnsBase->BaseParameters->OxiMassFracInlet);
        fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->OxiMassFracInlet));
        }
	for(i=0;i<CgnsBase->BaseParameters->NIngases;i++)
        {
        	printf("Enter mass fraction of inert gas[%d] at fuel stream\n",i+1);
                scanf("%lf",CgnsBase->BaseParameters->InertgasFuelInlet[i]);
                printf("Enter mass fraction of inert gas[%d] at oxygen stream\n",i+1);
                scanf("%lf",CgnsBase->BaseParameters->InertgasOxiInlet[i]);
        }
        if(CgnsBase->BaseParameters->enthalpyFuelInlet != NULL)
        {
	printf("Enter  enthalpy at fuel stream(hf)\n");
        scanf("%lf",CgnsBase->BaseParameters->enthalpyFuelInlet);
	fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->enthalpyFuelInlet));
        }
        if(CgnsBase->BaseParameters->enthalpyOxiInlet  != NULL)
        {
        printf("Enter  enthalpy at oxidant stream(ho)\n");
        scanf("%lf",CgnsBase->BaseParameters->enthalpyOxiInlet);
        fprintf(var2,"%lf\n",*(CgnsBase->BaseParameters->enthalpyOxiInlet));
        }


        nzones = CgnsBase->Number_of_zones;


        //CgnsBase->BaseParameters->InitialMassFracs = (double *)malloc(CgnsBase->BaseParameters->Nspecies*sizeof(double));
        CgnsBase->BaseParameters->InitialTemp = (double *)malloc(sizeof(double));
	*(CgnsBase->BaseParameters->InitialTemp) = temp;
        /*printf("Enter initial mass fractions of fuel, oxygen, inertgases, carbondioxide and water\n");
        for(i=0;i<CgnsBase->BaseParameters->Nspecies;i++)
        {
        	printf("Enter initial mass fraction of species[%d]\n",i+1);
                scanf("%lf",CgnsBase->BaseParameters->InitialMassFracs+i);
        }*/

        CgnsBase->BaseParameters->Stoich = (double *)malloc(sizeof(double));
	printf("Enter Stoichiometric Constant(O/F)\n");
        scanf("%lf",(CgnsBase->BaseParameters->Stoich));
      	//*(CgnsBase->BaseParameters->Stoich) = (*(CgnsBase->BaseParameters->Moles+1)*(*(CgnsBase->BaseParameters->MolWt+1)))/((*(CgnsBase->BaseParameters->Moles+0))*(*(CgnsBase->BaseParameters->MolWt+0)));
        //printf("\nStoich. constant = %lf\n",*(CgnsBase->BaseParameters->Stoich));
        if(CgnsBase->BaseParameters->Zst != NULL)
        {
        A = ((*(CgnsBase->BaseParameters->FuelMassFracInlet))/(*(CgnsBase->BaseParameters->OxiMassFracInlet)));
        B = 1.0 + (*(CgnsBase->BaseParameters->Stoich))*A;
         *(CgnsBase->BaseParameters->Zst) = 1.0/B;
        // printf("\nZst = %lf\n",*(CgnsBase->BaseParameters->Zst));
        }



	//Calculation of initial individual cps of species from initial temperature

                for(zonecount=0;zonecount<nzones;zonecount++)
       	 	{
			size = (CgnsBase->CgnsZone+zonecount)->Size;
                	NI = *(size+0);
                	NJ = *(size+1);
               		NK = *(size+2);

                	totalcells = (NI-1)*(NJ-1)*(NK-1);
                        ((z+zonecount)->Flowsol)->Tcoeff=(TurbCoeff *)malloc(sizeof(TurbCoeff));
                        (((z+zonecount)->Flowsol)->Tcoeff)->Cps = (double *)malloc((CgnsBase->BaseParameters->Nspecies)*totalcells*sizeof(double));
                        Cp = (((z+zonecount)->Flowsol)->Tcoeff)->Cps;
			for(j=0;j<CgnsBase->BaseParameters->Nspecies;j++)
                        {
				a=*(CgnsBase->BaseParameters->CpCoeff+j*8+0);                                            //nitesh
				b=*(CgnsBase->BaseParameters->CpCoeff+j*8+1);
				c=*(CgnsBase->BaseParameters->CpCoeff+j*8+2);
				d=*(CgnsBase->BaseParameters->CpCoeff+j*8+3);
				e=*(CgnsBase->BaseParameters->CpCoeff+j*8+4);
				f=*(CgnsBase->BaseParameters->CpCoeff+j*8+5);
				g=*(CgnsBase->BaseParameters->CpCoeff+j*8+6);
				h=*(CgnsBase->BaseParameters->CpCoeff+j*8+7);
				T=*(CgnsBase->BaseParameters->InitialTemp);

				*(CgnsBase->BaseParameters->Cp_ref+j) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0 + h*T*T*T*T*T*T*T/8.0;          //nitesh
				//*(CgnsBase->BaseParameters->Cp_ref+j) = a + b*T/1.0 + c*T*T/1.0 + d*T*T*T/1.0 + e*T*T*T*T/1.0 + f*T*T*T*T*T/1.0 + g*T*T*T*T*T*T/1.0 + h*T*T*T*T*T*T*T/1.0;          //nitesh
				//*(CgnsBase->BaseParameters->Cp_ref+2) = 1300.0;
				//printf("Cp_ref [%d] = %lf",j,*(CgnsBase->BaseParameters->Cp_ref+j));
				



				for(i=0;i<totalcells;i++)
                        	{
					//*(Cp+j*totalcells+i) = *(CgnsBase->BaseParameters->CpCoeff+j*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+j*3+1))*temp + (*(CgnsBase->BaseParameters->CpCoeff+j*3+2))*temp*temp;
					*(Cp+j*totalcells+i) = *(CgnsBase->BaseParameters->Cp_ref+j);
                        	}
                        }
			//printf("Iam at the end of cal. of cps in the domain\n");
                        numiface=(z+zonecount)->Number_of_Interfaces;

		        numBc=(z+zonecount)->Number_of_BC_Patches;
                        for(dir=1;dir<=6;dir++)
                        {
                                for(iface=0;iface<numiface;iface++)
                                {
                                        if(((z+zonecount)->Iface+iface)->Direction==dir)
                                        {
                                                totalcells = 0;
                                                Point_Range = ((z+zonecount)->Iface+iface)->Point_Range;
                                                if(dir==1||dir==2)
                                                totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
                                                if(dir==3||dir==4)
                                                totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
                                                if(dir==5||dir==6)
                                                totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

                                                ((((z+zonecount)->Iface+iface)->ZSlab)->Cps)=(double *)malloc((CgnsBase->BaseParameters->Nspecies)*totalcells*sizeof(double));

                                                Cpi = ((((z+zonecount)->Iface+iface)->ZSlab)->Cps);
                                                for(j=0;j<CgnsBase->BaseParameters->Nspecies;j++)
                                                {
							a=*(CgnsBase->BaseParameters->CpCoeff+j*8+0);                                  //nitesh
							b=*(CgnsBase->BaseParameters->CpCoeff+j*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+j*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+j*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+j*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+j*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+j*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+j*8+7);
							T=*(CgnsBase->BaseParameters->InitialTemp);
							//T=423.0;
							*(CgnsBase->BaseParameters->Cp_ref+j) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0 + h*T*T*T*T*T*T*T/8.0;          //nitesh
							
                                                        for(i=0;i<totalcells;i++)
                                                        {
                                                                //*(Cpi+j*totalcells+i) = *(CgnsBase->BaseParameters->CpCoeff+j*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+j*3+1))*temp + (*(CgnsBase->BaseParameters->CpCoeff+j*3+2))*temp*temp;
								*(Cpi+j*totalcells+i) = *(CgnsBase->BaseParameters->Cp_ref+j);           //nitesh
								
                                                        }
                                                }

                                        }
                                }

                                for(zBc=0;zBc<numBc;zBc++)
                                {
                                        if(((z+zonecount)->ZoneBc+zBc)->Direction==dir)
                                        {
                                                totalcells = 0;
                                                Point_Range = ((z+zonecount)->ZoneBc+zBc)->Point_Range;
                                                if(dir==1||dir==2)
                                                totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
                                                if(dir==3||dir==4)
                                                totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
                                                if(dir==5||dir==6)
                                                totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

                                                (((z+zonecount)->ZoneBc+zBc)->ZSlab)->Cps=(double *)malloc((CgnsBase->BaseParameters->Nspecies)*totalcells*sizeof(double));

                                                Cpi = (((z+zonecount)->ZoneBc+zBc)->ZSlab)->Cps;
                                                for(j=0;j<CgnsBase->BaseParameters->Nspecies;j++)
                                                {
							a=*(CgnsBase->BaseParameters->CpCoeff+j*8+0);                                 //nitesh
							b=*(CgnsBase->BaseParameters->CpCoeff+j*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+j*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+j*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+j*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+j*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+j*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+j*8+7);
							T=*(CgnsBase->BaseParameters->InitialTemp);
							//T=423.0;
							*(CgnsBase->BaseParameters->Cp_ref+j) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0 + h*T*T*T*T*T*T*T/8.0;          //nitesh
							
                                                        for(i=0;i<totalcells;i++)
                                                        {
                                                                //*(Cpi+j*totalcells+i) = *(CgnsBase->BaseParameters->CpCoeff+j*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+j*3+1))*temp + (*(CgnsBase->BaseParameters->CpCoeff+j*3+2))*temp*temp;
								*(Cpi+j*totalcells+i) = *(CgnsBase->BaseParameters->Cp_ref+j);            //nitesh
                                                       
                                                        }
                                                }

                                        }
                                }
                        }

        	}


       
      		  netCp = 1.0;    
		//Initialization of specific heat(net Cp) and MolWteq
		for(zonecount=0;zonecount<nzones;zonecount++)
       	 	{
			size = (CgnsBase->CgnsZone+zonecount)->Size;
                	NI = *(size+0);
                	NJ = *(size+1);
               		NK = *(size+2);

                	totalcells = (NI-1)*(NJ-1)*(NK-1);

                        ((CgnsBase->CgnsZone+zonecount)->Flowsol)->SpecificHeat = (double *)malloc(totalcells*sizeof(double));
                        ((CgnsBase->CgnsZone+zonecount)->PrevFlowsol)->SpecificHeat = (double *)malloc(totalcells*sizeof(double));
                        (((CgnsBase->CgnsZone+zonecount)->Flowsol)->Tcoeff)->MolWteq = (double *)malloc(totalcells*sizeof(double));
                        Cp =  ((CgnsBase->CgnsZone+zonecount)->Flowsol)->SpecificHeat;
                        MolWteq = (((CgnsBase->CgnsZone+zonecount)->Flowsol)->Tcoeff)->MolWteq;
			((CgnsBase->CgnsZone+zonecount)->Flowsol)->flagSpecificHeat = 1;

                        for(i=0;i<totalcells;i++)
                        {
				*(Cp+i) = netCp;
				 *(MolWteq+i) = *(CgnsBase->BaseParameters->MolWt+1); //initially oxidant Mol. Wt is assigned throughout the domain
                        }

			numiface=(z+zonecount)->Number_of_Interfaces;

		        numBc=(z+zonecount)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{
							if(((z+zonecount)->Iface+iface)->Direction==dir)
							{
								(((z+zonecount)->Iface+iface)->ZSlab)->SpecificHeat=(double *)initializeSlabVarCoeff(1,(z+zonecount),netCp,3,dir,((z+zonecount)->Iface+iface)->Point_Range);
								//printf("I am leaving Iface\n");
                                                                (((z+zonecount)->Iface+iface)->ZSlab)->MolWteq = (double *)initializeSlabVarCoeff(1,(z+zonecount),*(CgnsBase->BaseParameters->Moles+1),3,dir,((z+zonecount)->Iface+iface)->Point_Range);
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+zonecount)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+zonecount)->ZoneBc+zBc)->ZSlab)->SpecificHeat = (double *)initializeSlabVarCoeff(1,(z+zonecount),netCp,3,dir,((z+zonecount)->ZoneBc+zBc)->Point_Range);
							(((z+zonecount)->ZoneBc+zBc)->ZSlab)->MolWteq = (double *)initializeSlabVarCoeff(1,(z+zonecount),*(CgnsBase->BaseParameters->Moles+1),3,dir,((z+zonecount)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}


        	}

		for(i =CgnsBase->BaseParameters->NIngases;i<10;i++ ) //Making unused InertgasFuelInlet array and InertgasOxiInlet to NULL
                {
			CgnsBase->BaseParameters->InertgasFuelInlet[i] = NULL;
                        CgnsBase->BaseParameters->InertgasOxiInlet[i] = NULL;
                }


return(0);



}

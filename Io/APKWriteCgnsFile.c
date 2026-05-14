#include "APKWriteCgnsFile.h"

/**
This function WriteCgnsFile() generates the 3D output solution file in CGNS format.
 
This function uses 
the Mid-Level library functions (Document version 2.3.6 CGNS Version 2.2 ) the output file is compatible to 
CGNS version 2.2. The inputs to this function are pointer to the name of the file to be  written and the Root 
variable of the Data structure.
(1)	Checks for the flag for 2 dimensional flow 
	if flow is 2D skips further operations and calles the 2D function
	else
(2)     Opens the CGNS file in "WRITE" mode of the name for which  pointer has been supplied to this function.
(3)	Loops over all the zones and pulls out the appropriate data arrays using the Root variable 
	and currect links to the data locations, and writes the zonal data in the file.
(4)     closes the Loop and  CGNS file.
*/

int WriteCgnsFile(char* CgnsFileName,Root CgnsRoot)
{
 
	/* --------------- Variable Declaration --------------- */
	BCType_t type=8;
	BCDataType_t type1=2,type2=3;
	DataType_t soldtype=4;
	DataType_t ndt=2;
	int index_file,index_base,index_zone,index_grid,index_coord,index_bc,index_iface,indexbc,index_dset;
	/* Dummy Variables */
	int ier,nbc,base_ind,zone_ind,grid_ind,interface_ind,no,nlist,count;
	double *solution;
	int npnts=2,norind[]={-1,0,0};
	int *sizeptr,size_patch[1],dir;
	int datadim=1;
	int fluxsize[1];
	int size1,size2;
	PointSetType_t ptsettype=4;
	GridLocation_t location=3,locab=4;
	int sol_index,solarr_index,scal_ind,scalCoeff_ind;
	char floname[]="Solutions";
	char s1[]="Scalar1";
	char s2[]="Scalar2";
	char s3[]="Scalar3";
	char s4[]="Scalar4";
	char s5[]="Scalar5";
	char s6[]="Scalar6";
	char s7[]="Scalar7";
	char s8[]="Scalar8";
	char s9[]="Scalar9";
	char s10[]="Scalar10";
	char scoeff1[]="ScalarCoeff1";
	char scoeff2[]="ScalarCoeff2";
	char scoeff3[]="ScalarCoeff3";
	char scoeff4[]="ScalarCoeff4";
	char scoeff5[]="ScalarCoeff5";
	char scoeff6[]="ScalarCoeff6";
	char scoeff7[]="ScalarCoeff7";
	char scoeff8[]="ScalarCoeff8";
	char scoeff9[]="ScalarCoeff9";
	char scoeff10[]="ScalarCoeff10";
	char dau[]="VelocityX";
	char dav[]="VelocityY";
	char daw[]="VelocityZ";
	// char daute[]="VelocityXt";
       // char davte[]="VelocityYt";
        //char dawte[]="VelocityZt";


	char rho[]="Density";
	char energysource[]="EnergySource";
	char visc[]="ViscosityMolecular";
	char cond[]="ThermalConductivity";
	char spheat[]="SpecificHeatPressure";
	char dat[]="Temperature";
	char dap[]="Pressure";
	char daum[]="MasVelocityX";
	char davm[]="MasVelocityY";
	char dawm[]="MasVelocityZ";
	char fluxt[]="TotalHeatFluxOnBoundaryPatch";
	char Fluxs1[]="TotalFluxS1onBoundaryPatch";
	char Fluxs2[]="TotalFluxS2onBoundaryPatch";
	char Fluxs3[]="TotalFluxS3onBoundaryPatch";
	char Fluxs4[]="TotalFluxS4onBoundaryPatch";
	char Fluxs5[]="TotalFluxS5onBoundaryPatch";
	char Fluxs6[]="TotalFluxS6onBoundaryPatch";
	char Fluxs7[]="TotalFluxS7onBoundaryPatch";
	char Fluxs8[]="TotalFluxS8onBoundaryPatch";
	char Fluxs9[]="TotalFluxS9onBoundaryPatch";
	char Fluxs10[]="TotalFluxS10onBoundaryPatch";
	char gradt[]="Temperaturegradiant";
	char grads1[]="Scalar1gradiant";
	char grads2[]="Scalar2gradiant";
	char grads3[]="Scalar3gradiant";
	char grads4[]="Scalar4gradiant";
	char grads5[]="Scalar5gradiant";
	char grads6[]="Scalar6gradiant";
	char grads7[]="Scalar7gradiant";
	char grads8[]="Scalar8gradiant";
	char grads9[]="Scalar9gradiant";
	char grads10[]="Scalar10gradiant";
	char viscTurb[]="ViscosityTurbulent";
	char effviscTurb[]="EffectiveViscosity";
	char turbProdterm[]="TurbulentProductionTerm";
	char PhiK[]="PhiK";
	char PhiEps[]="PhiEps";
	char vortX[]="VorticityX";
	char vortY[]="VorticityY";
	char vortZ[]="VorticityZ";
	char dset[]="BCDSet";
	char MixtureFraction[] = "MixFrac";
	char FuelMassFraction[] = "FuelMassFrac";
	char OxiMassFrac[] = "OxiMassFrac";
	//    char InMassFrac[] = "InMassFrac";
	char ProdMassFrac[] = "ProdMassFrac";
	//name of boundary data set structures
	
	//Electrical Parameters
	char daphi[]="ElectricPotential";
	char daEfx[]="ElectricFieldX";
	char daEfy[]="ElectricFieldY";
	char daEfz[]="ElectricFieldZ";
	char daCharge[]="ChargeDensity";
	
	//LMMHD Parameters naveen
	char daphiMHD[]="ElectricPotentialMHD";
	char daMfx[]="MagneticFieldX";
	char daMfy[]="MagneticFieldY";
	char daMfz[]="MagneticFieldZ";
	char daCurrentx[]="CurrentDensityX";
	char daCurrenty[]="CurrentDensityY";
	char daCurrentz[]="CurrentDensityZ";
	char daTotMfx[]="TotalMagneticFieldX";//rajendra
	char daTotMfy[]="TotalMagneticFieldY";
	char daTotMfz[]="TotalMagneticFieldZ";
	char dapsi[]="PsiMHD";	
	//groundwater variables pawan  
	char dah[]="PressureHead";
	char dawa[]="WaterContent";
	char rel[]="UnSatPermeability";
	Base *base;
	Zone *zone;
	
	//sushil loop
	int nophase,nop;
	char *temp;
	char *tempx,*tempy,*tempz;
	char *P1,*P2,*P3,*tempx1;
	char tempvisc[]="Phaseviscosity";
	char tempden[] ="Phasedensity";
	char tempV[]="Vartheta";
	
 
 
 //For overall momentum
	char ovmx[]="momentum-X";
	char ovmy[]="momentum-Y";
	char ovmz[]="momentum-Z";
 
 
	char ovmxs[]="momentumstar-X";
	char ovmys[]="momentumstar-Y";
	char ovmzs[]="momentumstar-Z";
	
	temp=(char *)malloc(40*sizeof(char));
	tempx=(char *)malloc(40*sizeof(char));
	tempy=(char *)malloc(40*sizeof(char));
	tempz=(char *)malloc(40*sizeof(char));
	
	
	
	
	nop=CgnsRoot.CgnsBase->nphase;

	P1 = tempx;
	P2 =tempy;
	P3 =tempz;		

	
	//end 
 //  ####################################

	/* Opening the CGNS file in WRITE mode */	
	fluxsize[0]=1;
	ier = cg_open(CgnsFileName,MODE_WRITE,&index_file);
	CGNSErrorCheck(ier);	/* Checking for CGNS Error */
	
	/* fprintf(sum,"1) CGNS File Name : %s\n",CgnsFileName);*/

	no=(CgnsRoot.CgnsBase)->Number_of_zones;

	/* Loop Over Base Nodes for reading the nodes under each Base Node */
	for(index_base=1;index_base<=CgnsRoot.Number_of_bases;index_base++)
	{
		base = CgnsRoot.CgnsBase+(index_base-1);
		if(base->BaseParameters->Turb->numTurbScalars>=1)
		{
			if(base->flagTurbModel==2 || base->flagTurbModel==6)
			{
				strcpy(s1,"TKE");
				strcpy(s2,"Omega");
			}
			else if(base->flagTurbModel==8)
			{
				strcpy(s1,"ModEddyVisc");
			}
			else
			{
				strcpy(s1,"TKE");
				strcpy(s2,"Epsilon");
			}
		}
		
		if (((CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension) < ((CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension)) {
		   printf("\nWARNING: CGNS Mesh physical dimension is less than cell dimension. \nChanging Mesh physical dimension equal to cell dimesion.... \n");
		   (CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension = (CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension ;
		}
		ier = cg_base_write(index_file,(CgnsRoot.CgnsBase+(index_base-1))->basename,(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension,(CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension,&base_ind);
		CGNSErrorCheck(ier);
		/* Loop Over Zone Nodes for reading the nodes under each Zone Node */ 
		for(index_zone=1;index_zone<=(CgnsRoot.CgnsBase+(index_base-1))->Number_of_zones;index_zone++)
		{
			zone = base->CgnsZone+(index_zone-1);
			/*writing the Zone Name & Zone Size */
			ier =cg_zone_write(index_file,index_base,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->zonename,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Zonetype,&zone_ind);	
			CGNSErrorCheck(ier);


			/*##########################################GRID CO-ORDINATES WRITTING####################################*/
			/* ... Loop over all the Grid Coordinate Nodes for reading Grid Coordinates */
			for(index_grid=1;index_grid<=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes;index_grid++)
			{
				
				ier=cg_grid_write(index_file,index_base,index_zone,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->GridcoordinateName,&grid_ind);
				CGNSErrorCheck(ier);
					
				/* .......Writing the 1st Coordinate */      
				ier = cg_coord_write(index_file,index_base,index_zone,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype,  
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name1, 
						((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X),&index_coord );
				CGNSErrorCheck(ier);
				
				if(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates) >=2)
				{
					/* ....... 2nd Coordinate */
				
				
					/* .......Writting the 2nd Coordinate */      
					ier = cg_coord_write(index_file,index_base,index_zone,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype,  
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name2, 
						((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y),&index_coord);
					CGNSErrorCheck(ier); 
				}
			
				if(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates) >=3)
				{
					/*3rd Coordinate */
				
					/*Writing the 3rd Coordinate*/
					ier = cg_coord_write(index_file,index_base,index_zone,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype,  
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name3, 
						((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z),&index_coord);
					CGNSErrorCheck(ier);
				} 
			}/* Loop over the Grid Coordinate Nodes ends here*/  
			/*##########################################FLOW SOLUTION WRITTING####################################*/
			ier=cg_sol_write(index_file,index_base,index_zone,floname,location,&sol_index);
			CGNSErrorCheck(ier);
			
			if((CgnsRoot.CgnsBase+(index_base-1))->BaseParameters->flagRadiation==0)
			{
				CGNSErrorCheck(ier);
				if(CgnsRoot.CgnsBase->flagMPFV==0)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->u!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dau,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->u,&solarr_index);	
							CGNSErrorCheck(ier);
 //ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daute,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->utemp,&solarr_index);                                                        CGNSErrorCheck(ier);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daum,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->umass,&solarr_index);	
						CGNSErrorCheck(ier);
					}
					
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->v!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dav,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->v,&solarr_index);	
						CGNSErrorCheck(ier);
 //ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,davte,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vtemp,&solarr_index);                                                        CGNSErrorCheck(ier);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,davm,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vmass,&solarr_index);	
						CGNSErrorCheck(ier);
					}
					
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->w!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daw,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->w,&solarr_index);	
						CGNSErrorCheck(ier);
 //ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dawte,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->wtemp,&solarr_index);                                                        CGNSErrorCheck(ier);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dawm,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->wmass,&solarr_index);	
						CGNSErrorCheck(ier);
					}
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->P!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dap,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->P,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->H!=NULL &&(CgnsRoot.CgnsBase+(index_base-1))->flagGroundwater==1)
				{	
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dah,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->H,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if((CgnsRoot.CgnsBase+(index_base-1))->flagSaturation==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->WaterContent!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dawa,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->WaterContent,&solarr_index);	
						CGNSErrorCheck(ier);
					}
				
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->UnSatPermeability!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,rel,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->UnSatPermeability,&solarr_index);	
						CGNSErrorCheck(ier);
					}
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Flow_Potential!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Flow_Potential",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Flow_Potential,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->T!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dat,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->T,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->enthalpy!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"enthalpy",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->enthalpy,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Latent_Heat_Content!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Latent_Heat",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Latent_Heat_Content,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->f!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"liquid fraction",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->f,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MixFrac!=NULL) //Shyam
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,MixtureFraction,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MixFrac,&solarr_index);
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->FuelMassFrac!=NULL) //Shyam
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,FuelMassFraction,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->FuelMassFrac,&solarr_index);
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->OxiMassFrac!=NULL) //Shyam
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,OxiMassFrac,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->OxiMassFrac,&solarr_index);
					CGNSErrorCheck(ier);
				}
				/*if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->InMassFrac!=NULL) //Shyam
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,InMassFrac,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->InMassFrac,&solarr_index);
					CGNSErrorCheck(ier);
				}*/
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ProdMassFrac!=NULL) //Shyam
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,ProdMassFrac,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ProdMassFrac,&solarr_index);
					CGNSErrorCheck(ier);
				}
			
			
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagRho==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Rho!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,rho,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Rho,&solarr_index);	
						CGNSErrorCheck(ier);
					}
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagViscoeff==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Viscoeff!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,visc,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Viscoeff,&solarr_index);
						CGNSErrorCheck(ier);
					}
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagDiffcoeff==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Diffcoeff!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,cond,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Diffcoeff,&solarr_index);
						CGNSErrorCheck(ier);
					}
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagSpecificHeat==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->SpecificHeat!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,spheat,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->SpecificHeat,&solarr_index);
						CGNSErrorCheck(ier);
					}
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->flagHeatGeneration==1 && ((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->flagHeat_con_var==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Energysource!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,energysource,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Energysource,&solarr_index);
						CGNSErrorCheck(ier);
					}
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagTurb==1|| ((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagTurb==2 || ((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagTurb==3||((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagTurb==4)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TViscoeff!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,viscTurb,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TViscoeff,&solarr_index);
						CGNSErrorCheck(ier);
					}
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->EffViscoeff!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,effviscTurb,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->EffViscoeff,&solarr_index);
						CGNSErrorCheck(ier);
					}
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TProduction!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,turbProdterm,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TProduction,&solarr_index);
						CGNSErrorCheck(ier);
					}
					if(CgnsRoot.CgnsBase->flagTurbModel==4)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,PhiK,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->PhiK,&solarr_index);
						CGNSErrorCheck(ier);
						
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,PhiEps,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->PhiEps,&solarr_index);
						CGNSErrorCheck(ier);
					}
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vorticity_X!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,vortX,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vorticity_X,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vorticity_Y!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,vortY,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vorticity_Y,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vorticity_Z!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,vortZ,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vorticity_Z,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
			//*********************************** Electrical Parameters *************************//
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Potential!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daphi,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Potential,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoIon!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daCharge,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoIon,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->EfieldX!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daEfx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->EfieldX,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->EfieldY!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daEfy,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->EfieldY,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->EfieldZ!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daEfz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->EfieldZ,&solarr_index);	
					CGNSErrorCheck(ier);
				}

			}
			//*********************************** End of Electrical Parameters *******************//
			
			//Added by sushil.....for multiphase flow problems
			if(CgnsRoot.CgnsBase->flagMPFV==1)
			{
				//WriteCgnsFileMPFV( CgnsFileName, CgnsRoot,index_file, sol_index ,index_base,index_zone, solarr_index);
				for(nophase=0;nophase<nop;nophase++)
				{	

		
					sprintf(tempx,"Multiphasemomentum-X");
					sprintf(tempy,"Multiphasemomentum-Y");
					sprintf(tempz,"Multiphasemomentum-Z");
					sprintf(tempV,"Vartheta");
					sprintf(tempvisc,"Phaseviscosity");
					sprintf(tempden,"Phasedensity");

					sprintf(temp,"%d",nophase);
					strcat(tempx,temp);
		
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasex[nophase]!=NULL)	
					{
			
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,P1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasex[nophase],&solarr_index);	
						CGNSErrorCheck(ier);		
					}
		
					strcat(P2,temp);

					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasey[nophase]!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,P2,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasey[nophase],&solarr_index);	
						CGNSErrorCheck(ier);		
					}
					strcat(P3,temp);
		
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasez[nophase]!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,P3,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasez[nophase],&solarr_index);	
						CGNSErrorCheck(ier);		
					}
		
					tempx1=  strcat(tempV,temp);
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->VarthetaPhase[nophase]!=NULL)
					{
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,tempx1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->VarthetaPhase[nophase],&solarr_index);	
						CGNSErrorCheck(ier);		
					}
					tempx1=  strcat(tempden,temp);
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagRhoPhase[nophase]==1)
					{	
			
						if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoofPhase[nophase]!=NULL)
						{
				
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,tempx1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoofPhase[nophase],&solarr_index);	
							CGNSErrorCheck(ier);
				
						}
			
					}	
		
					tempx1=  strcat(tempvisc,temp);
		
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagViscoeffPhase[nophase]==1)	
					{
						if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ViscosityofPhase[nophase]!=NULL)
						{
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,tempx1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ViscosityofPhase[nophase],&solarr_index);	
							CGNSErrorCheck(ier);
						}
					}	
				}	
		
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumx!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumx,&solarr_index);	
					CGNSErrorCheck(ier);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmxs,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumxstar,&solarr_index);	
					CGNSErrorCheck(ier);
				}
	
		
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumy!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmy,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumy,&solarr_index);	
					CGNSErrorCheck(ier);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmys,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumystar,&solarr_index);	
					CGNSErrorCheck(ier);
				}
		
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumz!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumz,&solarr_index);	
					CGNSErrorCheck(ier);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmzs,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumzstar,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				
			}
			//End of multiphase flow loop..........		
			
			//Added by Adhi+shash  .....for Disperse phase flow problems
			if((CgnsRoot.CgnsBase->flagDroplet==1)||(CgnsRoot.CgnsBase->flagParticle==1))
			{
				WriteCgnsFileDisperse( CgnsFileName,CgnsRoot,index_file, sol_index ,index_base,index_zone, solarr_index);
			}
		
			//End of Disperse phase flow loop
			//*********************************** LMMHD Parameters *************************//
			if(CgnsRoot.CgnsBase->flagLMMHD==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->PotentialMHD!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daphiMHD,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->PotentialMHD,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldX!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daMfx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldX,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldY!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daMfy,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldY,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldZ!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daMfz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldZ,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityX!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daCurrentx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityX,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityY!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daCurrenty,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityY,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityZ!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daCurrentz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityZ,&solarr_index);	
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldX!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daTotMfx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldX,&solarr_index);	
					CGNSErrorCheck(ier);//rajendra
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldY!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daTotMfy,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldY,&solarr_index);	
					CGNSErrorCheck(ier);//rajendra
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldZ!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,daTotMfz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldZ,&solarr_index);	
					CGNSErrorCheck(ier);//rajendra
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Psi!=NULL)
				{
					cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,dapsi,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Psi,&solarr_index);	
					CGNSErrorCheck(ier);//rajendra
				}
			}
				
			//***********************************End of LMMHD Parameters *************************//
			for(scal_ind=0;scal_ind<10;scal_ind++)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[scal_ind]!=NULL)
				{
					switch(scal_ind)
					{
						case 0:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[0],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 1:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s2,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[1],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 2:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s3,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[2],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 3:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s4,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[3],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 4:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s5,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[4],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 5:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s6,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[5],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 6:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s7,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[6],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 7:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s8,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[7],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 8:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s9,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[8],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
						
						case 9:
							ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,s10,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[9],&solarr_index);
							CGNSErrorCheck(ier);
						
							break;
					}
						
				}
			}
				
			for(scalCoeff_ind=0;scalCoeff_ind<10;scalCoeff_ind++)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagScalcoeff[scalCoeff_ind]==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[scalCoeff_ind]!=NULL)
					{
						switch(scalCoeff_ind)
						{
							case 0:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[0],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 1:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff2,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[1],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 2:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff3,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[2],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 3:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff4,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[3],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 4:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff5,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[4],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 5:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff6,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[5],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 6:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff7,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[6],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 7:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff8,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[7],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 8:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff9,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[8],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
		
							case 9:
								ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,scoeff10,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[9],&solarr_index);
								CGNSErrorCheck(ier);
		
								break;
						}
					}
				}
			}
				
			
//#######################################GRID CONNECTIVITY WRITTING#########################################################
				
			/* ..... Loop over the Interface nodes for writting the data of all the interfaces in a zone */
			for(index_iface=1;index_iface<=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
			{
				/* ....... Getting the info for the Interface */ 
				ier = cg_1to1_write(index_file,index_base,index_zone,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Interface_Name, 
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Adjoining_Zone_Name
								,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range,
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range
								,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform,&interface_ind);  
				CGNSErrorCheck(ier);	       
			} 
			if(CgnsRoot.CgnsBase->flagMPFV==1)//Added by sushil....for multiphase flow.
			  {printf("\ncatch me if u can\n");
			 
				//#######################################BOUNDARY NODDE WRITTING#########################################################
				nbc= ((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))-> Number_of_BC_Patches;
				//sizeptr=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size;
				//. Loop over the Boundary Patches for getting info 
				for(index_bc = 1;index_bc<=nbc;index_bc++)
				{
						
					sizeptr=(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->Point_Range;
					dir=(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->Direction;
					switch(dir)
					{
						case 1:
							locab=5;
							size_patch[0]=(*(sizeptr+4))*(*(sizeptr+5));
							size1=*(sizeptr+4);
							size2=*(sizeptr+5);	
							break;
						
						case 2:
							locab=5;
							size_patch[0]=(*(sizeptr+4))*(*(sizeptr+5));//(*(sizeptr+3))*(*(sizeptr+4))*(*(sizeptr+5));//(*(sizeptr+1)-1)*(*(sizeptr+2)-1);
							size1=*(sizeptr+4);
							size2=*(sizeptr+5);
							break;
						
						case 3:
							locab=6;
							size_patch[0]=(*(sizeptr+3))*(*(sizeptr+5));//*/(*(sizeptr+3))*(*(sizeptr+4))*(*(sizeptr+5));//(*(sizeptr)-1)*(*(sizeptr+2)-1);
							size1=*(sizeptr+3);
							size2=*(sizeptr+5);
							break;
						
						case 4:
							locab=6;
							size_patch[0]=(*(sizeptr+3))*(*(sizeptr+5));//(*(sizeptr+3))*(*(sizeptr+4))*(*(sizeptr+5));//(*(sizeptr)-1)*(*(sizeptr+2)-1);
							size1=*(sizeptr+3);
							size2=*(sizeptr+5);
							break;
						
						case 5:
							locab=7;
							size_patch[0]=(*(sizeptr+3))*(*(sizeptr+4));//(*(sizeptr+3))*(*(sizeptr+4))*(*(sizeptr+5));//(*(sizeptr)-1)*(*(sizeptr+1)-1);
							size1=*(sizeptr+3);
							size2=*(sizeptr+4);
							break;
						
						case 6:
							locab=7;
							size_patch[0]=(*(sizeptr+3))*(*(sizeptr+4));//(*(sizeptr+3))*(*(sizeptr+4))*(*(sizeptr+5));//(*(sizeptr)-1)*(*(sizeptr+1)-1);
							size1=*(sizeptr+3);
							size2=*(sizeptr+4);
							break;
					}	 
					ier = cg_boco_write(index_file,index_base,index_zone,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->BCName,
							(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->Boundary_Tag,ptsettype,npnts,
							(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->Point_Range,&indexbc);
					CGNSErrorCheck(ier); 
					ier=cg_boco_normal_write(index_file,index_base,index_zone,indexbc,norind,0,ndt,&nlist);	
					CGNSErrorCheck(ier); 
					ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"end");
					CGNSErrorCheck(ier); 
					ier=cg_gridlocation_write(locab);
					CGNSErrorCheck(ier); 	
					
					ier=cg_dataset_write(index_file,index_base,index_zone,indexbc,dset,type,&index_dset);
					CGNSErrorCheck(ier); 	
					ier=cg_bcdata_write(index_file,index_base,index_zone,indexbc,index_dset,type1);
					CGNSErrorCheck(ier);
					ier=cg_bcdata_write(index_file,index_base,index_zone,indexbc,index_dset,type2);
					CGNSErrorCheck(ier);
						
					count=0;
					ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
					CGNSErrorCheck(ier); 
					
					solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->u);
					ier=cg_array_write(dau,soldtype,datadim,size_patch,solution);
					CGNSErrorCheck(ier);
					
					
					solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->v);
					ier=cg_array_write(dav,soldtype,datadim,size_patch,solution);
					CGNSErrorCheck(ier);
					
					
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->w!=NULL)
					{
						solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->w);
						ier=cg_array_write(daw,soldtype,datadim,size_patch,solution);
						CGNSErrorCheck(ier);
					}
			
					if(CgnsRoot.CgnsBase->flagGroundwater==0)
					{
						solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->P);
						ier=cg_array_write(dap,soldtype,datadim,size_patch,solution);
						CGNSErrorCheck(ier);
					}
				
					if((CgnsRoot.CgnsBase+(index_base-1))->flagGroundwater==1)
					{
						solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->H);
						ier=cg_array_write(dah,soldtype,datadim,size_patch,solution);
						CGNSErrorCheck(ier);
			
			
						if((CgnsRoot.CgnsBase+(index_base-1))->flagSaturation==1)
						{
							solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->UnSatPermeability);
							ier=cg_array_write(rel,soldtype,datadim,size_patch,solution);
							CGNSErrorCheck(ier);
						}
					}
						
					if(CgnsRoot.CgnsBase->flagTurbComb==0) //Shyam
					{
						if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->T!=NULL)
						{
							solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->T);
							ier=cg_array_write(dat,soldtype,datadim,size_patch,solution);
							CGNSErrorCheck(ier);
							ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
							CGNSErrorCheck(ier);
							solution=&(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradT));
							ier=cg_array_write(fluxt,soldtype,datadim,fluxsize,solution);
							CGNSErrorCheck(ier);
							solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradT);
							ier=cg_array_write(gradt,soldtype,datadim,size_patch,solution);
							CGNSErrorCheck(ier);
						}
					}
					
					for(scal_ind=0;scal_ind<10;scal_ind++)
					{
						if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[scal_ind]!=NULL)
						{
							switch(scal_ind)
							{
								case 0:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[0]);
									ier=cg_array_write(s1,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs1,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[0]);
									ier=cg_array_write(grads1,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 1:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[1]);
									ier=cg_array_write(s2,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs2,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[1]);
									ier=cg_array_write(grads2,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 2:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[2]);
									ier=cg_array_write(s3,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs3,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[2]);
									ier=cg_array_write(grads3,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 3:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[3]);
									ier=cg_array_write(s4,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs4,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[3]);
									ier=cg_array_write(grads4,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 4:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[4]);
									ier=cg_array_write(s5,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs5,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[4]);
									ier=cg_array_write(grads5,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 5:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[5]);
									ier=cg_array_write(s6,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs6,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[5]);
									ier=cg_array_write(grads6,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 6:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[6]);
									ier=cg_array_write(s7,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs7,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[6]);
									ier=cg_array_write(grads7,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 7:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[7]);
									ier=cg_array_write(s8,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs8,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[7]);
									ier=cg_array_write(grads8,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 8:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[8]);
									ier=cg_array_write(s9,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs9,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[8]);
									ier=cg_array_write(grads9,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
								
								case 9:
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Dirichlet,"end");
									CGNSErrorCheck(ier); 
									solution=ArrayDuplicator( size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->Scal[9]);
									ier=cg_array_write(s10,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									ier=cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",indexbc,"BCDataSet_t",1,"BCData_t",Neumann,"end");
									CGNSErrorCheck(ier); 
									solution=(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->SumNorGradScal+scal_ind));
									ier=cg_array_write(Fluxs10,soldtype,datadim,fluxsize,solution);
									CGNSErrorCheck(ier);
									solution=ArrayDuplicator(size1, size2,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+(index_bc-1))->ZSlab)->GradScal[9]);
									ier=cg_array_write(grads10,soldtype,datadim,size_patch,solution);
									CGNSErrorCheck(ier);
									break;
							}
							
						}
					}
					
				}
			}

			if((CgnsRoot.CgnsBase+(index_base-1))->BaseParameters->flagRadiation==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RadSolTemp!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"RADTemp",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RadSolTemp,&solarr_index);	
					CGNSErrorCheck(ier);
				}
			}
		}
	}
	
ier=cg_close(index_file);
return (0);
} 




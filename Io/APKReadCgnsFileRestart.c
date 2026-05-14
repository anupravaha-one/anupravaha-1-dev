#include "APKReadCgnsFileRestart.h"

/** This Function is called by main function if the user has started his computation form partialy converged solution.
 The partialy converged solutions are stored in file restart.cgns, this function takes the pointer to the file restart.cgns, and the Root variable of the data 
 structureas input.
 
 The function operates in following steps
 1) Opens the file restart.cgns in read mode.
 2) loops over the all the  zones 
 3) Reads the flow solution node of each zone (the partialy converged solution arras of variable)
 4) Stores the partialy converged solutuin arry using the Root variable with proper link to the flowsol node of every zone
 the entitr data structure is not read from this file only flow solutuon node is read rest of the information like grid co ordinates 
 and grid connectivity is read from the original cgns file that user entered at the beginnig.
 */
 
int ReadCgnsFileRestart(char* CgnsFileName, Root CgnsRoot)
{
 
/* --------------- Variable Declaration --------------- */
   DataType_t soldtype;
   int index_file,index_base,index_zone,index_flowsol=1,solarr_index=0;
   /* Dummy Variables */
   int ier,i,nzones;
   int rmin[4],rmax[4],scal_ind;
   char solname[33];
   Zone *zone;
/* --------------- ******************** --------------- */
/* Opening a Data File for writing the Summary of the CGNS File */

/* Opening the CGNS file in READ mode */	
   ier = cg_open(CgnsFileName,MODE_READ,&index_file);
   CGNSErrorCheck(ier);	/* Checking for CGNS Error */
  /* Loop Over Base Nodes for reading the nodes under each Base Node */
        for(index_base=1;index_base<=CgnsRoot.Number_of_bases;index_base++)
        {
			/* Loop Over Zone Nodes for reading the nodes under each Zone Node */ 
		nzones=CgnsRoot.CgnsBase->Number_of_zones;
		for(index_zone=1;index_zone<=nzones;index_zone++)
		{
			solarr_index=0;
			//reading flow solutions
			
			zone = (CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1);
			
			for(i=1;i<=(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension;i++)
			{
				rmin[i-1] = 1;
				rmax[i-1] = *(zone->Size+i-1)-1;	 
			} 

			if(zone->Flowsol->u!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);	
	
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->u);	
					CGNSErrorCheck(ier);
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->umass);	
				CGNSErrorCheck(ier);
			}
			if(zone->Flowsol->v!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->v);	
				CGNSErrorCheck(ier);
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->vmass);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->w!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->w);	
				CGNSErrorCheck(ier);
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->wmass);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->P!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->P);	
				CGNSErrorCheck(ier);

			}
			
			if(zone->Flowsol->Flow_Potential!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Flow_Potential);	
				CGNSErrorCheck(ier);

			}
			
			if(zone->Flowsol->T!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->T);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->enthalpy!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->enthalpy);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->Latent_Heat_Content!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Latent_Heat_Content);
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->f!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->f);
				CGNSErrorCheck(ier);

			}
			
			//Reading combustion variables
			if(zone->Flowsol->MixFrac!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->MixFrac);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->FuelMassFrac!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->FuelMassFrac);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->OxiMassFrac!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->OxiMassFrac);	
				CGNSErrorCheck(ier);

			}
			if(zone->Flowsol->ProdMassFrac!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->ProdMassFrac);	
				CGNSErrorCheck(ier);

			}
			
			
			if(zone->Flowsol->flagRho==1)
			{
				if(zone->Flowsol->Rho!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Rho);
					CGNSErrorCheck(ier);
				}
			}
			
			if(zone->Flowsol->flagViscoeff==1)
			{
				if(zone->Flowsol->Viscoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Viscoeff);
					CGNSErrorCheck(ier);
				}
			}
			
			if(zone->Flowsol->flagDiffcoeff==1)
			{
				if(zone->Flowsol->Diffcoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Diffcoeff);
					CGNSErrorCheck(ier);
				}
			}
			
			if(zone->Flowsol->flagSpecificHeat==1)
			{
				if(zone->Flowsol->SpecificHeat!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->SpecificHeat);
					CGNSErrorCheck(ier);
				}
			}
			
			if(zone->Flowsol->Energysource!=NULL && (zone->flagHeatGeneration==1 && zone->flagHeat_con_var==1))
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Energysource);
				CGNSErrorCheck(ier);
			}
			
			//Writing turbulence variables
			if(zone->Flowsol->flagTurb != 0)
			{
				if(zone->Flowsol->Tcoeff->TViscoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Tcoeff->TViscoeff);
					CGNSErrorCheck(ier);
				}
				if(zone->Flowsol->Tcoeff->EffViscoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Tcoeff->EffViscoeff);
					CGNSErrorCheck(ier);
				}
				if(zone->Flowsol->Tcoeff->TProduction!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Tcoeff->TProduction);
					CGNSErrorCheck(ier);
				}
				if(CgnsRoot.CgnsBase->flagTurbModel==4)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Tcoeff->PhiK);
					CGNSErrorCheck(ier);

					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Tcoeff->PhiEps);
					CGNSErrorCheck(ier);
				}
			}
			
			//Reading Vorticity components
			if(zone->Flowsol->vorticity_X!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->vorticity_X);
				CGNSErrorCheck(ier);
			}
			if(zone->Flowsol->vorticity_Y!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->vorticity_Y);
				CGNSErrorCheck(ier);
			}
			if(zone->Flowsol->vorticity_Z!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->vorticity_Z);
				CGNSErrorCheck(ier);
			}
			
			//*********************************** Electrical Parameters *************************//
			if(zone->Flowsol->Potential!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Potential);	
				CGNSErrorCheck(ier);
				
			}
			if(zone->Flowsol->RhoIon!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->RhoIon);	
				CGNSErrorCheck(ier);

			}
			
			//EfieldX
			if(zone->Flowsol->EfieldX!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->EfieldX);	
				CGNSErrorCheck(ier);
			}
				
			//EfieldY
			if(zone->Flowsol->EfieldY!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->EfieldY);	
				CGNSErrorCheck(ier);
			}
				
			//EfieldZ
			if(zone->Flowsol->EfieldZ!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->EfieldZ);	
				CGNSErrorCheck(ier);
			}
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->PotentialMHD!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->PotentialMHD);	
				CGNSErrorCheck(ier);
			}
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldX!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldX);	
				CGNSErrorCheck(ier);
			}

			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldY!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldY);	
				CGNSErrorCheck(ier);
			}

			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldZ!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MfieldZ);	
				CGNSErrorCheck(ier);
			}
		
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityX!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityX);	
				CGNSErrorCheck(ier);
			}

			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityY!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityY);	
				CGNSErrorCheck(ier);
			}

			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityZ!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ElecurrentdensityZ);	
				CGNSErrorCheck(ier);
			}
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldX!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldX);	
				CGNSErrorCheck(ier);
			}
				
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldY!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldY);	
				CGNSErrorCheck(ier);
			}
				
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldZ!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->TotalMfieldZ);	
				CGNSErrorCheck(ier);
			}

			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Psi!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Psi);	
				CGNSErrorCheck(ier);
			}
//*********************************** End of Electrical Parameters *******************//

			for(scal_ind=0;scal_ind<10;scal_ind++)
			{
				if(zone->Flowsol->Scal[scal_ind]!=NULL)
				{
						solarr_index=solarr_index+1;
						ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
						printf("\n %s",solname);
						ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Scal[scal_ind]);
						CGNSErrorCheck(ier);
				}
			}
			
			for(scal_ind=0;scal_ind<10;scal_ind++)
			{
				if(zone->Flowsol->flagScalcoeff[scal_ind]==1)
				{
					if(zone->Flowsol->Scalcoeff[scal_ind]!=NULL)
					{
							solarr_index=solarr_index+1;
							ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
							printf("\n %s",solname);
							ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,zone->Flowsol->Scalcoeff[scal_ind]);
							CGNSErrorCheck(ier);
					}
				}
			}
		}
	}
  ier=cg_close(index_file); 
  
  return(0);
}   

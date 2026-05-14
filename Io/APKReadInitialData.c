#include "APKReadInitialData.h"
/** This Function is called by main function if the user has started his computation form partialy converged solution.
 The partialy converged solutions are stored in file restart.cgns, this function takes the pointer to the file restart.cgns, and the Root variable of the data
 structureas input.

 The function operates in following steps
 1) Opens the file restart.cgns in read mode.
 2) loops over the all the  zones
 3) Reads the flow solution node of each zone (the partialy converged solution arras of variable)
 4) Stores the partialy converged solutuion array using the Root variable with proper link to the flowsol node of every zone
 the entitr data structure is not read from this file only flow solutuon node is read rest of the information like grid co ordinates
 and grid connectivity is read from the original cgns file that user entered at the beginnig.
 */

int ReadInitialData(char* CgnsFileName, Root CgnsRoot)
{

/* --------------- Variable Declaration --------------- */
   DataType_t soldtype;

   int index_file,index_base,index_zone,index_flowsol=1,solarr_index=0;
   /* Dummy Variables */

   int ier,i,nzones;
   int rmin[4],rmax[4];

   char solname[33];

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
			
			for(i=1;i<=(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension;i++)
			{
				rmin[i-1] = 1;
				rmax[i-1] = *(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+i-1)-1;	 
			} 

			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->u!=NULL)
			{
			solarr_index=solarr_index+1;
			ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
			printf("\n %s",solname);	

			ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->u);	
				CGNSErrorCheck(ier);
			solarr_index=solarr_index+1;
			ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
			printf("\n %s",solname);
			ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->umass);	
			CGNSErrorCheck(ier);
			}
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->v!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->v);	
				CGNSErrorCheck(ier);
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->vmass);	
				CGNSErrorCheck(ier);

			}
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->w!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->w);	
				CGNSErrorCheck(ier);
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->wmass);	
				CGNSErrorCheck(ier);

			}
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->P!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf(" %s\n",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->P);	
				CGNSErrorCheck(ier);

			}
			
			/*if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->T!=NULL)
			{
				solarr_index=solarr_index+1;
				ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
				printf("\n %s",solname);
				ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->T);	
				CGNSErrorCheck(ier);

			}
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagRho==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Rho!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Rho);
					CGNSErrorCheck(ier);
				}
			}
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagViscoeff==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Viscoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Viscoeff);
					CGNSErrorCheck(ier);
				}
			}
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagDiffcoeff==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Diffcoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Diffcoeff);
					CGNSErrorCheck(ier);
				}
			}
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagSpecificHeat==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->SpecificHeat!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->SpecificHeat);
					CGNSErrorCheck(ier);
				}
			}*/
			
			/*if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagTurb==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TViscoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TViscoeff);
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->EffViscoeff!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->EffViscoeff);
					CGNSErrorCheck(ier);
				}
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TProduction!=NULL)
				{
					solarr_index=solarr_index+1;
					ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
					printf("\n %s",solname);
					ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tcoeff->TProduction);
					CGNSErrorCheck(ier);
				}
			}*/
			
			/*for(scal_ind=0;scal_ind<10;scal_ind++)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[scal_ind]!=NULL)
				{
						solarr_index=solarr_index+1;
						ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
						printf("\n %s",solname);
						ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scal[scal_ind]);
						CGNSErrorCheck(ier);
				}
			}
			
			for(scal_ind=0;scal_ind<10;scal_ind++)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagScalcoeff[scal_ind]==1)
				{
					if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[scal_ind]!=NULL)
					{
							solarr_index=solarr_index+1;
							ier=cg_field_info(index_file,index_base,index_zone,index_flowsol,solarr_index,&soldtype,solname);
							printf("\n %s",solname);
							ier=cg_field_read(index_file,index_base,index_zone,index_flowsol,solname,soldtype,rmin,rmax,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Scalcoeff[scal_ind]);
							CGNSErrorCheck(ier);
					}
				}
			}*/
		}
	}
  ier=cg_close(index_file); 
  
  return(0);
}   

#include "APKWriteCgnsFile_vertex.h"

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

int WriteCgnsFile_vertex(char* CgnsFileName,Root CgnsRoot)
{
 
/* --------------- Variable Declaration --------------- */
	//BCType_t type=8;
	//BCDataType_t type1=2,type2=3;
	DataType_t soldtype=4;
	//DataType_t ndt=2;
	int index_file,index_base,index_zone,index_grid,index_coord;
	/* Dummy Variables */
	
	int ier,base_ind,zone_ind,grid_ind,no;
	//double *solution;
	//int npnts=2,norind[]={-1,0,0};
	//int *sizeptr,size_patch[1],dir;
	//int datadim=1;
	int fluxsize[1];
	//int size1,size2;
	//PointSetType_t ptsettype=4;
	//sushil
	int nophase,nop;
	char *temp;
	char *tempx,*tempy,*tempz;
	char *P1,*P2,*P3,*tempx1;
	temp=(char *)malloc(40*sizeof(char));
	tempx=(char *)malloc(40*sizeof(char));
	tempy=(char *)malloc(40*sizeof(char));
	tempz=(char *)malloc(40*sizeof(char));
	
	char tempV[]="Vartheta";
	P1 = tempx;
	P2 =tempy;
	P3 =tempz;
	nop = CgnsRoot.CgnsBase->nphase;
	//end
	///This location has been modified from Cellcenter to vertex for post processing compatability purpose.  Right now this is in testing mode.
	GridLocation_t location=2;
	int i,Coor_size;
	Zone *zone;
	Base *base;



	int sol_index,solarr_index;
	char floname[]="Solutions";

	/* Opening the CGNS file in WRITE mode */	
	fluxsize[0]=1;
	ier = cg_open(CgnsFileName,MODE_WRITE,&index_file);
	CGNSErrorCheck(ier);	/* Checking for CGNS Error */
	

	no=(CgnsRoot.CgnsBase)->Number_of_zones;

	/* Loop Over Base Nodes for reading the nodes under each Base Node */
	for(index_base=1;index_base<=CgnsRoot.Number_of_bases;index_base++)
	{
		base = CgnsRoot.CgnsBase+(index_base-1);
		ier = cg_base_write(index_file,base->basename,base->Cell_Dimension,base->Physical_Dimension,&base_ind);
		CGNSErrorCheck(ier);
		
		/* Loop Over Zone Nodes for reading the nodes under each Zone Node */ 
		for(index_zone=1;index_zone<=base->Number_of_zones;index_zone++)
		{
			zone = ((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1));

			/*writing the Zone Name & Zone Size */
			ier =cg_zone_write(index_file,index_base,zone->zonename,zone->Size,zone->Zonetype,&zone_ind);
			CGNSErrorCheck(ier);
		
			/*##########################################GRID CO-ORDINATES WRITTING####################################*/
			/* ... Loop over all the Grid Coordinate Nodes for reading Grid Coordinates */
			for(index_grid=1;index_grid<=zone->Number_of_Gridcoor_nodes;index_grid++)
			{
				ier=cg_grid_write(index_file,index_base,index_zone,(zone->CgnsGrid+index_grid-1)->GridcoordinateName,&grid_ind);
				CGNSErrorCheck(ier);
					
				/* .......Writing the 1st Coordinate */      
				ier = cg_coord_write(index_file,index_base,index_zone,(zone->CgnsGrid+index_grid-1)->Datatype,(zone->CgnsGrid+(index_grid-1))->Array_Name1,((zone->CgnsGrid+(index_grid-1))->X),&index_coord);
				CGNSErrorCheck(ier);
				
				if(((zone->CgnsGrid+index_grid-1)->Number_of_Coordinates) >=2)
				{
					/* ....... 2nd Coordinate */
					/* .......Writting the 2nd Coordinate */      
					ier = cg_coord_write(index_file,index_base,index_zone,(zone->CgnsGrid+index_grid-1)->Datatype,(zone->CgnsGrid+(index_grid-1))->Array_Name2,((zone->CgnsGrid+(index_grid-1))->Y),&index_coord);
					CGNSErrorCheck(ier); 
				}
				
				if(((zone->CgnsGrid+index_grid-1)->Number_of_Coordinates) >=3)
				{
					/*3rd Coordinate */
					/*Writing the 3rd Coordinate*/
					ier = cg_coord_write(index_file,index_base,index_zone,(zone->CgnsGrid+index_grid-1)->Datatype,(zone->CgnsGrid+(index_grid-1))->Array_Name3,((zone->CgnsGrid+(index_grid-1))->Z),&index_coord);
					CGNSErrorCheck(ier);
				}
			}/* Loop over the Grid Coordinate Nodes ends here*/
			
			
	/*##########################################FLOW SOLUTION WRITTING####################################*/
			ier=cg_sol_write(index_file,index_base,index_zone,floname,location,&sol_index);
			CGNSErrorCheck(ier);
			
			Coor_size = 1;
			
			for(i=0;i<3;i++)
			Coor_size = Coor_size*(*(zone->Size+i));
			zone->Flowsol->VAR_VERTEX = (double *)malloc(sizeof(double)*Coor_size);
			
			if(zone->Flowsol->u!=NULL)
			{
				Variable_AT_Vertex(base,zone,zone->Flowsol->u,1);
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"U velcity",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
			}
			
			if(zone->Flowsol->v!=NULL)
			{
				Variable_AT_Vertex(base,zone,zone->Flowsol->v,2);
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"V Velocity",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
			}
			
			if(zone->Flowsol->w!=NULL)
			{
				Variable_AT_Vertex(base,zone,zone->Flowsol->w,3);
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"W velocity",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
			}
			
			if(zone->Flowsol->T!=NULL)
			{
				Variable_AT_Vertex(base,zone,zone->Flowsol->T,6);
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Temperature",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
			}
			if(zone->Flowsol->Momentumx!=NULL)
			{
				Variable_AT_Vertex(base,zone,zone->Flowsol->Momentumx,31);
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Global Momentum_X",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
			}
			
			if(zone->Flowsol->Momentumy!=NULL)
			{
				Variable_AT_Vertex(base,zone,zone->Flowsol->Momentumy,32);
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Global Momentum_Y",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
			}
			if(CgnsRoot.CgnsBase->BaseParameters->flag2d == 0)
			{
				if(zone->Flowsol->Momentumz!=NULL)
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->Momentumz,33);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Global Momentum_Z",zone->Flowsol->VAR_VERTEX,&solarr_index);
					CGNSErrorCheck(ier);
				}
			}
			for(nophase=0;nophase<nop;nophase++)
			{	
				sprintf(tempx,"Multiphasemomentum-X");
				sprintf(tempy,"Multiphasemomentum-Y");
				sprintf(tempz,"Multiphasemomentum-Z");
				sprintf(tempV,"Vartheta");
				sprintf(temp,"%d",nophase);

				strcat(tempx,temp);
				if(zone->Flowsol->Momentumofphasex[nophase]!=NULL )	
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->Momentumofphasex[nophase],41+3*nophase);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,P1,zone->Flowsol->VAR_VERTEX,&solarr_index);
					CGNSErrorCheck(ier);	
				}
				strcat(P2,temp);
				if(zone->Flowsol->Momentumofphasey[nophase]!=NULL)	
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->Momentumofphasey[nophase],42+3*nophase);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,P2,zone->Flowsol->VAR_VERTEX,&solarr_index);
					CGNSErrorCheck(ier);	
				}
				if(CgnsRoot.CgnsBase->BaseParameters->flag2d == 0)
				{
					strcat(P3,temp);	
					if(zone->Flowsol->Momentumofphasez[nophase]!=NULL)	
					{
						Variable_AT_Vertex(base,zone,zone->Flowsol->Momentumofphasez[nophase],43+3*nophase);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,P3,zone->Flowsol->VAR_VERTEX,&solarr_index);
						CGNSErrorCheck(ier);	
					}
				}
				tempx1=  strcat(tempV,temp);	
			
				if(zone->Flowsol->VarthetaPhase[nophase]!=NULL)
				{	
			
					Variable_AT_Vertex(base,zone,zone->Flowsol->VarthetaPhase[nophase],34+nophase);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempV,zone->Flowsol->VAR_VERTEX,&solarr_index);
					CGNSErrorCheck(ier);
				}
				
				
			
			}
			if((base->flagParticle == 1) || (base->flagDroplet == 1))
			{
				if(zone->Flowsol->U_Disperse!=NULL)
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->U_Disperse,62);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Disperse Phase velocity-X",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
				}
			
				if(zone->Flowsol->V_Disperse!=NULL)
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->V_Disperse,63);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Disperse Phase velocity-Y",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
				}
			
				if(zone->Flowsol->W_Disperse!=NULL)
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->W_Disperse,64);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Disperse Phase velocity-Z",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
				}
				if(zone->Flowsol->VarthetaDisperse!=NULL)
				{
					Variable_AT_Vertex(base,zone,zone->Flowsol->Tdrop,65);
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Disperse Phase Normalized-Volume Fraction",zone->Flowsol->VAR_VERTEX,&solarr_index);
				CGNSErrorCheck(ier);
				}
				if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
				{
					if(zone->Flowsol->RhoDisperse!=NULL)
					{
						Variable_AT_Vertex(base,zone,zone->Flowsol->RhoDisperse,66);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Disperse Phase Density",zone->Flowsol->VAR_VERTEX,&solarr_index);
					CGNSErrorCheck(ier);
					}
					if(zone->Flowsol->DropFuelMassFrac!=NULL)
					{
						Variable_AT_Vertex(base,zone,zone->Flowsol->DropFuelMassFrac,68);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Fuel Mass Fraction",zone->Flowsol->VAR_VERTEX,&solarr_index);
					CGNSErrorCheck(ier);
					}
					if(zone->Flowsol->Tdrop!=NULL)
					{
						Variable_AT_Vertex(base,zone,zone->Flowsol->Tdrop,69);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"Temperature-Droplet",zone->Flowsol->VAR_VERTEX,&solarr_index);
						CGNSErrorCheck(ier);
					}
					if(zone->Flowsol->CPdrop!=NULL)
					{
						Variable_AT_Vertex(base,zone,zone->Flowsol->Tdrop,725);
						ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,"CP-Droplet",zone->Flowsol->VAR_VERTEX,&solarr_index);
						CGNSErrorCheck(ier);
					}
				}
			}

					free(zone->Flowsol->VAR_VERTEX);
		}
	}

	ier=cg_close(index_file);
	return (0);
}




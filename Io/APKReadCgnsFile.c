#include "APKReadCgnsFile.h"

/**
The function ReadCgnsFile() takes the character pointer to the name (with apropreate path of desired CGNS data 
file as input. It reads data from the CGNDS file and simultaneously allocates the memory to the data structure 
members, (Data structure is defined in the file Data_Structures_4.1.h), and stores the read data at its proper 
place. To read the data from a CGNS compatible file Mid Level Library function have been used in this function. 
All the  functions listed above are CGNS Mid-Level Library functions, typically start with a prefix cg_. 
Please reffer CGNS Mid-Level Library (Documennt version 2.3.6 CGNS version 2.2) to get a detailed discription 
of all the library functions. The alogorithem steps have been listed below -
(1) Opens the CGNS file in READ mode.
(2) Reads the  Base node of the CGNS file if cell dimention=2 skips futher 
    commands of ReadCgnsFile() and calls Read_2d_grid().
(3) Else if CellDimention =3  
Allocates the memory to all the Base nodes, Loops over all the base and stores base information. Allocates the 
memory to Zone nodes (sufficent for all the zones). In side the base loop control loops over all the zones 
Reads the zone information and stores it. Inside zone loop control loops over all the Grid coordinate nodes  
Allocates memory to grid co ordinate arrays stores the grid coordinates. Inside zone loop control loops over 
all the Grid connectivity nodes. Allocates necessary memory to store Grid connectivity information. Reads and 
stores the Grid connectivity information. Inside zone loop control loops over all the Boundary condition nodes. 
Allocates necessary memory to store boundary condition  information. Reads and stores the  boundary condition 
information closes the base and zone loop returns the root node of the data structure to main function.
*/

Root ReadCgnsFile(char* CgnsFileName)
{
 
/* --------------- Variable Declaration --------------- */
   
   FILE* sum;  /* File Pointer for the Summary File */
   FILE *var2,*var3;
   
   Root CgnsRoot;
   ZoneType_t ztype;
   DataType_t Dtype,Dtype1,nordt;
   BCType_t bctype;
   
   PointSetType_t ptsettype;
   int index_file,index_base,index_zone,index_grid,index_coord,index_bc,index_iface;
   /* Dummy Variables */
   int nzones,k;
   int ier,i,Coor_size,nbc,test[6],test123;
   int rmin[4],rmax[4];
   char Bond_Cond[33];
   int B[6];
   int A[6];
   void *norm;
   int npnts;
	int loop;
	double scale=1.0;
	
/* --------------- ******************** --------------- */
/* Opening a Data File for writing the Summary of the CGNS File */
   sum = fopen("Settings/Summary","w");
   
   fprintf(sum,"\t\t------------------------------------------------------------------------------- \n");
   fprintf(sum,"\t\t                                  SUMMARY  \n");
   fprintf(sum,"\t\t------------------------------------------------------------------------------- \n\n\n");

   fprintf(sum,"[A] SUMMARY OF THE FILE: \n");
/* Opening the CGNS file in READ mode */	
   ier = cg_open(CgnsFileName,MODE_READ,&index_file);
   CGNSErrorCheck(ier);	/* Checking for CGNS Error */
   fprintf(sum,"1) CGNS File Name : %s\n",CgnsFileName);

/* Reading the version of CGNS Library */
   ier = cg_version(index_file,&CgnsRoot.version);
   CGNSErrorCheck(ier);	/* Checking for CGNS Error */
   fprintf(sum,"2) CGNS MLL Version Number : %f\n",CgnsRoot.version);

/* Reading the Number of Base Nodes present under the Root Node */
   ier = cg_nbases(index_file,&CgnsRoot.Number_of_bases);   
   CGNSErrorCheck(ier);	/* Checking for CGNS Error */
   fprintf(sum,"3) Number of Base Nodes : %d\n",CgnsRoot.Number_of_bases);
/* A valid CGNS File must have atleast one base. The program will quit if
   input CGNS file has no bases. It'll print a warning before quiting */
   if(CgnsRoot.Number_of_bases == 0) {printf("\nNo_DataBase_Warning exiting from ReadCgnsFile\n");exit(1);}
   
//   printf("\nNumber of Base=%d,   Size of Base=%ld \n\n", CgnsRoot.Number_of_bases, sizeof(Base) );

/* Memory Allocation for Base Nodes
   + The Program will quit if sufficient memory is not available */
   CgnsRoot.CgnsBase = (Base *)malloc(CgnsRoot.Number_of_bases*sizeof(Base));
   if(CgnsRoot.CgnsBase == NULL) {printf("\nMemory Warning, exiting from ReadCgnsFile\n");exit(1);}

   fprintf(sum,"\n\n[B] SUMMARY FOR THE BASES:\n");
   
/* Loop Over Base Nodes for reading the nodes under each Base Node */
   for(index_base=1;index_base<=CgnsRoot.Number_of_bases;index_base++)
   {

	fprintf(sum,"\n<<<-------------------- Base Number %d ---------------------->>>\n",index_base); 


	/* Reading Base Name, Cell Dimension & Physical Dimension
	+ bname, cdim & physdim are dummy variables for Base Name,
	Cell Dimension & Physical Dimension resp.
	bname   - Name of the Base Node
	cdim    - Cell Dimension of the base i.e. Dimensionality
		of the cell in the mesh (e.g., 3 for volume cells, 2 for face cells).
	physdim - Physical Dimension of the base i.e. Number of
		coordinates required to define a node position
		(e.g. 1 for 1D, 2 for 2D, 3 for 3D) 
	+ Dummy Variables are used only to avoid lengthy statements */
	ier = cg_base_read(index_file,index_base,(CgnsRoot.CgnsBase+(index_base-1))->basename,&(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension,            
			&(CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension);CGNSErrorCheck(ier);
	
	//Added by chandrasekhar to change phy.dim 2 to 3 as phy dim read by Gridpro is 2
	//Note:This is to be done only when physical dimension read by GRIDPRO is 2.
	if(strcmp(((CgnsRoot.CgnsBase+(index_base-1))->basename),"BASE")== 0)
	{
		//if((CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension==2)
		//	(CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension=3;
		//printf("\nPhy Dim-%d\n",(CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension);
	}
	//End---chandrasekhar
	
	fprintf(sum,"1) Base Name          : %s\n",(CgnsRoot.CgnsBase+(index_base-1))->basename);
	fprintf(sum,"2) Cell Dimension     : %d\n",(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension);
	fprintf(sum,"3) Physical Dimension : %d\n",(CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension);
	//printf(" Cell Dimension     : %d\n",(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension);
	//printf(" Physical Dimension : %d\n",(CgnsRoot.CgnsBase+(index_base-1))->Physical_Dimension);
	
	
	
	/* Reading the number of Zone Nodes present under current Base Node
	+ nzones is a dummy variable for number of zones. It is used only
	to avoid lengthy statements
	+ A valid CGNS file must have atleast one zone node under each base             
	node. The program will quit if any base with no zones is found */
	ier = cg_nzones(index_file,index_base,&nzones);
	CGNSErrorCheck(ier);
	if(nzones == 0) {printf("\nNo_Zone_Warning exiting from ReadCgnsFile\n");exit(1);}
	(CgnsRoot.CgnsBase+(index_base-1))->Number_of_zones = nzones;
	fprintf(sum,"4) Number of Zones    : %d\n",nzones);
	(CgnsRoot.CgnsBase+(index_base-1))->CgnsZone = (Zone *)malloc(nzones*sizeof(Zone));
	if((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone == NULL) {printf("\n01Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
	
	((CgnsRoot.CgnsBase+(index_base-1))->BaseParameters)=(Parameters *)malloc(sizeof(Parameters));
	if((CgnsRoot.CgnsBase+(index_base-1))->BaseParameters == NULL) {printf("\n02Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
	
	if((CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension==3)
	{
		printf("\n\nEnter the scaling factor for the grid in all the directions: \nEnter 1.0 if no scaling is required\nFor ex. enter 0.001 for conversion of m into mm \nRemember that by default the solver is written for m dimenstions\nunless scaling factor is entered it will consider that the units are in meters:  ");
		scanf("%lf",&scale);	
                
                var2=fopen("Settings/coefficients.dat","a");    	
                var3=fopen("Settings/ProgramReadInputValues.dat","a");    
                fprintf(var2,"%lf\n",scale);	        
                fprintf(var3,"%lf \tEnter the scaling factor for the grid in all the directions \n\n",scale);
                fclose(var2);        
                fclose(var3);
                
		/* Updating the Main Database */
		(CgnsRoot.CgnsBase+(index_base-1))->BaseParameters-> flag2d=0;
		(CgnsRoot.CgnsBase+(index_base-1))->BaseParameters->AxisymAzimuthangle=NULL;
		/* Memory Allocation for Zone Nodes
		+ The Program will quit if sufficient memory is not available */
		

		fprintf(sum,"5) Summary of the Zone Data  \n");

		/* Loop Over Zone Nodes for reading the nodes under each Zone Node */ 
		for(index_zone=1;index_zone<=nzones;index_zone++)
		{

			fprintf(sum,"\n\tZone Number %d |----------------------------> \n",index_zone);

			/* ... Reading the Zone Type 
			+ Zone Type can either be structured or Unstructured i.e. it can
				have either structured or Unstructured grid
			+ This program can only read Structured Zone. So if Unstructured
				zone is encountered, the program will quit */    
			ier = cg_zone_type(index_file,index_base,index_zone,&ztype);
			CGNSErrorCheck(ier);

			if(ztype == Unstructured) {printf("\nInvalid Zone type found exiting from ReadCgnsFile\n");exit(1);}
			/* ... Updating the Main Database */
			((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Zonetype = ztype;
			fprintf(sum,"\t\ta) Zone Type                   : %d\n",ztype);    

			/* ... Memory Allocation for Zone Size 
			+ Zone Size is a 2D array of dimension [3 X Cell Dimension] which have Vertex Size,
				Cell Size and Vertex Size Boundary as its three rows. However,for Unstructured
				Grid the dimension of the array is [3 X 1]
			+ Vertex Size - For structured grids, it is the number of vertices in each
					index direction. For unstructured grids it is the total 
					number of vertices in the grid.
			+ Cell Size   - For structured grids, it is the number of cells in each
					index direction. For unstructured grids it is the total 
					number of cells in the grid.
			+ Vertex Size Boundary - If the nodes are sorted between internal & boundary 
					nodes then Vertex Size Boundary is the number of boundary nodes.
					Its value is zero if the nodes are unsorted. For structured grids 
					it is always taken as 0. For more details read the CGNS SIDS
					document (Page 45-46).
			+ The Program will quit if sufficient memory is not available */          

			((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size=(int*)malloc((3*(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension)*sizeof(int));
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size == NULL) {printf("\n03Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ... Reading the Zone Name & Zone Size */
			ier = cg_zone_read(index_file,index_base,index_zone,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->zonename,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size);	
			CGNSErrorCheck(ier);
			fprintf(sum,"\t\tb) Zone Name                   : %s\n",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->zonename);
			fprintf(sum,"\t\tc) Zone Vertex Size            : [%d\t%d\t%d]\n",*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size),
											*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+1),
											*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+2));
			fprintf(sum,"\t\td) Zone Cell Size              : [%d\t%d\t%d]\n",*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+3),
											*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+4),
											*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+5));
			fprintf(sum,"\t\te) Zone Boundary Vertex Size   : [%d\t%d\t%d]\n",*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+6),
											*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+7),
											*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+8));    





			/* ... Reading the Number of Grid Coordinate Nodes under the current Zone node.
			+ Time dependent meshes or moving meshes have time dependent grid coordinates
				So there is a provision of saving grid Coordinates at different time steps 
				or iterations under different Grid Coordinate Nodes. The first node always 
				represents the original or initial Grid Coordinates. 
				For stationary meshes there will only be one Grid Coordinate Node.*/ 

			fprintf(sum,"\t\tf) Details of Grid Coordinates : \n");  
			ier = cg_ngrids(index_file,index_base,index_zone,&((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes);
			CGNSErrorCheck(ier);
			fprintf(sum,"\t\t\t  i) Number of Grid Coordinate Nodes : %d\n",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes);
			/* ... Memory Allocation for Grid Coordinate Nodes
			+ The Program will quit if sufficient memory is not available */
			((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid=(GridCoordinate*)                
				malloc(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes*sizeof(GridCoordinate));                
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid == NULL) {printf("\n04Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
			/* ... Loop over all the Grid Coordinate Nodes for reading Grid Coordinates */
			for(index_grid=1;index_grid<=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes;index_grid++)
			{
			/* ..... Reading the Grid Coordinate Node Name */	      
				ier=cg_grid_read(index_file,index_base,index_zone,index_grid,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->GridcoordinateName);
				CGNSErrorCheck(ier);
				fprintf(sum,"\t\t\t ii) Grid Coordinate Node Name       : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->GridcoordinateName);
			/* ..... Reading the Number of Nodes storing Coordinates under current Grid Coordinate Nodes
				+ This number is equal to the physical dimension */	    
				ier = cg_ncoords(index_file,index_base,index_zone,&((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates));
				CGNSErrorCheck(ier);

				fprintf(sum,"\t\t\tiii) Number of Coordinates           : %d\n",((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates));
			/* ..... Finding the size of the Coordinate Arrays. 
				+ Coor_Size = NI*NJ*NK (For 3D)
				+ Coor_size = NI*NJ    (For 2D)
				+ Coor_size = NI       (For 1D) */
				Coor_size = 1;
				for(i=1;i<=(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension;i++)
				Coor_size = Coor_size*(*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+i-1)) ;

			/* ..... Finding the Maximum & Minimum Index of the Vertices in the Zone 
				+ Mimimum Index will always be 1 */              
				for(i=1;i<=(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension;i++)
				{
				rmin[i-1] = 1;
				rmax[i-1] = *(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Size+i-1);	 
				} 
			/* ....... 1st Coordinate */
				index_coord = 1;
				fprintf(sum,"\t\t\t\t* Details of the 1st Coordinate  :\n");
			/* ....... Reading the Array Name and the Datatype of the Coordinates for the 1st Coordinate */                
				ier=cg_coord_info(index_file,index_base,index_zone,index_coord,&Dtype,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name1);
				fprintf(sum,"\t\t\t\t\t + Name of the Coordinate Array : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name1);
				fprintf(sum,"\t\t\t\t\t + DataType of the Coordinates  : %d\n",Dtype);
				CGNSErrorCheck(ier);     
				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype = Dtype;  

			/* ....... Memory Allocation for the 1st Coordinate Array*/

				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X = (double*)malloc(Coor_size*sizeof(double));
				if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X == NULL) {printf("\n05Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ....... Reading the 1st Coordinate */      
				ier = cg_coord_read(index_file,index_base,index_zone,  
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name1, 
						Dtype,rmin,rmax,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X));
				
				CGNSErrorCheck(ier);

				fprintf(sum,"\t\t\t\t\t + Status (Success if read)     : Success\n");     

				if(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates) >=2)
				{
			/* ....... 2nd Coordinate */
				index_coord = 2;
				fprintf(sum,"\t\t\t\t* Details of the 2nd Coordinate  :\n");
			/* ....... Reading the Array Name and the Datatype of the Coordinates for the 1st Coordinate */                
				ier=cg_coord_info(index_file,index_base,index_zone,index_coord,&Dtype1,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name2);
				fprintf(sum,"\t\t\t\t\t + Name of the Coordinate Array : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name2);
				fprintf(sum,"\t\t\t\t\t + DataType of the Coordinates  : %d\n",Dtype);
				CGNSErrorCheck(ier);

				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype = Dtype1;
			/* ....... Memory Allocation for the 2nd Coordinate Array */

				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y = (double*)malloc(Coor_size*sizeof(double));
				if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y == NULL) {printf("\n06Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ....... Reading the 2nd Coordinate */      
				ier = cg_coord_read(index_file,index_base,index_zone,  
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name2, 
						Dtype,rmin,rmax,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y));
				CGNSErrorCheck(ier); 

				fprintf(sum,"\t\t\t\t\t + Status (Success if read)     : Success\n"); 

				}  

				if(((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates) >=3)
				{
			/* ....... 3rd Coordinate */
				index_coord = 3;
				fprintf(sum,"\t\t\t\t* Details of the 3rd Coordinate  :\n");
			/* ....... Reading the Array Name and the Datatype of the Coordinates for the 3rd Coordinate */                
				ier=cg_coord_info(index_file,index_base,index_zone,index_coord,&Dtype,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name3);
				fprintf(sum,"\t\t\t\t\t + Name of the Coordinate Array : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name3);
				fprintf(sum,"\t\t\t\t\t + DataType of the Coordinates  : %d\n",Dtype);
				CGNSErrorCheck(ier);     
				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype = Dtype;  

			/* ....... Memory Allocation for the 3rd Coordinate Array*/

				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z = (double*)malloc(Coor_size*sizeof(double));                   
				if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z == NULL) {printf("\n07Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ....... Reading the 3rd Coordinate */      
				ier = cg_coord_read(index_file,index_base,index_zone,  
						(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name3, 
						Dtype,rmin,rmax,((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z));
				CGNSErrorCheck(ier);        
				fprintf(sum,"\t\t\t\t\t + Status (Success if read)     : Success\n"); 	     
				}
				///Scaling the grid according to the user choice
				for(loop=0;loop<Coor_size;loop++)
				{
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X+loop)=*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X+loop)*scale;
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y+loop)=*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y+loop)*scale;
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z+loop)=*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z+loop)*scale;
				} 
			} /* Loop over the Grid Coordinate Nodes ends here */ 



			/* ..... Reading the Number of Interface Nodes  */
			ier = cg_n1to1(index_file,index_base,index_zone,&((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces);
			CGNSErrorCheck(ier);
			fprintf(sum,"\t\tg) Number of Interfaces        : %d\n",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces); 
			fprintf(sum,"\t\th) Summary of the Interfaces   : \n");
			/* ..... Memory Allocation for Interface */


			((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface = (Interface*)malloc(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces*sizeof(Interface));
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface == NULL) {printf("\n08Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
			/* ..... Loop over the Interface nodes for reading the data of all the interfaces in a zone */



			for(index_iface=1;index_iface<=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
			{

				/* ....... Memory Allocation for Point Range */
				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range = (int*)malloc(2*(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension*(sizeof(int)));
				if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range == NULL) {printf("\n09Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

				/* ....... Memory Allocation for Point Donor Range */
				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range = (int*)malloc(2*(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension*(sizeof(int)));
				if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range == NULL) {printf("\n10Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

				/* ....... Memory Allocation for Transformation Rule */
				(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform = (int*)malloc((CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension*(sizeof(int)));
				if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform == NULL) {printf("\n11Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

				/* ....... Getting the info for the Interface */ 
				ier = cg_1to1_read(index_file,index_base,index_zone,index_iface,(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Interface_Name, 
												(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Adjoining_Zone_Name,
												(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range,
												(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range,
												(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform);  

				CGNSErrorCheck(ier);	
				fprintf(sum,"\t\t\t * Interface Number %d \n",index_iface);
				fprintf(sum,"\t\t\t\t   i) Name of the Interface      : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Interface_Name);
				fprintf(sum,"\t\t\t\t  ii) Name of the Adjoining Zone : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Adjoining_Zone_Name);      

				fprintf(sum,"\t\t\t\t iii) Point Range                : [%d  %d  %d ] -> [%d  %d  %d ]\n",*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range,
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+1),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+2),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+3),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+4),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+5));
				fprintf(sum,"\t\t\t\t  iv) Point Donor Range          : [%d  %d  %d ] -> [%d  %d  %d ]\n",*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range,
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+1),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+2),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+3),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+4),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+5));
				fprintf(sum,"\t\t\t\t   v) Transformation Rule        : [%d  %d  %d ]\n"               ,*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform,
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform+1),
															*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform+2));

			} 


			/* ..... Reading the Boundary nodes of the zone */

			ier = cg_nbocos(index_file,index_base,index_zone,&nbc);  
			CGNSErrorCheck(ier);
			 //printf("%d....",nbc);
			/* ..... Memory Allocation for Boundary */
			k=6-((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces;
			/*printf("%d..........",k);*/
			((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc = (ZoneBoundary*)malloc(nbc*sizeof(ZoneBoundary));
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc == NULL) {printf("\n12Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
			/* ..... Loop over the Boundary Patches for getting info */

			int No_Ofbc=0;

			for(index_bc = 1;index_bc<=nbc;index_bc++) 
			{
				/* ....... Memory Allocation for Point Range */


				ier = cg_boco_info(index_file,index_base,index_zone,index_bc,Bond_Cond,&bctype,&ptsettype,&npnts,B,B,&nordt,B);
				CGNSErrorCheck(ier); 

				ier=cg_boco_read(index_file,index_base,index_zone,index_bc,A,norm);
				CGNSErrorCheck(ier); 


				for(index_iface=1;index_iface<=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
				{
					test[index_iface-1]=0;
					for(i=0;i<6;i++)
					{


						if( A[i]==*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+i))
							test[index_iface-1]=1+test[index_iface-1];

						else
							test[index_iface-1]=0+test[index_iface-1];
					}


				}     
				test123=0;

				for(index_iface=1;index_iface<=((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
				{	if(test[index_iface-1]==6)
					{
						test123=6;
						break;
					}
					else
					{
						test123=0;	
					}
				}		 




				if(test123==0)
				{ 	

					for(i=0;i<33;i++)
					(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->BCName[i]=Bond_Cond[i];

					(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range = (int*)malloc(2*(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension*(sizeof(int)));
					if((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range == NULL) {printf("\n13Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range)=*A;
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+1)=*(A+1);
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+2)=*(A+2);
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+3)=*(A+3);
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+4)=*(A+4);
					*((((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+5)=*(A+5);
					(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Boundary_Tag=bctype;

					fprintf(sum,"\t\t\t * Boundary Patch Number %d\n",index_bc);
					fprintf(sum,"\t\t\t\t   i) Name of the Boundary Patch : %s\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->BCName);
					fprintf(sum,"\t\t\t\t  ii) Type Of Boundary           : %d\n",(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Boundary_Tag);
					fprintf(sum,"\t\t\t\t iii) Point set Type             : %d\n",ptsettype);
					fprintf(sum,"\t\t\t\t  iv) Point Range                : [%d  %d  %d] -> [%d  %d  %d]\n",*(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range,A[1],A[2],A[3],A[4],A[5]);
					No_Ofbc=No_Ofbc+1;
				}


			}


				((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_BC_Patches = No_Ofbc;
				fprintf(sum,"\t\ti) Number of Zone Boundaries     : %d\n",((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Number_of_BC_Patches);
		}

	
	}
	 else
	 {
	 	(CgnsRoot.CgnsBase+(index_base-1))->Cell_Dimension=3;
		
		read_2d_grid(nzones,CgnsRoot.CgnsBase,index_file,index_base,sum);
	 } 
  
  
  }
  fclose(sum);
ier=cg_close(index_file); 
 return (CgnsRoot);
}
/*##########################CODE FOR 2D to 3D conversion#############################*/

/*!
/detail This function Read_2d_grid() is called by ReadCgnsFile() if the entered grid is a 2-D grid this is called after the base node is read.
input parameters to this function are nomber of zones in the base, the pointer to the structure base, index of the file tht was opened in ReadCgnsFile()
and pointer to summry file. The algorithmic steps for reading a 2-D file are same as a 3-D file. The diffrance between tho functions  ReadCgnsFile() and 
 Read_2d_grid() is that the function Read_2d_grid() apart from reading the data form 2-D CGNS file converts the 2-D data in to 3-D data and then stores it in the data structure.
*/

int read_2d_grid(int nzones, Base *base,int index_file,int index_base,FILE *sum)
{
        FILE *var2,*var3;		
	int dummy_size[6],dummy_pr[4],dummy_pdr[4],dummy_tr[2],dummyA[4],max_ind[2],min_ind[2];
	double *ptrY,*ptrX,*ptrZ,theta,grad;
	char kmin[]="Kmin";
	char kmax[]="Kmax";

	ZoneType_t ztype;
   DataType_t Dtype,Dtype1,nordt;
   BCType_t bctype;
   
   PointSetType_t ptsettype;
   int index_zone,index_grid,index_coord,index_bc,index_iface,domainflag;
   /* Dummy Variables */
   int k,j;
   int ier,i,Coor_size,nbc,test[6],test123;
   int rmin[4],rmax[4];
   char Bond_Cond[33];
   int B[6];
   int A[6];
   void *norm;
   int npnts,loop;
	double scale=1.0;

		fprintf(sum,"5) Summary of the Zone Data  \n");
		
		/* Loop Over Zone Nodes for reading the nodes under each Zone Node */ 
		base->BaseParameters->flag2d=1;
                
                do {
		printf("\n\nYOU HAVE ENTERED A TWO DIMENTIONAL GRID!!\n");
		printf("Please note that this solver is an inherent 3d solver.\nHence your 2d domain is going to change in to a 3d one\n");
		printf("ENTER:  1 if you are solving for a 2-DIMENTIONAL FLOW\n");
		printf("ENTER:  2 if you are solving for an AXI-SYMMETRIC flow(without swirl):");
                scanf("%d", &domainflag);
                }while(!(domainflag==1 || domainflag==2));

		printf("\n\nEnter the scaling factor for the grid in all the directions: \nEnter 1.0 if no scaling is required\nFor ex. enter 0.001 if the mesh is created in millimeters (mm)\nRemember that by default the solver is written for meter (m) dimensions\nUnless scaling factor is entered it will consider that the units are in meters:  ");
		scanf("%lf",&scale);                
                
                var2=fopen("Settings/coefficients.dat","a");    	
                var3=fopen("Settings/ProgramReadInputValues.dat","a");                  
                fprintf(var2,"%d\n",domainflag);                  
                fprintf(var2,"%lf\n",scale);	        
                fprintf(var3,"%d \tEnter 1 for 2-DIMENTIONAL FLOW, 2 for AXI-SYMMETRIC flow \n\n",domainflag);                
                fprintf(var3,"%lf \tEnter the scaling factor for the grid in all the directions \n\n",scale);
                fclose(var2);        
                fclose(var3);
                
		if(domainflag==2)
		{
			printf("YOU HAVE CHOSEN AXI-SYMMETRIC FLOW!!\n\n");
			printf("PLEASE NOTICE THAT Y DIRECTION IN YOUR GRID WILL BE TREATED AS DEFAULT RADIAL DIRECTION :-\n");
			printf("MAKE SURE THAT:- \n"); 
			printf("(1)YOUR GRID IS COMPATIBLE TO ABOVE STATEMENT\n(2)THERE IS NO NEGATIVE Y CO-ORDINATE IN YOUR GRID\n");
			base->BaseParameters->AxisymAzimuthangle=(double*)malloc(sizeof(double));	
		
		}
		else if(domainflag==1)
		{
			printf("\nYOU HAVE CHOSEN 2-DIMENTIONAL FLOW!!\n");
			base->BaseParameters->AxisymAzimuthangle=NULL;
		}
		for(index_zone=1;index_zone<=nzones;index_zone++)
		{

			fprintf(sum,"\n\tZone Number %d |----------------------------> \n",index_zone);

			/* ... Reading the Zone Type 
			+ Zone Type can either be structured or Unstructured i.e. it can
				have either structured or Unstructured grid
			+ This program can only read Structured Zone. So if Unstructured
				zone is encountered, the program will quit */    
			ier = cg_zone_type(index_file,index_base,index_zone,&ztype);
			CGNSErrorCheck(ier);

			if(ztype == Unstructured) {printf("\nInvalid Zone Type found exiting from ReadCgnsFile\n");exit(1);}
			/* ... Updating the Main Database */
			(base->CgnsZone+(index_zone-1))->Zonetype = ztype;
			fprintf(sum,"\t\ta) Zone Type                   : %d\n",ztype);    

			
			(base->CgnsZone+(index_zone-1))->Size=(int*)malloc(9*sizeof(int));
			if((base->CgnsZone+(index_zone-1))->Size == NULL) {printf("\n14Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ... Reading the Zone Name & Zone Size */
			ier = cg_zone_read(index_file,index_base,index_zone,(base->CgnsZone+(index_zone-1))->zonename,dummy_size);	
			CGNSErrorCheck(ier);
			
			*((base->CgnsZone+(index_zone-1))->Size)=dummy_size[0];
			*((base->CgnsZone+(index_zone-1))->Size+1)=dummy_size[1];
			*((base->CgnsZone+(index_zone-1))->Size+2)=2;
			*((base->CgnsZone+(index_zone-1))->Size+3)=dummy_size[2];
			*((base->CgnsZone+(index_zone-1))->Size+4)=dummy_size[3];
			*((base->CgnsZone+(index_zone-1))->Size+5)=1;
			*((base->CgnsZone+(index_zone-1))->Size+6)=dummy_size[4];
			*((base->CgnsZone+(index_zone-1))->Size+7)=dummy_size[5];
			*((base->CgnsZone+(index_zone-1))->Size+8)= 0;    
			fprintf(sum,"\t\tb) Zone Name                   : %s\n",(base->CgnsZone+(index_zone-1))->zonename);
			fprintf(sum,"\t\tc) Zone Vertex Size            : [%d\t%d\t%d]\n",*((base->CgnsZone+(index_zone-1))->Size),
											*((base->CgnsZone+(index_zone-1))->Size+1),
											*((base->CgnsZone+(index_zone-1))->Size+2));
			fprintf(sum,"\t\td) Zone Cell Size              : [%d\t%d\t%d]\n",*((base->CgnsZone+(index_zone-1))->Size+3),
											*((base->CgnsZone+(index_zone-1))->Size+4),
											*((base->CgnsZone+(index_zone-1))->Size+5));
			fprintf(sum,"\t\te) Zone Boundary Vertex Size   : [%d\t%d\t%d]\n",*((base->CgnsZone+(index_zone-1))->Size+6),
											*((base->CgnsZone+(index_zone-1))->Size+7),
											*((base->CgnsZone+(index_zone-1))->Size+8));    
			fprintf(sum,"\t\tf) Details of Grid Coordinates : \n");  
			ier = cg_ngrids(index_file,index_base,index_zone,&(base->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes);
			CGNSErrorCheck(ier);
			fprintf(sum,"\t\t\t  i) Number of Grid Coordinate Nodes : %d\n",(base->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes);
			/* ... Memory Allocation for Grid Coordinate Nodes
			+ The Program will quit if sufficient memory is not available */
			(base->CgnsZone+(index_zone-1))->CgnsGrid=(GridCoordinate*)                
				malloc((base->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes*sizeof(GridCoordinate));                
			if((base->CgnsZone+(index_zone-1))->CgnsGrid == NULL) {printf("\n15Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
			/* ... Loop over all the Grid Coordinate Nodes for reading Grid Coordinates */
			for(index_grid=1;index_grid<=(base->CgnsZone+(index_zone-1))->Number_of_Gridcoor_nodes;index_grid++)
			{
			/* ..... Reading the Grid Coordinate Node Name */	      
				ier=cg_grid_read(index_file,index_base,index_zone,index_grid,((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->GridcoordinateName);
				CGNSErrorCheck(ier);
				fprintf(sum,"\t\t\t ii) Grid Coordinate Node Name       : %s\n",((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->GridcoordinateName);
			/* ..... Reading the Number of Nodes storing Coordinates under current Grid Coordinate Nodes
				+ This number is equal to the physical dimension */	    
				ier = cg_ncoords(index_file,index_base,index_zone,&(((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates));
				CGNSErrorCheck(ier);
				
				//Added by chandrasekhar as number of coord-nodes  read by GridPro is 2.
				if(strcmp((base->basename),"BASE")== 0)
				{
					if(((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates==2)
						((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates=3;
				}
				//End
				
				fprintf(sum,"\t\t\tiii) Number of Coordinates           : %d\n",((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates);
				
			/* ..... Finding the size of the Coordinate Arrays. 
				+ Coor_Size = NI*NJ*NK (For 3D)
				+ Coor_size = NI*NJ    (For 2D)
				+ Coor_size = NI       (For 1D) */
				Coor_size = 1;
				for(i=1;i<= base->Cell_Dimension;i++)
				Coor_size = Coor_size*(*((base->CgnsZone+(index_zone-1))->Size+i-1)) ;

			/* ..... Finding the Maximum & Minimum Index of the Vertices in the Zone 
				+ Mimimum Index will always be 1 */              
				for(i=1;i<=base->Cell_Dimension;i++)
				{
					rmin[i-1] = 1;
					rmax[i-1] = *((base->CgnsZone+(index_zone-1))->Size+i-1);	 
				} 
			/* ....... 1st Coordinate */
				index_coord = 1;
				fprintf(sum,"\t\t\t\t* Details of the 1st Coordinate  :\n");
			/* ....... Reading the Array Name and the Datatype of the Coordinates for the 1st Coordinate */                
				ier=cg_coord_info(index_file,index_base,index_zone,index_coord,&Dtype,((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name1);
				fprintf(sum,"\t\t\t\t\t + Name of the Coordinate Array : %s\n",((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name1);
				fprintf(sum,"\t\t\t\t\t + DataType of the Coordinates  : %d\n",Dtype);
				CGNSErrorCheck(ier);     
				((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Datatype = Dtype;  

			/* ....... Memory Allocation for the 1st Coordinate Array*/

				((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X = (double*)malloc(Coor_size*sizeof(double));
				if(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X == NULL) {printf("\n16Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ....... Reading the 1st Coordinate */      
				ier = cg_coord_read(index_file,index_base,index_zone,  
						((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name1, 
						Dtype,rmin,rmax,((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X);
				CGNSErrorCheck(ier);

				fprintf(sum,"\t\t\t\t\t + Status (Success if read)     : Success\n");     
				
				
			/* ....... 2nd Coordinate */
				index_coord = 2;
				fprintf(sum,"\t\t\t\t* Details of the 2nd Coordinate  :\n");
			/* ....... Reading the Array Name and the Datatype of the Coordinates for the 1st Coordinate */                
				ier=cg_coord_info(index_file,index_base,index_zone,index_coord,&Dtype1,((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name2);
				fprintf(sum,"\t\t\t\t\t + Name of the Coordinate Array : %s\n",((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name2);
				fprintf(sum,"\t\t\t\t\t + DataType of the Coordinates  : %d\n",Dtype);
				CGNSErrorCheck(ier);

				
			/* ....... Memory Allocation for the 2nd Coordinate Array */

				((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y = (double*)malloc(Coor_size*sizeof(double));
				if(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y == NULL) {printf("\n17Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

			/* ....... Reading the 2nd Coordinate */      
				ier = cg_coord_read(index_file,index_base,index_zone,  
						((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name2, 
						Dtype,rmin,rmax,((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y);
				CGNSErrorCheck(ier); 

				fprintf(sum,"\t\t\t\t\t + Status (Success if read)     : Success\n"); 

				///Scaling the grid according to the user choice
				for(loop=0;loop<Coor_size;loop++)
				{
					*(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X+loop)=*(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X+loop)*scale;
					*(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y+loop)=*(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y+loop)*scale;
				} 

				ptrY=((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Y;
				ptrX=((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->X;
				for(j=0;j<dummy_size[1];j++)
				{
					for(i=0;i<dummy_size[0];i++)
					{
						*(ptrX+dummy_size[0]*dummy_size[1]+dummy_size[0]*j+i)=*(ptrX+dummy_size[0]*j+i);	
				 		*(ptrY+dummy_size[0]*dummy_size[1]+dummy_size[0]*j+i)=*(ptrY+dummy_size[0]*j+i);
				 	}
				}
				 ((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z = (double*)malloc(Coor_size*sizeof(double));                   
				if(((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z == NULL) {printf("\n18Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

                                if((((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Number_of_Coordinates >=3) && (base->Physical_Dimension==3))
				{
				/* ....... 3rd Coordinate */
					index_coord = 3;
					fprintf(sum,"\t\t\t\t* Details of the 3rd Coordinate  :\n");
				/* ....... Reading the Array Name and the Datatype of the Coordinates for the 3rd Coordinate */                
					ier=cg_coord_info(index_file,index_base,index_zone,index_coord,&Dtype,((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name3);
					CGNSErrorCheck(ier); 
					char name[]="CoordinateZ";
					strcpy(((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name3,name);
					fprintf(sum,"\t\t\t\t\t + Name of the Coordinate Array : %s\n",((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Array_Name3);
					fprintf(sum,"\t\t\t\t\t + DataType of the Coordinates  : %d\n",Dtype);
				}
				else
				{
					
					char name[]="CoordinateZ";
					strcpy(((base->CgnsZone+(index_zone-1))->CgnsGrid+index_grid-1)->Array_Name3,name); 
				}
				 ptrZ=((base->CgnsZone+(index_zone-1))->CgnsGrid+(index_grid-1))->Z; 
				if(domainflag==1)
				{
					
					for(k=0;k<2;k++)
					{
						for(j=0;j<dummy_size[1];j++)
						{
							for(i=0;i<dummy_size[0];i++)
							{
								if(k==0)
								*(ptrZ+dummy_size[0]*dummy_size[1]*k+dummy_size[0]*j+i)=0.0;
								else	
								*(ptrZ+dummy_size[0]*dummy_size[1]*k+dummy_size[0]*j+i)=1.0;
							}
						}
					}
				}
				else if(domainflag==2)
				{
					
					theta=0.1;
					grad=tan(theta/2.00);
					*(base->BaseParameters->AxisymAzimuthangle)=theta;
					for(k=0;k<2;k++)
					{
						for(j=0;j<dummy_size[1];j++)
						{
							for(i=0;i<dummy_size[0];i++)
							{
								if(*(ptrY+dummy_size[0]*j+i)<0e0)
								{
									printf("\nGRID CONTAINS NEGATIVE Y CO-ORDINATE\n");
									exit(1);
								}
								if(k==0)
								*(ptrZ+dummy_size[0]*dummy_size[1]*k+dummy_size[0]*j+i)=-1.0*grad*(*(ptrY+dummy_size[0]*j+i));
								else	
								*(ptrZ+dummy_size[0]*dummy_size[1]*k+dummy_size[0]*j+i)=grad*(*(ptrY+dummy_size[0]*j+i));
							}
						}
					}
				}
				
			} /* Loop over the Grid Coordinate Nodes ends here */ 

			

			/* ..... Reading the Number of Interface Nodes  */
			ier = cg_n1to1(index_file,index_base,index_zone,&(base->CgnsZone+(index_zone-1))->Number_of_Interfaces);
			CGNSErrorCheck(ier);
			fprintf(sum,"\t\tg) Number of Interfaces        : %d\n",(base->CgnsZone+(index_zone-1))->Number_of_Interfaces); 
			fprintf(sum,"\t\th) Summary of the Interfaces   : \n");
			/* ..... Memory Allocation for Interface */


			(base->CgnsZone+(index_zone-1))->Iface = (Interface*)malloc((base->CgnsZone+(index_zone-1))->Number_of_Interfaces*sizeof(Interface));
			if((base->CgnsZone+(index_zone-1))->Iface == NULL) {printf("\n19Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
			/* ..... Loop over the Interface nodes for reading the data of all the interfaces in a zone */



			for(index_iface=1;index_iface<=(base->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
			{

				/* ....... Memory Allocation for Point Range */
				((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range = (int*)malloc(2*3*(sizeof(int)));
				if(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range == NULL) {printf("\n20Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

				/* ....... Memory Allocation for Point Donor Range */
				((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range = (int*)malloc(2*3*(sizeof(int)));
				if(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range == NULL) {printf("\nMemory_Warning exiting from ReadCgnsFile\n");exit(1);}

				/* ....... Memory Allocation for Transformation Rule */
				((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform = (int*)malloc(3*(sizeof(int)));
				if(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform == NULL) {printf("\n22Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

				/* ....... Getting the info for the Interface */ 
				ier = cg_1to1_read(index_file,index_base,index_zone,index_iface,((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Interface_Name, 
												((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Adjoining_Zone_Name,
												dummy_pr,dummy_pdr,dummy_tr);
												  

				CGNSErrorCheck(ier);	
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range)=dummy_pr[0];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+1)=dummy_pr[1];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+2)=1;
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+3)=dummy_pr[2];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+4)=dummy_pr[3];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+5)=2;
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range)=dummy_pdr[0];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+1)=dummy_pdr[1];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+2)=1;
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+3)=dummy_pdr[2];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+4)=dummy_pdr[3];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+5)=2;
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform)=dummy_tr[0];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform+1)=dummy_tr[1];
				*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform+2)=3;
				fprintf(sum,"\t\t\t * Interface Number %d \n",index_iface);
				fprintf(sum,"\t\t\t\t   i) Name of the Interface      : %s\n",((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Interface_Name);
				fprintf(sum,"\t\t\t\t  ii) Name of the Adjoining Zone : %s\n",((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Adjoining_Zone_Name);      

				fprintf(sum,"\t\t\t\t iii) Point Range                : [%d  %d  %d ] -> [%d  %d  %d ]\n",*((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range,
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+1),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+2),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+3),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+4),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+5));
				fprintf(sum,"\t\t\t\t  iv) Point Donor Range          : [%d  %d  %d ] -> [%d  %d  %d ]\n",*((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range,
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+1),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+2),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+3),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+4),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Donor_Range+5));
				fprintf(sum,"\t\t\t\t   v) Transformation Rule        : [%d  %d  %d ]\n"               ,*((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform,
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform+1),
															*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Transform+2));

			} 


			/* ..... Reading the Boundary nodes of the zone */

			ier = cg_nbocos(index_file,index_base,index_zone,&nbc);  
			CGNSErrorCheck(ier);
			//printf("%d....",nbc);
			/* ..... Memory Allocation for Boundary */
			k=4-((base->CgnsZone+(index_zone-1))->Number_of_Interfaces);
			//printf("k=%d..........",k);
			(base->CgnsZone+(index_zone-1))->ZoneBc = (ZoneBoundary*)malloc((k+2)*sizeof(ZoneBoundary));
			if((base->CgnsZone+(index_zone-1))->ZoneBc == NULL) {printf("\n23Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
			/* ..... Loop over the Boundary Patches for getting info */

			int No_Ofbc=0;

			for(index_bc = 1;index_bc<=nbc;index_bc++) 
			{
				/* ....... Memory Allocation for Point Range */


				ier = cg_boco_info(index_file,index_base,index_zone,index_bc,Bond_Cond,&bctype,&ptsettype,&npnts,B,B,&nordt,B);
				CGNSErrorCheck(ier); 

				ier=cg_boco_read(index_file,index_base,index_zone,index_bc,dummyA,norm);
				CGNSErrorCheck(ier); 
				if(index_bc==1)
				{
					if (dummyA[0]<dummyA[2])
					{
						min_ind[0]=dummyA[0];
						max_ind[0]=dummyA[2];
					}
					else
					{
						min_ind[0]=dummyA[2];
						max_ind[0]=dummyA[0];
					}
					if(dummyA[1]<dummyA[3])
					{
						min_ind[1]=dummyA[1];
						max_ind[1]=dummyA[3];	
					}
					else
					{
						min_ind[1]=dummyA[3];
						max_ind[1]=dummyA[1];		
					}
				}
				else
				{
					if(min_ind[0]>dummyA[0])
					min_ind[0]=dummyA[0];
					if(min_ind[0]>dummyA[2])
					min_ind[0]=dummyA[2];
					if(min_ind[1]>dummyA[1])
					min_ind[1]=dummyA[1];
					if(min_ind[1]>dummyA[3])
					min_ind[1]=dummyA[3];
					
					if(max_ind[0]<dummyA[0])
					max_ind[0]=dummyA[0];
					if(max_ind[0]<dummyA[2])
					max_ind[0]=dummyA[2];
					if(max_ind[1]<dummyA[1])
					max_ind[1]=dummyA[1];
					if(max_ind[1]<dummyA[3])
					max_ind[1]=dummyA[3];
					
				}
				A[0]=dummyA[0];
				A[1]=dummyA[1];
				A[2]=1;
				A[3]=dummyA[2];
				A[4]=dummyA[3];
				A[5]=2;
				
				for(index_iface=1;index_iface<=(base->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
				{
					test[index_iface-1]=0;
					for(i=0;i<6;i++)
					{
						if( A[i]==*(((base->CgnsZone+(index_zone-1))->Iface+(index_iface-1))->Point_Range+i))
							test[index_iface-1]=1+test[index_iface-1];

						else
							test[index_iface-1]=0+test[index_iface-1];
					}
				}     
				test123=0;

				for(index_iface=1;index_iface<=(base->CgnsZone+(index_zone-1))->Number_of_Interfaces;index_iface++)
				{	if(test[index_iface-1]==6)
					{
						test123=6;
						break;
					}
					else
					{
						test123=0;	
					}
				}		 
				if(test123==0)
				{ 	

					for(i=0;i<33;i++)
					((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->BCName[i]=Bond_Cond[i];

					((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range = (int*)malloc(2*3*(sizeof(int)));
					if(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range == NULL) {printf("\n24Memory_Warning exiting from ReadCgnsFile\n");exit(1);}

					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range)=*A;
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+1)=*(A+1);
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+2)=*(A+2);
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+3)=*(A+3);
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+4)=*(A+4);
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+5)=*(A+5);
					((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Boundary_Tag=bctype;

					fprintf(sum,"\t\t\t * Boundary Patch Number %d\n",index_bc);
					fprintf(sum,"\t\t\t\t   i) Name of the Boundary Patch : %s\n",((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->BCName);
					fprintf(sum,"\t\t\t\t  ii) Type Of Boundary           : %d\n",((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Boundary_Tag);
					fprintf(sum,"\t\t\t\t iii) Point set Type             : %d\n",ptsettype);
					fprintf(sum,"\t\t\t\t  iv) Point Range                : [%d  %d  %d] -> [%d  %d  %d]\n",*((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range,A[1],A[2],A[3],A[4],A[5]);
					No_Ofbc=No_Ofbc+1;
				}
			}
				
				//Added by chandrasekhar
			//	if(strcmp(base->basename,"BASE") == 0 || strcmp(base->basename,"BASE#1") == 0)//For Grids generated by GRIDPRO and ICEMCFD
				{
					min_ind[0]=1;
					min_ind[1]=1;
					max_ind[0]=dummy_size[0];
					max_ind[1]=dummy_size[1];
				}
				//End 
				
				strcpy(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->BCName,kmin);
				strcpy(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->BCName,kmax);
				for(i=0;i<2;i++)
				{
				((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+i))->Point_Range = (int*)malloc(2*3*(sizeof(int)));
				if(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+i))->Point_Range == NULL) {printf("\n25Memory_Warning exiting from ReadCgnsFile\n");exit(1);}
				}
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range)=min_ind[0];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+1)=min_ind[1];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+2)=1;
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+3)=max_ind[0];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+4)=max_ind[1];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+5)=1;
				// KMinBc is assigned a BoundaryTag of 13 ie outflow. This KMin Bc is generated during conversion of 2d grid to 3d grid;
					((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Boundary_Tag=13;//bctype;
					
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Point_Range)=min_ind[0];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Point_Range+1)=min_ind[1];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Point_Range+2)=2;
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Point_Range+3)=max_ind[0];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Point_Range+4)=max_ind[1];
					*(((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Point_Range+5)=2;
				// KMaxBc is assigned a BoundaryTag of 13 ie outflow. This KMaxBc is generated during conversion of 2d grid to 3d grid;
					((base->CgnsZone+(index_zone-1))->ZoneBc+(No_Ofbc+1))->Boundary_Tag=13;//bctype;
					
					fprintf(sum,"\t\t\t * Boundary Patch Number %d\n",index_bc+1);
					fprintf(sum,"\t\t\t\t   i) Name of the Boundary Patch : %s\n",((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->BCName);
					fprintf(sum,"\t\t\t\t  ii) Type Of Boundary           : %d\n",((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Boundary_Tag);
					fprintf(sum,"\t\t\t\t iii) Point set Type             : %d\n",ptsettype);
					fprintf(sum,"\t\t\t\t  iv) Point Range                : [%d  %d  %d] -> [%d  %d  %d]\n",*((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range,
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+1),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+2),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+3),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+4),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc)->Point_Range+5));
																
					fprintf(sum,"\t\t\t * Boundary Patch Number %d\n",index_bc+2);
					fprintf(sum,"\t\t\t\t   i) Name of the Boundary Patch : %s\n",((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->BCName);
					fprintf(sum,"\t\t\t\t  ii) Type Of Boundary           : %d\n",((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Boundary_Tag);
					fprintf(sum,"\t\t\t\t iii) Point set Type             : %d\n",ptsettype);
					fprintf(sum,"\t\t\t\t  iv) Point Range                : [%d  %d  %d] -> [%d  %d  %d]\n",*((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Point_Range,
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Point_Range+1),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Point_Range+2),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Point_Range+3),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Point_Range+4),
																*(((base->CgnsZone+(index_zone-1))->ZoneBc+No_Ofbc+1)->Point_Range+5));											
				   (base->CgnsZone+(index_zone-1))->Number_of_BC_Patches = No_Ofbc+2;
				  fprintf(sum,"\t\ti) Number of Zone Boundaries     : %d\n",(base->CgnsZone+(index_zone-1))->Number_of_BC_Patches);
				  	
		}


return(0);
}

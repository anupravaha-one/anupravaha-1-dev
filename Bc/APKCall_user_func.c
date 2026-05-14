#include "APKCall_user_func.h"

/**
The function Call_user_func() is used for those boundary patches where user wants to use functions. 
Call_user_func() calls the desired function out of fifteen user_func_ ..(), for a perticular BC patch.

The outer switch-case is  for six Boundary patches, each case in turn contains switch 
and fifteen cases for calling user_func_...(). Thus ensuring the complit genrality out of 15 odd functions,
user can use any function for any bc patch of a particular zone. Note the this also allows uses of a perticular 
function on more than one BC patch for a particular zone. 
*/


int Call_user_func(int ind,double *dum,Zone *z,ZoneBoundary *zonebc,double time)
{
/*(funcchoice,((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct,CgnsZone,(CgnsZone->ZoneBc+BC_index),0.0);*/
  int direct_index,ni,nj,nk,i,j,k,*ptrsize;
  int count,nami,po;
  int c1,c2,c3,c4,c5;
	double *crdX,*crdY,*crdZ;
	double cen_x,cen_y,cen_z;
	//double MeanDistance,half;
	//double centroid[3];
	char name1[]="CoordinateX";
	char name2[]="CoordinateY";
	char name3[]="CoordinateZ";
	FILE *sol;    
        sol=fopen("SOLUTIONS/readu4bc.dat","r");
	
        if(sol==NULL) {
	  printf("\nError opening file SOLUTIONS/readu4bc.dat to read....probably it does not exist....Trying to create....\n");
          sol=fopen("SOLUTIONS/readu4bc.dat","w");
          if(sol==NULL) {printf("\nFailed....Exiting....\n"); exit(12);} else {printf("\nCreated.....\n");}
	}
        fclose(sol);

        sol=fopen("SOLUTIONS/readu4bc.dat","r");
        if(sol==NULL)
	{
	 printf("File opening failed in Call_user_func.c\n\n");
	 printf("Could not open file SOLUTIONS/readu4bc.dat\n");
	 printf("Exiting ....\n");
	 exit(12);
	}	
	/*This following routine separates the coordinates of all the vertices according to the corresponding 
	coordinate name and stores it in crdX,crdY,crdZ*/
		c1=strcmp(z->CgnsGrid->Array_Name1,name1);
		c2=strcmp(z->CgnsGrid->Array_Name1,name2);
		c3=strcmp(z->CgnsGrid->Array_Name1,name3);
		c4=strcmp(z->CgnsGrid->Array_Name2,name2);
		c5=strcmp(z->CgnsGrid->Array_Name2,name3);
		if(c1==0)
		{
			crdX=z->CgnsGrid->X;	
			if(c4==0)
			{
			crdY=z->CgnsGrid->Y;
			crdZ=z->CgnsGrid->Z;	
			}
			else
			{
			crdY=z->CgnsGrid->Z;
			crdZ=z->CgnsGrid->Y;	
			}
		}
		
		 if(c2==0)
		{
			crdX=z->CgnsGrid->Y;
			if(c5==0)
			{
				crdY=z->CgnsGrid->Z;
				crdZ=z->CgnsGrid->X;
			}
			else
			{
				crdY=z->CgnsGrid->X;
				crdZ=z->CgnsGrid->Z;
			}
		}
		if(c3==0)
		{
			crdX=z->CgnsGrid->Z;
			if(c4==0)
			{
				crdY=z->CgnsGrid->Y;
				crdZ=z->CgnsGrid->X;
			}
			else
			{
				crdY=z->CgnsGrid->X;
				crdZ=z->CgnsGrid->Y;
			}
		}
	
	/*Direction Index is the number which  indicates the direction of the boundary face*/
	/*Please refer direction.h for the direction values*/
	direct_index=zonebc->Direction;
	ni=*(z->Size);
	nj=*(z->Size+1);
	nk=*(z->Size+2); 
	ptrsize=z->Size;

		/*Looping over all the boundary patches according to their corresponding directions*/
	switch(direct_index)
	{
		case 1://face Imax
			i=ni-1;//<--- This is nowhere used.
			switch(ind)
			{
				case 1:
						/* The looping for cells on the face Imax starts here*/
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
							/*'cen_x''cen_y''cen_z' are the face centroids*/
							/*face_center is the subroutine which calculates the centroid of the face
							which takes the coordinates of all the vertices on that face as input and 
							coordinates of face centroid as output*/
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
							/*'user_func_1' will take these values of face centroid to the particular 
							mathematical function as prescribed by the user while entering 
							the boundary condition*/
								*(dum+(k*(nj-1))+j)=user_func_1(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
							}
						}
				
				break;
					
				case 2:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_2(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				
				case 3:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								
								*(dum+(k*(nj-1))+j)=user_func_3(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				
				case 4:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_4(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				
				case 5:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								
								*(dum+(k*(nj-1))+j)=user_func_5(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;

				case 6:
				     
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_6(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				
				case 7:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_7(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 8:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								
								*(dum+(k*(nj-1))+j)=user_func_8(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				
				case 9:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_9(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 10:
				     		for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_10(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 11:
				     	
						for(k=0;k<nk-1;k++)
						{	
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_11(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 12:
				
					
						for(k=0;k<nk-1;k++)
						{	
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_12(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 13:
				
					
						for(k=0;k<nk-1;k++)
						{	
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_13(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 14:
						for(k=0;k<nk-1;k++)
						{	
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_14(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 15:
						for(k=0;k<nk-1;k++)
						{	
							for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_15(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 16:
						for(k=0;k<nk-1;k++)
						{	
							for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)+(ni-1)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_16(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
	                        case 20:
                                         	nami = (nk-1)*(nj-1);
                                         	for (po=0;po<nami;po++)
						  {                                          	
								
							      fscanf(sol,"%lf",(dum+po));
							      //    printf("%lf\n",*(dum+po));
                 					
						  }
					
				break;
			}
		break;
				
		case 2://  face Imin
			i=0;
			switch(ind)
			{
				case 1: 
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_1(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				
				case 2:
				        
						for(k=0;k<nk-1;k++)
						{	for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_2(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 3:
				        
						for(k=0;k<nk-1;k++)
						{	for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_3(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				
				case 4:
				        for(k=0;k<nk-1;k++)
						{	for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_4(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 5:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								
								*(dum+(k*(nj-1))+j)=user_func_5(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 6:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_6(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 7:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_7(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 8:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								
								*(dum+(k*(nj-1))+j)=user_func_8(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 9:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_9(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 10:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_10(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 11:
				     	
					
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_11(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 12:
				    	
					
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_12(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 13:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_13(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 14:
						for(k=0;k<nk-1;k++)
						{
							for(j=0;j<nj-1;j++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_14(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 15:
						for(k=0;k<nk-1;k++)
						{	for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_15(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				
				case 16:
						for(k=0;k<nk-1;k++)
						{	for(j=0;j<nj-1;j++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(j)),direct_index,ni,nj,nk);
								*(dum+(k*(nj-1))+j)=user_func_16(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
                                case 20:       
	                                       nami = (nk-1)*(nj-1);
                                               	for (po=0;po<nami;po++)
						  {                                          	
								
							      fscanf(sol,"%lf",(dum+po));
							      //   printf("%lf\n",*(dum+po));
                 					
						  }
				
				break;
				
				
			}	
		
		break;
		case 4:/*  face Jmin */
			j=0;
			
			switch(ind)/*Switch according to function index entered by the user*/
			{
				case 1:	
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_1(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;				
				case 2:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_2(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 3:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_3(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 4:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_4(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 5:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								
								*(dum+(k*(ni-1))+i)=user_func_5(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 6:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_6(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 7:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_7(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 8:
						for(k=0;k<nk-1;k++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								
								*(dum+(k*(ni-1))+i)=user_func_8(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 9:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_9(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 10:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{	
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_10(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 11:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_11(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 12:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_12(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 13:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_13(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 14:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_14(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 15:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_15(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 16:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*0+(i)),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_16(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
						
				break;
                                case 20:
						nami = (nk-1)*(ni-1);
                                          	for (po=0;po<nami;po++)
						  {                                          	
								
							      fscanf(sol,"%lf",(dum+po));
							      //	      printf("%lf\n",*(dum+po));
                 					
						  }
				
					
				break;
			}
		break;
			
		case 3:/*  face jmax */
		 	j=nj-1;
			switch(ind)/*Switch according to function index entered by the user*/
			{
				case 1:	
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								 cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_1(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 2:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_2(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 3:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_3(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 4:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_4(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 5:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								
								*(dum+(k*(ni-1))+i)=user_func_5(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
							}
						}
				
				break;
				case 6:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_6(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 7:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{		
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_7(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 8:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								
								*(dum+(k*(ni-1))+i)=user_func_8(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 9:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_9(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
							
                 					}
						}
				break;
				case 10:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_10(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 11:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_11(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;		
				case 12:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_12(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
							}
						}
				break;
				case 13:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_13(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
							
                 					}
						}
				
				break;
				case 14:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_14(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				
				break;
				case 15:
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_15(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
								
				case 16:	
						for(k=0;k<nk-1;k++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(k)+ni*(nj-1)+i),direct_index,ni,nj,nk);
								*(dum+(k*(ni-1))+i)=user_func_16(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
                                case 20:
                                              	nami = (nk-1)*(ni-1);
				         	for (po=0;po<nami;po++)
						  {                                          	
								
							      fscanf(sol,"%lf",(dum+po));
							      //	      printf("%lf\n",*(dum+po));
                 					
						  }
				
					
				break;
			}
		break;
		 	
		case 6:/* The looping for cells on the face Kmin starts here*/
		 	k=0;
			switch(ind)/*Switch according to function index entered by the user*/		
			{
				case 1:
					count=0;
					for(j=0;j<nj-1;j++)
					{
						for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_1(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}

					}
				break;	
				case 2:
					for(j=0;j<nj-1;j++)
					{	for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_2(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}
					}
				
				break;
				case 3:
					for(j=0;j<nj-1;j++)
					{	for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_3(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}
					}
				
				break;
				case 4:
					for(j=0;j<nj-1;j++)
					{	for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_4(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}
					}
				break;	
				case 5:
					for(j=0;j<nj-1;j++)
					{	for(i=0;i<ni-1;i++)
						{		
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_5(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}
					}
				
				break;
				case 6:
					for(j=0;j<nj-1;j++)
					{	for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_6(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);

						}
					}
				
				break;
				case 7:
					for(j=0;j<nj-1;j++)
					{	for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_7(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);

						}
					}
				
				break;
				case 8:
					for(j=0;j<nj-1;j++)
					{
						for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);

							*(dum+((j*(ni-1))+i))=user_func_8(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}
					}
				
				break;
				case 9:
					for(j=0;j<nj-1;j++)
					{
						for(i=0;i<ni-1;i++)
						{
							cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
							*(dum+((j*(ni-1))+i))=user_func_9(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
						}
					}
				break;
				case 10:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_10(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
								
                 					}
						}
				break;
				case 11:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_11(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
							}
						}
				break;
				case 12:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_12(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 13:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_13(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 14:
						for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_14(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 15:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_15(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 16:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(0)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_16(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
                                case 20:
						nami = (ni-1)*(nj-1);	
	                                         for (po=0;po<nami;po++)
						  {                                          	
								
							      fscanf(sol,"%lf",(dum+po));
							      //	      printf("%lf\n",*(dum+po));
                 					
						  }
				
					
				break;
			}	
		break;
			
		case 5:/*Face Kmin*/
			k=nk-1;
			switch(ind)/*Switch according to function index entered by the user*/
			{
				case 1:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_1(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 2:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_2(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 3:
				        	for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_3(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 4:
				        	for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_4(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 5:
						for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_5(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 6:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_6(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 7:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_7(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 8:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								
								*(dum+((j*(ni-1))+i))=user_func_8(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 9:
						for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_9(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 10:
						for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_10(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 11:
						for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_11(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 12:
						for(j=0;j<nj-1;j++)
						{
							for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_12(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 13:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_13(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 14:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_14(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				
				break;
				case 15:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_15(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
				case 16:
						for(j=0;j<nj-1;j++)
						{	for(i=0;i<ni-1;i++)
							{
								cen_x=face_center((crdX+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_y=face_center((crdY+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								cen_z=face_center((crdZ+(ni*nj)*(nk-1)+ni*(j)+i),direct_index,ni,nj,nk);
								*(dum+((j*(ni-1))+i))=user_func_16(cen_x,cen_y,cen_z,time,z,i,j,k,zonebc);
                 					}
						}
				break;
                                case 20:
				              	nami = (ni-1)*(nj-1);
                                                
                                                  for (po=0;po<nami;po++)
						  {                                          	
								
							      fscanf(sol,"%lf",(dum+po));
							      //	      printf("%lf\n",*(dum+po));
                 					
						  }
				
				
				break;
			}
		break;
	
		
	}
fclose(sol);
return (0);
}


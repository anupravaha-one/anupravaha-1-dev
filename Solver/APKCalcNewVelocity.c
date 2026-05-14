#include "APKCalcNewVelocity.h"




double CalcNewVelocity(Base *base)
{

  
      	double *volume,vol,*poten_source;


	int NI,NJ,NK,k,j,i,indexP,*ptrSize;
	int numZones,nZones;

        int c1,c2,c3,c4,c5,face;
	double *crdX,*crdY,*crdZ;
	double F[6],cx[6],cy[6],cz[6],centroid[3];



	char name1[]="CoordinateX";
	char name2[]="CoordinateY";
	char name3[]="CoordinateZ";




	Zone *zone;



		nZones=(base->Number_of_zones);



  
     
      
           
       
    	for(numZones=0;numZones<nZones;numZones++)
         	{
	        	zone = (base->CgnsZone+numZones);
		
	        	ptrSize= zone->Size;

			c1=strcmp(zone->CgnsGrid->Array_Name1,name1);
			c2=strcmp(zone->CgnsGrid->Array_Name1,name2);
			c3=strcmp(zone->CgnsGrid->Array_Name1,name3);
			c4=strcmp(zone->CgnsGrid->Array_Name2,name2);
			c5=strcmp(zone->CgnsGrid->Array_Name2,name3);
		
		if(c1==0)
		{
			crdX=zone->CgnsGrid->X;	
			if(c4==0)
			{
			crdY=zone->CgnsGrid->Y;
			crdZ=zone->CgnsGrid->Z;	
			}
			else
			{
			crdY=zone->CgnsGrid->Z;
			crdZ=zone->CgnsGrid->Y;	
			}
		}
		
		 if(c2==0)
		{
			crdX=zone->CgnsGrid->Y;
			if(c5==0)
			{
				crdY=zone->CgnsGrid->Z;
				crdZ=zone->CgnsGrid->X;
			}
			else
			{
				crdY=zone->CgnsGrid->X;
				crdZ=zone->CgnsGrid->Z;
			}
		}
		if(c3==0)
		{
			crdX=zone->CgnsGrid->Z;
			if(c4==0)
			{
				crdY=zone->CgnsGrid->Y;
				crdZ=zone->CgnsGrid->X;
			}
			else
			{
				crdY=zone->CgnsGrid->X;
				crdZ=zone->CgnsGrid->Y;
			}
		}


	        	NI=*(ptrSize+0);
	           	NJ=*(ptrSize+1);
	           	NK=*(ptrSize+2);

	         
    
                       

	                        volume=MakeDomainArray(zone,5);
                          
	
			
                  
				
	           	for(k=0;k<NK-1;k++)
	          	for(j=0;j<NJ-1;j++)
	          	for(i=0;i<NI-1;i++)
	        	{
		        	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		           	
                        	vol = *(volume+indexP);
			
				//	printf("\n reizoku oou");
                               
					poten_source = zone->Flowsol->Massflux+indexP*6;

			face = 0;  // west face
				F[face]  = ( *(poten_source));                             
			        cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			     
			face = 1;  // east face
				F[face]  = ( *(poten_source+1));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 2;  // south face
			     	F[face]  = ( *(poten_source+2));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
		 	        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
		       face = 3;  // north face
			 	F[face]  = ( *(poten_source+3));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			       
			face = 4;  // bottom face
			 	F[face]  = ( *(poten_source+4));
		                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 5;  // top face
			 	F[face]  = ( *(poten_source+5));
			        cx[face] = face_center((crdX+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
 

			    	
		           
		     
				*((zone->Flowsol)->utemp+indexP) =((F[0]*cx[0])+(F[1]*cx[1])+(F[2]*cx[2])+(F[3]*cx[3])+(F[4]*cx[4])+(F[5]*cx[5]))/vol;
			
                
				*((zone->Flowsol)->vtemp+indexP) = ((F[0]*cy[0])+(F[1]*cy[1])+(F[2]*cy[2])+(F[3]*cy[3])+(F[4]*cy[4])+(F[5]*cy[5]))/vol;

				 if(base->BaseParameters->flag2d==0)
				   {
					*((zone->Flowsol)->wtemp+indexP) =  ((F[0]*cz[0])+(F[1]*cz[1])+(F[2]*cz[2])+(F[3]*cz[3])+(F[4]*cz[4])+(F[5]*cz[5]))/vol;
				   }
                               
				

			
			}
		 	
				}
         
     

   
  
       return 0;
	}

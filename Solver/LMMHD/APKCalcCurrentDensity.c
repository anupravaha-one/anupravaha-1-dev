#include "APKCalcCurrentDensity.h"
/**
This function is used to calculate the current density at each cell centroid using Ohm's law
*/
double CalcCurrentDensity(Base *base)
{

        double dBx_dY,dBx_dZ,dBy_dX,dBy_dZ,dBz_dX,dBz_dY;
      	double *Bx,*By,*Bz,*volume,*SurfArea,vol,*phi,*poten_source,*beta,*diffusionCoeffdomain;
	double *slab_Bx[6],*slab_By[6],*slab_Bz[6],*slabVolume[6],*gradPhiSlab[6];

	int NI,NJ,NK,k,j,i,indexP,index[3],*ptrSize;
	int numZones,nZones,flagelectcond,flagpermeability;
	double B_x, B_y, B_z, velPu, velPv, velPw;
	double  sigma, Dphi[3],magperm,diffdia;
	double CoeffN[6],potDiff[6];
        int c1,c2,c3,c4,c5,face;
	double *crdX,*crdY,*crdZ;
	double F[6],cx[6],cy[6],cz[6],centroid[3];



	char name1[]="CoordinateX";
	char name2[]="CoordinateY";
	char name3[]="CoordinateZ";



	//double J_x,J_y,J_z;
	Zone *zone;
	
	FILE *bp;
	FILE *mp;
	FILE *kp;
	FILE *hp;
	bp=fopen("SOLUTIONS/uxbds.dat","w");
	mp=fopen("SOLUTIONS/minusgradp.dat","w");
        kp=fopen("SOLUTIONS/mgradppuxbds.dat","w");
        hp=fopen("SOLUTIONS/currentdensity.dat","w");
	nZones=(base->Number_of_zones);
	/*J_x=(zone->Flowsol)->ElecurrentdensityX;
	J_y=(zone->Flowsol)->ElecurrentdensityY;
	J_z=(zone->Flowsol)->ElecurrentdensityZ;*/



    if(base->BaseParameters->flag2d==1)
     {
       if(base->flaginduced==1)
           {
            	for(numZones=0;numZones < nZones;numZones++)
		{
			zone=(base->CgnsZone+numZones);
		
			{
				ptrSize= zone->Size;
	
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);

                   		flagpermeability=zone->Flowsol->flagpermeability;
	                     
				volume=MakeDomainArray(zone,5);
	
                               
		           	Bx=MakeDomainArray(zone,56);
				By=MakeDomainArray(zone,57);
			
				MakeSlabPtrArray(zone,slab_Bx,56);
				MakeSlabPtrArray(zone,slab_By,57);
	
				SurfArea=zone->ZGeom->SurfaceArea;
	
				MakeSlabPtrArray(zone,slabVolume,5);
	
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index[0]=k,index[1]=j,index[2]=i;
	
					indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;
	
                                           magperm = *(zone->Flowsol->Permeability+indexP*flagpermeability);
					dBx_dY=CalcPressureGrad(Bx,slab_Bx,volume,slabVolume,SurfArea,index,ptrSize,2);
					dBy_dX=CalcPressureGrad(By,slab_By,volume,slabVolume,SurfArea,index,ptrSize,1);
					
					vol = *(volume+indexP);
					
					*((zone->Flowsol)->ElecurrentdensityX+indexP)=0.0;
					*((zone->Flowsol)->ElecurrentdensityY+indexP)=0.0;
					*((zone->Flowsol)->ElecurrentdensityZ+indexP)=(-1.0/(magperm*vol))*(dBy_dX-dBx_dY);
				}
			}

		}
           }
       else{
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

	         	flagelectcond=zone->Flowsol->flageleccond;
    
                        beta = zone->ZGeom->Beta;

		           	phi=MakeDomainArray(zone,17);
                         	MakeSlabPtrArray(zone,gradPhiSlab,17);
	                        volume=MakeDomainArray(zone,5);
                          	SurfArea=zone->ZGeom->SurfaceArea;
	
				MakeSlabPtrArray(zone,slabVolume,5);
                  	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,17);

	           	for(k=0;k<NK-1;k++)
	          	for(j=0;j<NJ-1;j++)
	          	for(i=0;i<NI-1;i++)
	        	{
		        	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		           	index[0]=k,index[1]=j,index[2]=i;
                         	vol = *(volume+indexP);
                                diffdia = Diffusion2dWithConstDiffCoeff_opti(index,phi,volume,ptrSize,gradPhiSlab,slabVolume,diffusionCoeffdomain,beta,potDiff,CoeffN);
                                FindCellCentre(ptrSize,(crdX+(NI*NJ)*(k)+NI*(j)+(i)),(crdY+(NI*NJ)*(k)+NI*(j)+(i)),(crdZ+(NI*NJ)*(k)+NI*(j)+(i)),centroid);
				poten_source = zone->Flowsol->PotentialMHDpoissonRHS+indexP*6;

			face = 0;  // west face
				F[face]  = (-potDiff[face]+ *(poten_source));                               
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
                                    
			face = 1;  // east face
				F[face]  = (-potDiff[face]+ *(poten_source+1));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 2;  // south face
			     	F[face]  = (-potDiff[face]+ *(poten_source+2));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
		 	        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 3;  // north face
			 	F[face]  = (-potDiff[face]+ *(poten_source+3));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			       
			face = 4;  // bottom face
			 	F[face]  = (-potDiff[face]+ *(poten_source+4));
		                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 5;  // top face
			 	F[face]  = (-potDiff[face]+ *(poten_source+5));
			        cx[face] = face_center((crdX+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
 

			          	
		           	sigma = *(zone->PrevFlowsol->Eleconductivity+indexP*flagelectcond);
		     
				*((zone->Flowsol)->ElecurrentdensityX+indexP) = sigma*((F[0]*cx[0])+(F[1]*cx[1])+(F[2]*cx[2])+(F[3]*cx[3])+(F[4]*cx[4])+(F[5]*cx[5]))/vol;
			
                
				*((zone->Flowsol)->ElecurrentdensityY+indexP) = sigma*((F[0]*cy[0])+(F[1]*cy[1])+(F[2]*cy[2])+(F[3]*cy[3])+(F[4]*cy[4])+(F[5]*cy[5]))/vol;
				
                   
				*((zone->Flowsol)->ElecurrentdensityZ+indexP) =  sigma*((F[0]*cz[0])+(F[1]*cz[1])+(F[2]*cz[2])+(F[3]*cz[3])+(F[4]*cz[4])+(F[5]*cz[5]))/vol;

                                fprintf(mp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",-potDiff[0],-potDiff[1],-potDiff[2],-potDiff[3],-potDiff[4],-potDiff[5]);
 
				fprintf(bp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",	*(poten_source),*(poten_source+1),*(poten_source+2),*(poten_source+3),*(poten_source+4),*(poten_source+5));
				fprintf(kp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",-potDiff[0]+(*(poten_source)),-potDiff[1]+(*(poten_source+1)),-potDiff[2]+(*(poten_source+2)),-potDiff[3]+(*(poten_source+3)),-potDiff[4]+(*(poten_source+4)),-potDiff[5]+(*(poten_source+5)));
		
				fprintf(hp,"%lf \t %lf \t %lf \n",(potDiff[0]-(*(poten_source))-potDiff[1]+(*(poten_source+1)))*0.5,(potDiff[2]-(*(poten_source+2))-potDiff[3]+(*(poten_source+3)))*0.5,(potDiff[4]-(*(poten_source+4))-potDiff[5]+(*(poten_source+5)))*0.5);



	//*((zone->Flowsol)->Elecurrentdensity+indexP)= sqrt(J_x*J_x+J_y*J_y+J_z*J_z);
	            	}
		 	
				}
       }  
     }

    else{
         
          if(base->flaginduced==1)
	  {
          for(numZones=0;numZones < nZones;numZones++)
		{
			zone=(base->CgnsZone+numZones);
		
			{
				ptrSize= zone->Size;
	
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);

                   		flagpermeability=zone->Flowsol->flagpermeability;
	                       
				volume=MakeDomainArray(zone,5);
	
                                Bx=MakeDomainArray(zone,56);
				By=MakeDomainArray(zone,57);
				Bz=MakeDomainArray(zone,58);

				MakeSlabPtrArray(zone,slab_Bx,56);
				MakeSlabPtrArray(zone,slab_By,57);
	                	MakeSlabPtrArray(zone,slab_Bz,58);

				SurfArea=zone->ZGeom->SurfaceArea;
	
				MakeSlabPtrArray(zone,slabVolume,5);
	
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index[0]=k,index[1]=j,index[2]=i;
	
					indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;
                                      	 magperm = *(zone->Flowsol->Permeability+indexP*flagpermeability);
					dBx_dY=CalcPressureGrad(Bx,slab_Bx,volume,slabVolume,SurfArea,index,ptrSize,2);
					dBx_dZ=CalcPressureGrad(Bx,slab_Bx,volume,slabVolume,SurfArea,index,ptrSize,3);
					dBy_dX=CalcPressureGrad(By,slab_By,volume,slabVolume,SurfArea,index,ptrSize,1);
					dBy_dZ=CalcPressureGrad(By,slab_By,volume,slabVolume,SurfArea,index,ptrSize,3);
					dBz_dX=CalcPressureGrad(Bz,slab_Bz,volume,slabVolume,SurfArea,index,ptrSize,1);
					dBz_dY=CalcPressureGrad(Bz,slab_Bz,volume,slabVolume,SurfArea,index,ptrSize,2);
	
					
					vol = *(volume+indexP);
					
					*((zone->Flowsol)->ElecurrentdensityX+indexP)=(-1.0/(magperm*vol))*(dBz_dY-dBy_dZ);;
					*((zone->Flowsol)->ElecurrentdensityY+indexP)=(-1.0/(magperm*vol))*(dBx_dZ-dBz_dX);;
					*((zone->Flowsol)->ElecurrentdensityZ+indexP)=(-1.0/(magperm*vol))*(dBy_dX-dBx_dY);
				}
			}

		}
	  }
          else
          {
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
  
                              beta = zone->ZGeom->Beta;

		           	phi=MakeDomainArray(zone,17);
                         	MakeSlabPtrArray(zone,gradPhiSlab,17);
	                        volume=MakeDomainArray(zone,5);
                          	SurfArea=zone->ZGeom->SurfaceArea;
	
				MakeSlabPtrArray(zone,slabVolume,5);
                            	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,17);

	        	flagelectcond=zone->Flowsol->flageleccond;

	            	for(k=0;k<NK-1;k++)
	         	for(j=0;j<NJ-1;j++)
	        	for(i=0;i<NI-1;i++)
	        	{
		         	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		          	index[0]=k,index[1]=j,index[2]=i;
                        
                            	vol = *(volume+indexP);
                               diffdia = Diffusion3dWithConstDiffCoeff_opti(index,phi,volume,ptrSize,gradPhiSlab,slabVolume,diffusionCoeffdomain,beta,potDiff);
                               FindCellCentre(ptrSize,(crdX+(NI*NJ)*(k)+NI*(j)+(i)),(crdY+(NI*NJ)*(k)+NI*(j)+(i)),(crdZ+(NI*NJ)*(k)+NI*(j)+(i)),centroid);
                        	poten_source = zone->Flowsol->PotentialMHDpoissonRHS+indexP*6;
                  	face = 0;  // west face
				F[face]  = (-potDiff[face]+ *(poten_source));                               
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),2,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
                                    
			face = 1;  // east face
				F[face]  = (-potDiff[face]+ *(poten_source+1));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i+1)),2,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 2;  // south face
			     	F[face]  = (-potDiff[face]+ *(poten_source+2));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
		 	        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),4,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 3;  // north face
			 	F[face]  = (-potDiff[face]+ *(poten_source+3));
                                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j+1)+(i)),4,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			       
			face = 4;  // bottom face
			 	F[face]  = (-potDiff[face]+ *(poten_source+4));
		                cx[face] = face_center((crdX+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k)+NI*(j)+(i)),6,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
			face = 5;  // top face
			 	F[face]  = (-potDiff[face]+ *(poten_source+5));
			        cx[face] = face_center((crdX+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
			        cx[face] = cx[face]-centroid[0];
			        cy[face] = face_center((crdY+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
                                cy[face] = cy[face]-centroid[1];
			        cz[face] = face_center((crdZ+(NI*NJ)*(k+1)+NI*(j)+(i)),6,NI,NJ,NK);
                                cz[face] = cz[face]-centroid[2];
 

                               
		         	sigma = *(zone->PrevFlowsol->Eleconductivity+indexP*flagelectcond);

                          	*((zone->Flowsol)->ElecurrentdensityX+indexP) = sigma*((F[0]*cx[0])+(F[1]*cx[1])+(F[2]*cx[2])+(F[3]*cx[3])+(F[4]*cx[4])+(F[5]*cx[5]))/vol;
			
                
				*((zone->Flowsol)->ElecurrentdensityY+indexP) = sigma*((F[0]*cy[0])+(F[1]*cy[1])+(F[2]*cy[2])+(F[3]*cy[3])+(F[4]*cy[4])+(F[5]*cy[5]))/vol;
				
                   
				*((zone->Flowsol)->ElecurrentdensityZ+indexP) =  sigma*((F[0]*cz[0])+(F[1]*cz[1])+(F[2]*cz[2])+(F[3]*cz[3])+(F[4]*cz[4])+(F[5]*cz[5]))/vol;


                                fprintf(mp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",-potDiff[0],-potDiff[1],-potDiff[2],-potDiff[3],-potDiff[4],-potDiff[5]);
 
				fprintf(bp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",	*(poten_source),*(poten_source+1),*(poten_source+2),*(poten_source+3),*(poten_source+4),*(poten_source+5));
				fprintf(kp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",-potDiff[0]+(*(poten_source)),-potDiff[1]+(*(poten_source+1)),-potDiff[2]+(*(poten_source+2)),-potDiff[3]+(*(poten_source+3)),-potDiff[4]+(*(poten_source+4)),-potDiff[5]+(*(poten_source+5)));
		
				fprintf(hp,"%lf \t %lf \t %lf \n",(potDiff[0]-(*(poten_source))-potDiff[1]+(*(poten_source+1)))*0.5,(potDiff[2]-(*(poten_source+2))-potDiff[3]+(*(poten_source+3)))*0.5,(potDiff[4]-(*(poten_source+4))-potDiff[5]+(*(poten_source+5)))*0.5);







		              
			//*((zone->Flowsol)->Elecurrentdensity+indexP)= sqrt(J_x*J_x+J_y*J_y+J_z*J_z);
	        	}
			
        	}
	  }
        }
    fclose(bp);
     fclose(mp);
fclose(kp);
fclose(hp);
	return 0;	
	}

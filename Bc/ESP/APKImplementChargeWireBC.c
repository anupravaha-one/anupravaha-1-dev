#include "APKImplementChargeWireBC.h"
/**
This function ImplementChargeWireBC() is used for implementing/correcting the charge density on wire
after specified number of timesteps as specified in the input file
/bin/Settings/ESPInput/WireChargeUpdateData.dat, this '.dat' file
contains information about the relaxation coefficients Omega and Gamma, Number
of timesteps after which charge density has to be updated, Radius of wire (meters),
coefficients of Peek's formula (A and B) and amount of tolerance between Peek's
and numerically calculated electric fields respectively.

The charge density on wire will be updated after every 'N' timesteps and only once
for each zone in that particular timestep. So inorder to avoid multiple correction
of charge density in each G-S iterations of a particular timestep, the variable
"ZoneWireTag" is defined., if this is set to "0", then it implies that, we have to
update charge density in the present timestep, once charge density is updated for
this zone in this iteration, it is immediately set to "1" so that it is not
repeated again.

Remember, however that the most general case of wire direction is SOUTH.
*/

int ImplementChargeWireBC(Base* base ,int direction,double *flowsol,int vartag,Zone *z,int size,Dataset *dset,double *dummy)
{
	int ni,nj,nk,i,j,k,*ptrSize,indexP,index[3],Out;
	double p,q,r,Area,*vols,*beta,*surf;
	double *slabVariable[6],*slabVolume[6],Fdiff[6];
	double Ewire,Epeek,*variable;
	double Tolerance, RhoStar, RhoWire, Gamma = 0.0, Omega = 0.0, errorT = 0.0, error = 0.0;
	double Radius = 1.0e-4, A=32.3e05, B=0.846e05, diffdia, CoeffN[6];
	FILE *sav,*Tstep;
	int stp, stp1;
	
	Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","r");
	fscanf(Tstep,"%d",&stp);
	fclose(Tstep);
	
	sav=fopen("Settings/ESPInput/WireChargeUpdateData.dat","r");
	fscanf(sav,"%lf",&Omega);	  //Omega value for updating wire charge
	fscanf(sav,"%lf",&Gamma);	  //Gamma value for updating wire charge
	fscanf(sav,"%d",&stp1); 	  //Number of timesteps after which updation has to take place
	fscanf(sav,"%lf",&Radius);    //Radius of wire
	fscanf(sav,"%lf",&A);         //First coefficient of Peek's formula
	fscanf(sav,"%lf",&B);         //Second coefficient of Peek's formula
	fscanf(sav,"%lf",&Tolerance); //Amount of tolerance between Epeek and Ewire
	fclose(sav);
	
	//printf("\nTimestep = %d\tOmega = %lf\tGamma = %lf\tCounter = %d\tRad = %e\tA = %lf\n", stp,Omega,Gamma,stp1,Radius,A);
	
	stp = (stp%stp1); 
	if ((stp==0) && (z->ZoneWireTag == 0))
	{
		z->ZoneWireTag = 1;
		//conv=fopen("Settings/ESPInput/IonChargeConvergency.dat","w");
		ni=*(z->Size+0);
		nj=*(z->Size+1);
		nk=*(z->Size+2); 
		
		variable=MakeDomainArray(z,8); //Electric Potential
		MakeSlabPtrArray(z,slabVariable,8); //Electric Potential		
		MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/
		
		Epeek = (A+B*sqrt(1.0/Radius))*1.0e-3;//This is in Kvolts per meter
		
		vols= z->ZGeom->Volume;
		beta= z->ZGeom->Beta;
		surf= z->ZGeom->SurfaceArea;
		ptrSize=z->Size;

		switch(direction)
		{ 
		
			case 1://east
				errorT = 0.0;error = 0.0;Out = 0;Ewire = 0.0;
				for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						i=ni-2;//last but one cell
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p=*(surf+(indexP)*6*3+3*1+0);
						q=*(surf+(indexP)*6*3+3*1+1);
						r=*(surf+(indexP)*6*3+3*1+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = fabs(Fdiff[1])/Area;
						errorT = errorT + ((Epeek - Ewire)*(Epeek - Ewire)/(Ewire*Ewire));
					}
				
				errorT = sqrt(errorT/((nj-1)*(nk-1)));
				
				if(errorT>Tolerance)
				{
					printf ("\nPercentage Error = %e\nI am updating Ion Charge..........\n",errorT);
					for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							i=ni-2;//last but one cell
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							index[0]=k,index[1]=j,index[2]=i;
							diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
							p=*(surf+(indexP)*6*3+3*1+0);
							q=*(surf+(indexP)*6*3+3*1+1);
							r=*(surf+(indexP)*6*3+3*1+2);
							Area = sqrt(p*p+q*q+r*r);
							Ewire = fabs(Fdiff[1])/Area;
							RhoWire = *(dummy+((nj-1)*k+j));
							RhoStar = RhoWire*pow((Ewire/Epeek),Gamma);
							*(dummy+((nj-1)*k+j))=(1-Omega)*RhoWire+Omega*RhoStar;
							*(dset->funct+((nj-1)*k+j)) = *(dummy+((nj-1)*k+j));

						}
						base->ChargeConv = 1;
					//fprintf(conv,"%d",1);//1 - Implies Electric field not converged
				}//end if
				else
				{
					//printf("\nI have not updated Ion Charge\n");
					printf("\nNo Error in this Iteration\n");
					base->ChargeConv = 0;
					//fprintf(conv,"%d",0);//0 - Implies Electric field converged on the wire
				}

			break;
				
			case 2://west
				errorT = 0.0;error = 0.0;Out = 0;Ewire = 0.0;
				for(k=0;k<nk-1;k++)
					for(j=0;j<nj-1;j++)
					{
						i=0;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = fabs(Fdiff[0])/Area;
						errorT = errorT + ((Epeek - Ewire)*(Epeek - Ewire)/(Ewire*Ewire));
					}
					
				errorT = sqrt(errorT/((nj-1)*(nk-1)));	
				
				if(errorT>Tolerance)
				{
					printf ("\nPercentage Error = %e\nI am updating Ion Charge..........\n",errorT);
					for(k=0;k<nk-1;k++)
					for(j=0;j<nj-1;j++)
					{
						i=0;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Ewire = fabs(Fdiff[0])/Area;
						RhoWire = *(dummy+((nj-1)*k+j));
						RhoStar = RhoWire*pow((Ewire/Epeek),Gamma);
						*(dummy+((nj-1)*k+j))=(1-Omega)*RhoWire+Omega*RhoStar;
						*(dset->funct+((nj-1)*k+j)) = *(dummy+((nj-1)*k+j));
					}
					base->ChargeConv = 1;
					//fprintf(conv,"%d",1);//1 - Implies Electric field not converged
				}//end if
				else
				{
					//printf("\nI have not updated Ion Charge\n");
					printf("\nNo Error in this Iteration\n");
					base->ChargeConv = 0;
					//fprintf(conv,"%d",0);//0 - Implies Electric field converged on the wire
				}
			break;

			case 3://north
				errorT = 0.0;error = 0.0;Out = 0;Ewire = 0.0;
				for(k=0;k<nk-1;k++) 
				{
					for(i=0;i<ni-1;i++)
					{
						j = nj-2;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p = *(surf+(indexP)*6*3+3*3+0);
						q = *(surf+(indexP)*6*3+3*3+1);
						r = *(surf+(indexP)*6*3+3*3+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = fabs(Fdiff[3])/Area;
						errorT = errorT + ((Epeek - Ewire)*(Epeek - Ewire)/(Epeek*Epeek));
					}
				}

				errorT = sqrt(errorT/((ni-1)*(nk-1)));
				
				if(errorT>Tolerance)
				{
					printf ("\nPercentage Error = %e\nI am updating Ion Charge..........\n",errorT);
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						j = nj-2;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p = *(surf+(indexP)*6*3+3*3+0);
						q = *(surf+(indexP)*6*3+3*3+1);
						r = *(surf+(indexP)*6*3+3*3+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = fabs(Fdiff[3])/Area;
						RhoWire = *(dummy+k*(ni-1)+i);
						RhoStar = RhoWire*pow((Ewire/Epeek),Gamma);
						*(dummy+k*(ni-1)+i)=(1-Omega)*RhoWire+Omega*RhoStar;
						*(dset->funct+(k*(ni-1)+i)) = *(dummy+(k*(ni-1)+i));
						
					}
					base->ChargeConv = 1;base->ChargeConv = 1;
					//fprintf(conv,"%d",1);//1 - Implies Electric field not converged
				}//end if
				else
				{
					//printf("\nI have not updated Ion Charge\n");
					printf("\nNo Error in this Iteration\n");
					base->ChargeConv = 0;
					//fprintf(conv,"%d",0);//0 - Implies Electric field converged on the wire
				}
			break;

			//Most general Direction is case 4/South face
			case 4://South
			errorT = 0.0;error = 0.0;Out = 0;Ewire = 0.0;
			for (k=0;k<nk-1;k++) 
			{
				for (i=0;i<ni-1;i++) 
				{
					j=0;
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					index[0]=k,index[1]=j,index[2]=i;
					diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
					p=*(surf+(indexP)*6*3+3*2+0);
					q=*(surf+(indexP)*6*3+3*2+1);
					r=*(surf+(indexP)*6*3+3*2+2);
					Area = sqrt(p*p+q*q+r*r);
					Ewire = 1.0*Fdiff[2]/Area;
					errorT = errorT + ((Epeek - Ewire)*(Epeek - Ewire)/(Epeek*Epeek));
				}
			}
			
			errorT = sqrt(errorT/((ni-1)*(nk-1)));
			
			if(errorT>Tolerance)
			{
				printf ("\nPercentage Error = %e\nI am updating Ion Charge..........\n",errorT);
				for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						j=0;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p=*(surf+(indexP)*6*3+3*2+0);
						q=*(surf+(indexP)*6*3+3*2+1);
						r=*(surf+(indexP)*6*3+3*2+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = 1.0*Fdiff[2]/Area;
						RhoWire = *(dset->funct+(k*(ni-1)+i));
						RhoStar = RhoWire*pow((Ewire/Epeek),Gamma);
						*(dummy+k*(ni-1)+i)=(1-Omega)*RhoWire+Omega*RhoStar;
						*(dset->funct+(k*(ni-1)+i)) = *(dummy+(k*(ni-1)+i));

					}
					base->ChargeConv = 1;
					//fprintf(conv,"%d",1);//1 - Implies Electric field not converged
			}//end if
			else
			{
				//printf("\nI have not updated Ion Charge\n");
				printf("\nNo Error between Peek's and Calculated electric fields in this Iteration\n");
				base->ChargeConv = 0;
				//fprintf(conv,"%d",0);//0 - Implies Electric field converged on the wire
			}
			
			break;
			
			case 5://top
			errorT = 0.0;error = 0.0;Out = 0;Ewire = 0.0;
			for(j=0;j<nj-1;j++) 
			{
				for(i=0;i<ni-1;i++)
				{
					k=nk-2;
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					index[0]=k,index[1]=j,index[2]=i;
					diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
					p=*(surf+(indexP)*6*3+3*5+0);
					q=*(surf+(indexP)*6*3+3*5+1);
					r=*(surf+(indexP)*6*3+3*5+2);
					Area = sqrt(p*p+q*q+r*r);
					Ewire = fabs(Fdiff[5])/Area;
					errorT = errorT + ((Epeek - Ewire)*(Epeek - Ewire)/(Ewire*Ewire));
				}
			}
			
			errorT = sqrt(errorT/((ni-1)*(nj-1)));
			
			if(errorT>Tolerance)
			{
				printf ("\nPercentage Error = %e\nI am updating Ion Charge..........\n",errorT);
				for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						k=nk-2;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p=*(surf+(indexP)*6*3+3*5+0);
						q=*(surf+(indexP)*6*3+3*5+1);
						r=*(surf+(indexP)*6*3+3*5+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = fabs(Fdiff[5])/Area;
						RhoWire = *(dummy+((ni-1)*j+i));
						RhoStar = RhoWire*pow((Ewire/Epeek),Gamma);
						*(dummy+((ni-1)*j+i))=(1-Omega)*RhoWire+Omega*RhoStar;
						*(dset->funct+((ni-1)*j+i)) = *(dummy+((ni-1)*j+i));
					}
					base->ChargeConv = 1;
					//fprintf(conv,"%d",1);//1 - Implies Electric field not converged
			}//end if
			else
			{
				//printf("\nI have not updated Ion Charge\n");
				printf("\nNo Error in this Iteration\n");
				base->ChargeConv = 0;
				//fprintf(conv,"%d",0);//0 - Implies Electric field converged on the wire
			}
			break;
				
			case 6://bottom
			errorT = 0.0;error = 0.0;Out = 0;Ewire = 0.0;
			for(j=0;j<nj-1;j++) 
			{
				for(i=0;i<ni-1;i++)
				{
					k=0;
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					index[0]=k,index[1]=j,index[2]=i;
					diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
					p=*(surf+(indexP)*6*3+3*4+0);
					q=*(surf+(indexP)*6*3+3*4+1);
					r=*(surf+(indexP)*6*3+3*4+2);
					Area = sqrt(p*p+q*q+r*r);
					Ewire = fabs(Fdiff[4])/Area;
					errorT = errorT + ((Epeek - Ewire)*(Epeek - Ewire)/(Ewire*Ewire));
				}
			}
			
			errorT = sqrt(errorT/((ni-1)*(nj-1)));
			
			if(errorT>Tolerance)
			{
				printf ("\nPercentage Error = %e\nI am updating Ion Charge..........\n",errorT);
				for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						k=0;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						index[0]=k,index[1]=j,index[2]=i;
						diffdia =Diffusion2d(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff,CoeffN);
						p=*(surf+(indexP)*6*3+3*4+0);
						q=*(surf+(indexP)*6*3+3*4+1);
						r=*(surf+(indexP)*6*3+3*4+2);
						Area = sqrt(p*p+q*q+r*r);
						Ewire = fabs(Fdiff[4])/Area;
						RhoWire = *(dummy+((ni-1)*j+i));
						RhoStar = RhoWire*pow((Ewire/Epeek),Gamma);
						*(dummy+((ni-1)*j+i))=(1-Omega)*RhoWire+Omega*RhoStar;
						*(dset->funct+((ni-1)*j+i)) = *(dummy+((ni-1)*j+i));
					}
					base->ChargeConv = 1;
					//fprintf(conv,"%d",1);//1 - Implies Electric field not converged
			}//end if
			else
			{
				//printf("\nI have not updated Ion Charge\n");
				printf("\nNo Error in this Iteration\n");
				base->ChargeConv = 0;
				//fprintf(conv,"%d",0);//0 - Implies Electric field converged on the wire
			}
			break;
			
			
		}
		//fclose(conv);
	}//End if
	
	
	return (stp);
}//End of main


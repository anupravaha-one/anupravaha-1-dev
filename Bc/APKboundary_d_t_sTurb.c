#include "APKboundary_d_t_sTurb.h"

/**
The function boundary_d_t_sTurb() actually Impliments the priscribed Turbulent boundary conditions on respective boundary patchs 
of all the zones. The input parameters to this function are Pointer of the base, direction index of the 
boundary patch (give ref), pointer to the varible arry of flowsol node for which BC is being implimented, 
variable tag of the variable, pointer to the zone, cell size of the boundary patch, pointer 
to the Dataset member of data structure, pointer to the variable array in boundary slab respectivly.      

This function is divided in to two major parts by an if statement, the division is based on whether the 
boundary condition to be implimented is constant or a function.

The outer most if loop is used for verious BC types viz. Drichilet, Nueman, Robin, Orlanski ( put referance).

The boundary conditions which include the gradiant of the variable ie. Nueman, Robin, Orlanski are further 
subdivided in to Non zero gradiant with orthogonal grid at the  boundary, Zero gradiant with orthogonal grid 
at the  boundary, Non zero or Zero gradiant with non-orthogonal grid at the boundary. 

In the case of non orthogonal grid, the discritized BC equatuion contains naighbouring boundary points. 
This makes the bounary conditin calculation implicit and hence the BC points are calculated iterativly. 

In all tha cases discussed above, for retrival of current data from the data structure, switch cases is used for 
six direction of a zone.
 */

double boundary_d_t_sTurb(Base* base,ZoneBoundary *zBc,int direction,double *flowsol,int vartag,Zone *z,int size,Dataset *dset,double *dummy)
{
	int ni,nj,nk,size_ind,nobc,numbc;
	int flagLogKEps;
	double p,q,r,Area,qo,res,diffCoef,Uav=1.0;
	double *vols,*beta,*surf,*deltax;
	int i,j,k;
	int *ptrsize,indexP,slabIndex=0;
	double *slabVelocity[6],*slabVolume[6],Fdiff[6];
	int index[3];
	double digonal,sum=0.00,Diff,delta_t,*slabPreviousTimeStep;
	double diagonal,relaxation=1.0;
	double Gamma,Alpha,BetaCoeff,RHS,convcri;
	
	int flagMPFV;//To implement the wall function treatment for the case of MPFV--chandrasekhar
	flagMPFV=base->flagMPFV;
	
	flagLogKEps = base->flagLogKEps;
	
	ni=*(z->Size+0);
	nj=*(z->Size+1);
	nk=*(z->Size+2); 

	delta_t=base->BaseParameters->TimeStep;
	numbc=z->Number_of_BC_Patches;

	relaxation=*(base->BaseParameters->OMEGA+0);
	if(vartag == 4)
	{
		diffCoef=1.0;
		relaxation=*(base->BaseParameters->OMEGA+1);
	}
	else if(vartag == 6)
	{
		relaxation=*(base->BaseParameters->OMEGA+2);
	}
	else if(vartag == 600)
		relaxation=*(base->BaseParameters->OMEGA+4);
	else if(vartag == 601)
		relaxation=*(base->BaseParameters->OMEGA+5);
	
	ptrsize=z->Size;
	
	if(dset->identity[0]==0)// Constant Value
	{
		if(dset->Bc_Type==10)//NON-ORTHOGONAL GRID Neumann BC
		{
			//Remember Inward flux is positive and outward flux is negative 
			qo=*(dset->Data_Array);
			MakeSlabPtrArray(z,slabVelocity,vartag);/*Called for making slabs of u */
			MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;

			convcri = find_convcri(vartag,base);
			do
			{
				slabIndex=0;
				sum=0.0;
				switch(direction)
				{
					case 1://east
						i=ni-2;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*1+0);
							q=*(surf+(indexP)*6*3+3*1+1);
							r=*(surf+(indexP)*6*3+3*1+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=ni-2;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*1+0);
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=(qo*Area)+diffCoef*Fdiff[1];//Diff;
							}
							else
								res=(qo*Area)-diffCoef*Fdiff[1];//Diff;
							
							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
					
					case 2://west
						i=0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*0+0);
							q=*(surf+(indexP)*6*3+3*0+1);
							r=*(surf+(indexP)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=0;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*0+0));
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=qo*Area-diffCoef*Fdiff[0];
							}
							else
								res=qo*Area+diffCoef*Fdiff[0];
							
							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
				
					case 3://north
						j = nj-2;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p = *(surf+(indexP)*6*3+3*3+0);
							q = *(surf+(indexP)*6*3+3*3+1);
							r = *(surf+(indexP)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=nj-2,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*3+1);
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=(qo*Area)+diffCoef*Fdiff[3];//Diff;
							}
							else
								res=(qo*Area)-diffCoef*Fdiff[3];//Diff;
							
							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
					
					case 4://south	
						j=0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*2+0);
							q=*(surf+(indexP)*6*3+3*2+1);
							r=*(surf+(indexP)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=0,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*2+1));
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=qo*Area-diffCoef*Fdiff[2];
							}
							else
								res=qo*Area+diffCoef*Fdiff[2];
				
							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
					
					case 5://top
						k=nk-2;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*5+0);
							q=*(surf+(indexP)*6*3+3*5+1);
							r=*(surf+(indexP)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=nk-2,index[1]=j,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*5+2);
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=(qo*Area)+diffCoef*Fdiff[5];//Diff;
							}
							else
								res=(qo*Area)-diffCoef*Fdiff[5];//Diff;
				
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
					
					case 6://bottom
						k=0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*4+0);
							q=*(surf+(indexP)*6*3+3*4+1);
							r=*(surf+(indexP)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=0,index[1]=j,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*4+2));
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=qo*Area-diffCoef*Fdiff[4];
							}
							else
								res=qo*Area+diffCoef*Fdiff[4];
				
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
							}
							sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
				}
				
				if(sum > 1.0e8)
				{
					printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
					exit(24);
				}
			}while(sum > convcri);
		}
		else if(dset->Bc_Type==0)//ORTHOGONAL GRID Homogeneous Neumann BC
		{
			switch(direction)
			{
				case 1:
					i=ni-2;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						*(dummy+k*(nj-1)+j)=*(flowsol+indexP);
// 						if(zBc->Boundary_Tag==22&&base->masstag==0&&(vartag==44||vartag==45||vartag==46))
// 						{
// 							*(dummy+k*(nj-1)+j)=*(dummy+k*(nj-1)+j)/(*(zBc->ZSlab->VarthetaPhase[1]+indexP)* 1.225);
// 						}
					}
					break;
				
				case 2:
					i=0;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						*(dummy+k*(nj-1)+j)=*(flowsol+indexP);
// 						if(zBc->Boundary_Tag==22&&base->masstag==0&&(vartag==44||vartag==45||vartag==46))
// 						{
// 							*(dummy+k*(nj-1)+j)=*(dummy+k*(nj-1)+j)/(*(zBc->ZSlab->VarthetaPhase[1]+indexP)* 1.225);
// 						}
					}
					break;
				
				case 3:
					j=nj-2;
					for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						*(dummy+k*(ni-1)+i)=*(flowsol+indexP);
// 						if(zBc->Boundary_Tag==22&&base->masstag==0&&(vartag==44||vartag==45||vartag==46))
// 						{
// 							*(dummy+k*(ni-1)+i)=*(dummy+k*(ni-1)+i)/(*(zBc->ZSlab->VarthetaPhase[1]+indexP)* 1.225);
// 						}
					}
					break;
				
				case 4:
					j=0;
					for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						*(dummy+k*(ni-1)+i)=*(flowsol+indexP);
// 						if(zBc->Boundary_Tag==22&&base->masstag==0&&(vartag==44||vartag==45||vartag==46))
// 						{
// 							*(dummy+k*(ni-1)+i)=*(dummy+k*(ni-1)+i)/(*(zBc->ZSlab->VarthetaPhase[1]+indexP)* 1.225);
// 						}
					}
					break;
				
				case 5:
					k=nk-2;
					for (j=0;j<nj-1;j++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						*(dummy+j*(ni-1)+i)=*(flowsol+indexP);
// 						if(zBc->Boundary_Tag==22&&base->masstag==0&&(vartag==44||vartag==45||vartag==46))
// 						{
// 							*(dummy+j*(ni-1)+i)=*(dummy+j*(ni-1)+i)/(*(zBc->ZSlab->VarthetaPhase[1]+indexP)* 1.225);
// 						}
					}
					break;
				
				case 6:
					k=0;
					for (j=0;j<nj-1;j++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						*(dummy+j*(ni-1)+i)=*(flowsol+indexP);
// 						if(zBc->Boundary_Tag==22&&base->masstag==0&&(vartag==44||vartag==45||vartag==46))
// 						{
// 							*(dummy+j*(ni-1)+i)=*(dummy+j*(ni-1)+i)/(*(zBc->ZSlab->VarthetaPhase[1]+indexP)* 1.225);
// 						}
					}
					break;
			}
			//if(vartag==31)
				//printf("\n *(dummy)=%lf",*(dummy));
		}
		else if(dset->Bc_Type==1)//ORTHOGONAL GRID non-homogeneous Neumann BC
		{
			vols=z->ZGeom->Volume;
			beta=z->ZGeom->Beta;
			surf=z->ZGeom->SurfaceArea;
			qo=*(dset->Data_Array);
			switch(direction)
			{
				case 1:
					i=ni-2;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*1+0);
						q=*(surf+(indexP)*6*3+3*1+1);
						r=*(surf+(indexP)*6*3+3*1+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*1+0))*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
				
				case 2://west
					i=0;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*0+0))*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
			
				case 3://north
					j=nj-2;
					for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*3+0);
						q=*(surf+(indexP)*6*3+3*3+1);
						r=*(surf+(indexP)*6*3+3*3+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=*(beta+(indexP)*6*3+3*3+1)*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
			
				case 4://south
					j=0;
					for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*2+0);
						q=*(surf+(indexP)*6*3+3*2+1);
						r=*(surf+(indexP)*6*3+3*2+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*2+1))*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
			
				case 5://top
					k=nk-2;
					for (j=0;j<nj-1;j++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*5+0);
						q=*(surf+(indexP)*6*3+3*5+1);
						r=*(surf+(indexP)*6*3+3*5+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=*(beta+(indexP)*6*3+3*5+2)*diffCoef;
						if(digonal <= 0.0)
						digonal=-digonal;

						*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;

				case 6://bottom
					k=0;
					for (j=0;j<nj-1;j++) 
					for (i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*4+0);
						q=*(surf+(indexP)*6*3+3*4+1);
						r=*(surf+(indexP)*6*3+3*4+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*4+2))*diffCoef;
						if(digonal <= 0.0)
						digonal=-digonal;

						*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
			}
		}
		else if(dset->Bc_Type==3)//Homogeneous Dirichlet BC
		{
			if((flagLogKEps == 1) && ((vartag==600)||(vartag==601)))
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind) = WallFunctionDirichlet(z,zBc,vartag,size_ind,*(dset->Data_Array),flagMPFV);
					*(dummy+size_ind) = log(*(dummy+size_ind));
				}
			}
			else
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind) = WallFunctionDirichlet(z,zBc,vartag,size_ind,*(dset->Data_Array),flagMPFV);
				}
			}
		}
		else if(dset->Bc_Type==2)//Dirichlet BC
		{
			if((flagLogKEps == 1) && ((vartag==600)||(vartag==601)))
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind) = WallFunctionDirichlet(z,zBc,vartag,size_ind,*(dset->Data_Array),flagMPFV);
					*(dummy+size_ind) = log(*(dummy+size_ind));
					//printf("Boundary Dirichlet value= %lf vartag=%d Dir=%d\n",*(dummy+size_ind),vartag,direction);
				}
			}
			else
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind) = WallFunctionDirichlet(z,zBc,vartag,size_ind,*(dset->Data_Array),flagMPFV);
					
				}
			}
		}
		else if(dset->Bc_Type==41)// Robin BC for non-orthogonal grid
		{
			Gamma=*(dset->Data_Array+2);
			Alpha=*(dset->Data_Array+0);
			BetaCoeff=*(dset->Data_Array+1);
			
			MakeSlabPtrArray(z,slabVelocity,vartag);//Called for making slabs of u 
			MakeSlabPtrArray(z,slabVolume,5);//Called for making slabs of Volumes

			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;
			convcri = find_convcri(vartag,base);
			do
			{
				slabIndex=0;
				sum=0.0;
				switch(direction)
				{
					case 1://east
						i=ni-2;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*1+0);
							q=*(surf+(indexP)*6*3+3*1+1);
							r=*(surf+(indexP)*6*3+3*1+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=ni-2;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=Gamma*Area-BetaCoeff*Fdiff[1]-Alpha*Area*(*(dummy+((nj-1)*k+j)));//Diff;
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*1+0))+Alpha*Area;

							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
					
					case 2://west
						i=0;
						for(k=0;k<nk-1;k++)
						for(j=0;j<nj-1;j++)
						{
							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*0+0);
							q=*(surf+(indexP)*6*3+3*0+1);
							r=*(surf+(indexP)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=0;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=Gamma*Area-BetaCoeff*Fdiff[0]-Alpha*Area*(*(dummy+((nj-1)*k+j)));
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*0+0))+Alpha*Area;

							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
			
					case 3://north
						j = nj-2;
						for(k=0;k<nk-1;k++) 
						for(i=0;i<ni-1;i++)
						{
							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p = *(surf+(indexP)*6*3+3*3+0);
							q = *(surf+(indexP)*6*3+3*3+1);
							r = *(surf+(indexP)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=nj-2,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=(Gamma*Area)-BetaCoeff*Fdiff[3]-Alpha*Area*(*(dummy+((ni-1)*k+i)));//Diff;
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*3+1))+Alpha*Area;

							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
					
					case 4://south	
						j=0;
						for(k=0;k<nk-1;k++) 
						for(i=0;i<ni-1;i++)
						{
							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*2+0);
							q=*(surf+(indexP)*6*3+3*2+1);
							r=*(surf+(indexP)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=0,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=Gamma*Area-BetaCoeff*Fdiff[2]-Alpha*Area*(*(dummy+((ni-1)*k+i)));
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*2+1))+Alpha*Area;

							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(digonal);
							sum=sum+res*res;

						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
					
					case 5://top
						k=nk-2;
						for(j=0;j<nj-1;j++) 
						for(i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*5+0);
							q=*(surf+(indexP)*6*3+3*5+1);
							r=*(surf+(indexP)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=nk-2,index[1]=j,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=(Gamma*Area)-BetaCoeff*Fdiff[5]-Alpha*Area*(*(dummy+((ni-1)*j+i)));
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*5+2))+Alpha*Area;

							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
					
					case 6://bottom
						k=0;
						for(j=0;j<nj-1;j++) 
						for(i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							index[0]=0,index[1]=j,index[2]=i;
							p=*(surf+(indexP)*6*3+3*4+0);
							q=*(surf+(indexP)*6*3+3*4+1);
							r=*(surf+(indexP)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
	
							res=Gamma*Area-BetaCoeff*Fdiff[4]-Alpha*Area*(*(dummy+((ni-1)*j+i)));
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*4+2))+Alpha*Area;
	
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
				}
					
				if(sum > 1.0e8)
				{
					printf("Robin BC: variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
					exit(24);
				}
			}while(sum > convcri);
		}
		else if(dset->Bc_Type==40)// Robin BC for orthogonal grid
		{
			Gamma=*(dset->Data_Array+2);
			Alpha=*(dset->Data_Array+0);
			BetaCoeff=*(dset->Data_Array+1);
			
			MakeSlabPtrArray(z,slabVelocity,vartag);//Called for making slabs of u 
			MakeSlabPtrArray(z,slabVolume,5);//Called for making slabs of Volumes
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;
			
			switch(direction)
			{
				case 1://east
					i=ni-2;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*1+0);
						q=*(surf+(indexP)*6*3+3*1+1);
						r=*(surf+(indexP)*6*3+3*1+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=j,index[2]=ni-2;

						Fdiff[1]=(*(beta+(indexP)*6*3+3*1+0))*(*(flowsol+indexP));
						res=Gamma*Area+BetaCoeff*Fdiff[1];//Diff;
						digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*1+0))+Alpha*Area;

						*(dummy+((nj-1)*k+j))=res/(digonal);
					}
					break;
				
				case 2://west
					i=0;
					for(k=0;k<nk-1;k++)
					for(j=0;j<nj-1;j++)
					{
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=j,index[2]=0;

						Fdiff[0]=(*(beta+(indexP)*6*3+3*0+0))*(*(flowsol+indexP));
						res=Gamma*Area-BetaCoeff*Fdiff[0];
						digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*0+0))+Alpha*Area;

						*(dummy+((nj-1)*k+j))=res/(digonal);
					}
					break;
				
				case 3://north
					j = nj-2;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p = *(surf+(indexP)*6*3+3*3+0);
						q = *(surf+(indexP)*6*3+3*3+1);
						r = *(surf+(indexP)*6*3+3*3+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=nj-2,index[2]=i;

						Fdiff[3]=(*(beta+(indexP)*6*3+3*3+1))*(*(flowsol+indexP));
						res=(Gamma*Area)+BetaCoeff*Fdiff[3];
						digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*3+1))+Alpha*Area;

						*(dummy+((ni-1)*k+i))=res/(digonal);
					}
					break;
				
				case 4://south	
					j=0;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*2+0);
						q=*(surf+(indexP)*6*3+3*2+1);
						r=*(surf+(indexP)*6*3+3*2+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=0,index[2]=i;

						Fdiff[2]=(*(beta+(indexP)*6*3+3*2+1))*(*(flowsol+indexP));
						res=Gamma*Area-BetaCoeff*Fdiff[2];
						digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*2+1))+Alpha*Area;

						*(dummy+((ni-1)*k+i))=res/(digonal);
					}
					break;
				
				case 5://top
					k=nk-2;
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*5+0);
						q=*(surf+(indexP)*6*3+3*5+1);
						r=*(surf+(indexP)*6*3+3*5+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=nk-2,index[1]=j,index[2]=i;

						Fdiff[5]=(*(beta+(indexP)*6*3+3*5+2))*(*(flowsol+indexP));
						res=(Gamma*Area)+BetaCoeff*Fdiff[5];
						digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*5+2))+Alpha*Area;

						*(dummy+((ni-1)*j+i))=res/(digonal);
					}
					break;
				
				case 6://bottom
					k=0;
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*4+0);
						q=*(surf+(indexP)*6*3+3*4+1);
						r=*(surf+(indexP)*6*3+3*4+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=0,index[1]=j,index[2]=i;

						Fdiff[4]=(*(beta+(indexP)*6*3+3*4+2))*(*(flowsol+indexP));
						res=Gamma*Area-BetaCoeff*Fdiff[4];
						digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*4+2))+Alpha*Area;

						*(dummy+((ni-1)*j+i))=res/(digonal);
					}
					break;
			}
		}
		else if(dset->Bc_Type==5)//Orlanski BC <<< Not working need to be worked out ... Narsimha 8th Feb, 2008
		{
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;
			switch(direction)
			{
				case 1:
					deltax= z->ZGeom->deltaX;
						
					for(nobc=0;nobc<numbc;nobc++)
					{ 
						if(direction==(z->ZoneBc+nobc)->Direction)
							slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
					}
						
					for(k=0;k<nk-1;k++)
					{
						for(j=0;j<nj-1;j++)
						{
							Fdiff[1]=*(flowsol+(ni-1)*(nj-1)*k+(ni-1)*j+(ni-2))/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*j+(ni-2))*6*3+3*1+0));
							RHS=*(slabPreviousTimeStep+((nj-1)*k+j))+(Uav*delta_t*Fdiff[1]);
							diagonal=1.0+Uav*delta_t/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*j+(ni-2))*6*3+3*1+0));
							*(dummy+((nj-1)*k+j))=RHS/diagonal;
						}
					}
					break;
				case 2:
					deltax= z->ZGeom->deltaX;
						
					for(nobc=0;nobc<numbc;nobc++)
					{ 
						if(direction==(z->ZoneBc+nobc)->Direction)
							slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
					}
						
					for (k=0;k<nk-1;k++) 
					{
						for (j=0;j<nj-1;j++)
						{
							Fdiff[0]=*(flowsol+(ni-1)*(nj-1)*k+(ni-1)*j)/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*j+0)*6*3+3*0+0));
							RHS=*(slabPreviousTimeStep+((nj-1)*k+j))+(Uav*delta_t*Fdiff[0]);//Diff;
							diagonal=1+Uav*delta_t/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*j+0)*6*3+3*0+0));
							*(dummy+((nj-1)*k+j))=RHS/diagonal;
						}
					}
					break;
				case 3:
					deltax= z->ZGeom->deltaX;	
						
					for(nobc=0;nobc<numbc;nobc++)
					{ 
						if(direction==(z->ZoneBc+nobc)->Direction)
							slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
					}
						
					for (k=0;k<nk-1;k++) 
					{
						for (i=0;i<ni-1;i++)
						{
							Fdiff[3]=*(flowsol+(ni-1)*(nj-1)*k+(ni-1)*(nj-2)+i)/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*(nj-2)+i)*6*3+3*3+1));
							RHS=*(slabPreviousTimeStep+(ni-1)*k+i)+(Uav*delta_t*Fdiff[3]);//Diff;
							diagonal=1+Uav*delta_t/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*(nj-2)+i)*6*3+3*3+1));
							*(dummy+((ni-1)*k+i))=RHS/diagonal;
						}
					}
						
					break;
				case 4:
					for(nobc=0;nobc<numbc;nobc++)
					{ 
						if(direction==(z->ZoneBc+nobc)->Direction)
							slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
					}
						
					for (k=0;k<nk-1;k++) 
					{
						for (i=0;i<ni-1;i++)
						{
							
							Fdiff[2]=*(flowsol+(ni-1)*(nj-1)*k+(ni-1)*0+i)/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*0+i)*6*3+3*2+1));
							RHS=*(slabPreviousTimeStep+(ni-1)*k+i)+(Uav*delta_t*Fdiff[2]);//Diff;
							diagonal=1+Uav*delta_t/(*(deltax+((ni-1)*(nj-1)*k+(ni-1)*0+i)*6*3+3*2+1));
							*(dummy+((ni-1)*k+i))=RHS/diagonal;
						}
					}
					break;	
					
				case 5:
					for(nobc=0;nobc<numbc;nobc++)
					{ 
						if(direction==(z->ZoneBc+nobc)->Direction)
							slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
					}
						
					for (j=0;j<nj-1;j++) 
					{
						for (i=0;i<ni-1;i++)
						{
							Fdiff[5]=*(flowsol+(ni-1)*(nj-1)*(nk-2)+(ni-1)*j+i)/(*(deltax+((ni-1)*(nj-1)*(nk-2)+(ni-1)*j+i)*6*3+3*5+2));
							RHS=*(slabPreviousTimeStep+((ni-1)*j+i))+(Uav*delta_t*Fdiff[5]);//Diff;
							digonal=1+Uav*delta_t/(*(deltax+((ni-1)*(nj-1)*(nk-2)+(ni-1)*j+i)*6*3+3*5+2));
							*(dummy+((ni-1)*j+i))=RHS/diagonal;
						}
					}
					break;
				case 6:
					for(nobc=0;nobc<numbc;nobc++)
					{ 
						if(direction==(z->ZoneBc+nobc)->Direction)
							slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
					}
						
					for (j=0;j<nj-1;j++) 
					{
						for (i=0;i<ni-1;i++)
						{
							
							Fdiff[4]=*(flowsol+(ni-1)*(nj-1)*0+(ni-1)*j+i)/(*(deltax+((ni-1)*(nj-1)*0+(ni-1)*j+i)*6*3+3*4+2));
							RHS=*(slabPreviousTimeStep+((ni-1)*j+i))+(Uav*delta_t*Fdiff[4]);//Diff;
							digonal=1+Uav*delta_t/(*(deltax+((ni-1)*(nj-1)*0+(ni-1)*j+i)*6*3+3*4+2));
							*(dummy+((ni-1)*j+i))=RHS/diagonal;
						}
					}
					break;
			}
		}
	}
	if(dset->identity[0]==1)// If Given Function Then it will Activate
	{
		if(dset->Bc_Type==2)//Drichlet BC
		{
			if((flagLogKEps == 1) && ((vartag==600)||(vartag==601)))
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind) = WallFunctionDirichlet(z,zBc,vartag,size_ind,*(dset->funct+size_ind),flagMPFV);
					*(dummy+size_ind) = log(*(dummy+size_ind));
				}
			}
			else
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind) = WallFunctionDirichlet(z,zBc,vartag,size_ind,*(dset->funct+size_ind),flagMPFV);
					
				}
			}
		}
		else if(dset->Bc_Type==10)//NON-ORTHOGONAL GRID Neumann BC <<< Implementation is to be verified ... Narsimha 8th Feb, 2008
		{
			MakeSlabPtrArray(z,slabVelocity,vartag);/*Called for making slabs of u */
			MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;

			convcri = find_convcri(vartag,base);
			do
			{
				slabIndex=0;
				sum=0.0;
				switch(direction)
				{
					case 1://east
						i=ni-2;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							qo=*((dset->funct)+(nj-1)*k+j);
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*1+0);
							q=*(surf+(indexP)*6*3+3*1+1);
							r=*(surf+(indexP)*6*3+3*1+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=ni-2;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							if (dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(nj-1)*k+j);

							if (dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(nj-1)*k+j);

							if (dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(nj-1)*k+j);

							digonal=*(beta+(indexP)*6*3+3*1+0);
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=(qo*Area)+diffCoef*Fdiff[1];//Diff;
							}
							else
								res=(qo*Area)-diffCoef*Fdiff[1];//Diff;
								
							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
					
					case 2://west
						i=0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							qo=*((dset->funct)+(nj-1)*k+j);
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*0+0);
							q=*(surf+(indexP)*6*3+3*0+1);
							r=*(surf+(indexP)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=0;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*0+0));
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=qo*Area-diffCoef*Fdiff[0];
							}
							else
								res=qo*Area+diffCoef*Fdiff[0];
							
							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
					
					case 3://north
						j = nj-2;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							qo=*((dset->funct)+(ni-1)*k+i);
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p = *(surf+(indexP)*6*3+3*3+0);
							q = *(surf+(indexP)*6*3+3*3+1);
							r = *(surf+(indexP)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=nj-2,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*3+1);
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=(qo*Area)+diffCoef*Fdiff[3];//Diff;
							}
							else
								res=(qo*Area)-diffCoef*Fdiff[3];//Diff;

							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
					
					case 4://south	
						j=0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							qo=*((dset->funct)+(ni-1)*k+i);
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*2+0);
							q=*(surf+(indexP)*6*3+3*2+1);
							r=*(surf+(indexP)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=0,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*2+1));
							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=qo*Area-diffCoef*Fdiff[2];
							}
							else
								res=qo*Area+diffCoef*Fdiff[2];
				
							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
					
					case 5://top
						k=nk-2;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							qo=*((dset->funct)+(ni-1)*j+i);
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*5+0);
							q=*(surf+(indexP)*6*3+3*5+1);
							r=*(surf+(indexP)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=nk-2,index[1]=j,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*5+2);

							if(digonal <= 0.0)
							{
								digonal=-digonal;
								res=(qo*Area)+diffCoef*Fdiff[5];//Diff;
							}
							else
							res=(qo*Area)-diffCoef*Fdiff[5];//Diff;
				
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
					
					case 6://bottom
						k=0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							qo=*((dset->funct)+(ni-1)*j+i);
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							p=*(surf+(indexP)*6*3+3*4+0);
							q=*(surf+(indexP)*6*3+3*4+1);
							r=*(surf+(indexP)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=0,index[1]=j,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*4+2));
							if(digonal <= 0.0)
							{
							digonal=-digonal;
							res=qo*Area-diffCoef*Fdiff[4];
							}
							else
							res=qo*Area+diffCoef*Fdiff[4];
				
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
				}
					
				if(sum > 1.0e8)
				{
					printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
					exit(24);
				}
			}while(sum > convcri);
		}
		else if(dset->Bc_Type==1)//ORTHOGONAL GRID non-homogeneous Neumann BC
		{
			vols=z->ZGeom->Volume;
			beta=z->ZGeom->Beta;
			surf=z->ZGeom->SurfaceArea;
			switch(direction)
			{
				case 1:
					i=ni-2;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						qo=*((dset->funct)+(nj-1)*k+j);
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*1+0);
						q=*(surf+(indexP)*6*3+3*1+1);
						r=*(surf+(indexP)*6*3+3*1+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*1+0))*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
				
				case 2://west
					i=0;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						qo=*((dset->funct)+(nj-1)*k+j);
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*0+0))*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
				
				case 3://north
					j=nj-2;
					for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						qo=*((dset->funct)+(ni-1)*k+i);
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*3+0);
						q=*(surf+(indexP)*6*3+3*3+1);
						r=*(surf+(indexP)*6*3+3*3+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=*(beta+(indexP)*6*3+3*3+1)*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
				
				case 4://south
					j=0;
					for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{
						qo=*((dset->funct)+(ni-1)*k+i);
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*2+0);
						q=*(surf+(indexP)*6*3+3*2+1);
						r=*(surf+(indexP)*6*3+3*2+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*2+1))*diffCoef;
						if(digonal <= 0.0)
							digonal=-digonal;

						*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
				
				case 5://top
					k=nk-2;
					for (j=0;j<nj-1;j++) 
					for (i=0;i<ni-1;i++)
					{
						qo=*((dset->funct)+(ni-1)*j+i);
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*5+0);
						q=*(surf+(indexP)*6*3+3*5+1);
						r=*(surf+(indexP)*6*3+3*5+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=*(beta+(indexP)*6*3+3*5+2)*diffCoef;
						if(digonal <= 0.0)
						digonal=-digonal;

						*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
					
				case 6://bottom
					k=0;
					for (j=0;j<nj-1;j++) 
					for (i=0;i<ni-1;i++)
					{
						qo=*((dset->funct)+(ni-1)*j+i);
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
						p=*(surf+(indexP)*6*3+3*4+0);
						q=*(surf+(indexP)*6*3+3*4+1);
						r=*(surf+(indexP)*6*3+3*4+2);
						Area=sqrt(p*p+q*q+r*r);
						digonal=(*(beta+(indexP)*6*3+3*4+2))*diffCoef;
						if(digonal <= 0.0)
						digonal=-digonal;

						*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo*Area)/digonal;
					}
					break;
			}
		}
		else if(dset->Bc_Type==41||dset->Bc_Type==71)//Robin OR Time dependent BC for NON-ORTHOGONAL GRID
		{
			MakeSlabPtrArray(z,slabVelocity,vartag);//Called for making slabs of u 
			MakeSlabPtrArray(z,slabVolume,5);//Called for making slabs of Volumes

			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;
			convcri = find_convcri(vartag,base);
			do
			{
				slabIndex=0;
				sum=0.0;
				switch(direction)
				{
					case 1://east
						i=ni-2;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							if(dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(nj-1)*k+j); 
	
							if(dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(nj-1)*k+j);
	
							if(dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(nj-1)*k+j);
	
							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*1+0);
							q=*(surf+(indexP)*6*3+3*1+1);
							r=*(surf+(indexP)*6*3+3*1+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=ni-2;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=Gamma*Area-BetaCoeff*Fdiff[1]-Alpha*Area*(*(dummy+((nj-1)*k+j)));//Diff;
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*1+0))+Alpha*Area;

							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
	
					case 2://west
						i=0;
						for(k=0;k<nk-1;k++)
						for(j=0;j<nj-1;j++)
						{
							if(dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(nj-1)*k+j); 

							if(dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(nj-1)*k+j);

							if(dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(nj-1)*k+j);

							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*0+0);
							q=*(surf+(indexP)*6*3+3*0+1);
							r=*(surf+(indexP)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=j,index[2]=0;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=Gamma*Area-BetaCoeff*Fdiff[0]-Alpha*Area*(*(dummy+((nj-1)*k+j)));
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*0+0))+Alpha*Area;

							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						break;
	
					case 3://north
						j = nj-2;
						for(k=0;k<nk-1;k++) 
						for(i=0;i<ni-1;i++)
						{
							if(dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(ni-1)*k+i); 

							if(dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(ni-1)*k+i);

							if(dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(ni-1)*k+i);

							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p = *(surf+(indexP)*6*3+3*3+0);
							q = *(surf+(indexP)*6*3+3*3+1);
							r = *(surf+(indexP)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=nj-2,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=(Gamma*Area)-BetaCoeff*Fdiff[3]-Alpha*Area*(*(dummy+((ni-1)*k+i)));//Diff;
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*3+1))+Alpha*Area;

							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
	
					case 4://south	
						j=0;
						for(k=0;k<nk-1;k++) 
						for(i=0;i<ni-1;i++)
						{
							if(dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(ni-1)*k+i); 

							if(dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(ni-1)*k+i);

							if(dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(ni-1)*k+i);

							indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*2+0);
							q=*(surf+(indexP)*6*3+3*2+1);
							r=*(surf+(indexP)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k,index[1]=0,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=Gamma*Area-BetaCoeff*Fdiff[2]-Alpha*Area*(*(dummy+((ni-1)*k+i)));
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*2+1))+Alpha*Area;

							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						break;
	
					case 5://top
						k=nk-2;
						for(j=0;j<nj-1;j++) 
						for(i=0;i<ni-1;i++)
						{
							if(dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(ni-1)*j+i); 

							if(dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(ni-1)*j+i);

							if(dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(ni-1)*j+i);

							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*5+0);
							q=*(surf+(indexP)*6*3+3*5+1);
							r=*(surf+(indexP)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=nk-2,index[1]=j,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);

							res=(Gamma*Area)-BetaCoeff*Fdiff[5]-Alpha*Area*(*(dummy+((ni-1)*j+i)));
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*5+2))+Alpha*Area;

							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
	
					case 6://bottom
						k=0;
						for(j=0;j<nj-1;j++) 
						for(i=0;i<ni-1;i++)
						{
							if(dset->identity[1]==0)
								Alpha=*(dset->Data_Array+0);
							else
								Alpha=*(dset->funct+(ni-1)*j+i); 
	
							if(dset->identity[2]==0)
								BetaCoeff=*(dset->Data_Array+1);
							else
								BetaCoeff=*(dset->funct_beta+(ni-1)*j+i);
	
							if(dset->identity[3]==0)
								Gamma=*(dset->Data_Array+2);
							else
								Gamma=*(dset->funct_gama+(ni-1)*j+i);
	
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							p=*(surf+(indexP)*6*3+3*4+0);
							q=*(surf+(indexP)*6*3+3*4+1);
							r=*(surf+(indexP)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=0,index[1]=j,index[2]=i;
							Diff=DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
	
							res=Gamma*Area-BetaCoeff*Fdiff[4]-Alpha*Area*(*(dummy+((ni-1)*j+i)));
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*4+2))+Alpha*Area;
	
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						break;
				}	

				if(sum > 1.0e8)
				{
					printf("Robin/TimeDependent BC: variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
					exit(24);
				}
			}while(sum > convcri);
		}
		else if(dset->Bc_Type==40||dset->Bc_Type==70)//Robin OR Time dependent BC for ORTHOGONAL GRID
		{
			MakeSlabPtrArray(z,slabVelocity,vartag);//Called for making slabs of u 
			MakeSlabPtrArray(z,slabVolume,5);//Called for making slabs of Volumes
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;
			
			switch(direction)
			{
				case 1://east
					i=ni-2;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
						if (dset->identity[1]==0)
							Alpha=*(dset->Data_Array+0);
						else
							Alpha=*(dset->funct+(nj-1)*k+j);
	
						if (dset->identity[2]==0)
							BetaCoeff=*(dset->Data_Array+1);
						else
							BetaCoeff=*(dset->funct_beta+(nj-1)*k+j);
	
						if (dset->identity[3]==0)
							Gamma=*(dset->Data_Array+2);
						else
							Gamma=*(dset->funct_gama+(nj-1)*k+j);

						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*1+0);
						q=*(surf+(indexP)*6*3+3*1+1);
						r=*(surf+(indexP)*6*3+3*1+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=j,index[2]=ni-2;

						Fdiff[1]=(*(beta+(indexP)*6*3+3*1+0))*(*(flowsol+indexP));
						res=Gamma*Area+BetaCoeff*Fdiff[1];//Diff;
						digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*1+0))+Alpha*Area;

						*(dummy+((nj-1)*k+j))=res/(digonal);
					}
					break;
				
				case 2://west
					i=0;
					for(k=0;k<nk-1;k++)
					for(j=0;j<nj-1;j++)
					{
						if (dset->identity[1]==0)
							Alpha=*(dset->Data_Array+0);
						else
							Alpha=*(dset->funct+(nj-1)*k+j);

						if (dset->identity[2]==0)
							BetaCoeff=*(dset->Data_Array+1);
						else
							BetaCoeff=*(dset->funct_beta+(nj-1)*k+j);

						if (dset->identity[3]==0)
							Gamma=*(dset->Data_Array+2);
						else
							Gamma=*(dset->funct_gama+(nj-1)*k+j);

						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=j,index[2]=0;

						Fdiff[0]=(*(beta+(indexP)*6*3+3*0+0))*(*(flowsol+indexP));
						res=Gamma*Area-BetaCoeff*Fdiff[0];
						digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*0+0))+Alpha*Area;

						*(dummy+((nj-1)*k+j)) = res/(digonal);
					}
					break;
				
				case 3://north
					j = nj-2;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						if (dset->identity[1]==0)
							Alpha=*(dset->Data_Array+0);
						else
							Alpha=*(dset->funct+(ni-1)*k+i);

						if (dset->identity[2]==0)
							BetaCoeff=*(dset->Data_Array+1);
						else
							BetaCoeff=*(dset->funct_beta+(ni-1)*k+i);

						if (dset->identity[3]==0)
							Gamma=*(dset->Data_Array+2);
						else
							Gamma=*(dset->funct_gama+(ni-1)*k+i);

						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p = *(surf+(indexP)*6*3+3*3+0);
						q = *(surf+(indexP)*6*3+3*3+1);
						r = *(surf+(indexP)*6*3+3*3+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=nj-2,index[2]=i;

						Fdiff[3]=(*(beta+(indexP)*6*3+3*3+1))*(*(flowsol+indexP));
						res=(Gamma*Area)+BetaCoeff*Fdiff[3];
						digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*3+1))+Alpha*Area;

						*(dummy+((ni-1)*k+i)) = res/(digonal);
					}
					break;
				
				case 4://south	
					j=0;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						if (dset->identity[1]==0)
							Alpha=*(dset->Data_Array+0);
						else
							Alpha=*(dset->funct+(ni-1)*k+i);

						if (dset->identity[2]==0)
							BetaCoeff=*(dset->Data_Array+1);
						else
							BetaCoeff=*(dset->funct_beta+(ni-1)*k+i);

						if (dset->identity[3]==0)
							Gamma=*(dset->Data_Array+2);
						else
							Gamma=*(dset->funct_gama+(ni-1)*k+i);

						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*2+0);
						q=*(surf+(indexP)*6*3+3*2+1);
						r=*(surf+(indexP)*6*3+3*2+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=0,index[2]=i;

						Fdiff[2]=(*(beta+(indexP)*6*3+3*2+1))*(*(flowsol+indexP));
						res=Gamma*Area-BetaCoeff*Fdiff[2];
						digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*2+1))+Alpha*Area;

						*(dummy+((ni-1)*k+i)) = res/(digonal);
					}
					break;
	
				case 5://top
					k=nk-2;
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						if (dset->identity[1]==0)
						Alpha=*(dset->Data_Array+0);
						else
						Alpha=*(dset->funct+(ni-1)*j+i);

						if (dset->identity[2]==0)
						BetaCoeff=*(dset->Data_Array+1);
						else
						BetaCoeff=*(dset->funct_beta+(ni-1)*j+i);

						if (dset->identity[3]==0)
						Gamma=*(dset->Data_Array+2);
						else
						Gamma=*(dset->funct_gama+(ni-1)*j+i);
		
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*5+0);
						q=*(surf+(indexP)*6*3+3*5+1);
						r=*(surf+(indexP)*6*3+3*5+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=nk-2,index[1]=j,index[2]=i;

						Fdiff[5]=(*(beta+(indexP)*6*3+3*5+2))*(*(flowsol+indexP));
						res=(Gamma*Area)+BetaCoeff*Fdiff[5];
						digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*5+2))+Alpha*Area;

						*(dummy+((ni-1)*j+i)) = res/(digonal);
					}
					break;
				
				case 6://bottom
					k=0;
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						if (dset->identity[1]==0)
						Alpha=*(dset->Data_Array+0);
						else
						Alpha=*(dset->funct+(ni-1)*j+i);

						if (dset->identity[2]==0)
						BetaCoeff=*(dset->Data_Array+1);
						else
						BetaCoeff=*(dset->funct_beta+(ni-1)*j+i);

						if (dset->identity[3]==0)
						Gamma=*(dset->Data_Array+2);
						else
						Gamma=*(dset->funct_gama+(ni-1)*j+i);
	
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						p=*(surf+(indexP)*6*3+3*4+0);
						q=*(surf+(indexP)*6*3+3*4+1);
						r=*(surf+(indexP)*6*3+3*4+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=0,index[1]=j,index[2]=i;

						Fdiff[4]=(*(beta+(indexP)*6*3+3*4+2))*(*(flowsol+indexP));
						res=Gamma*Area-BetaCoeff*Fdiff[4];
						digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*4+2))+Alpha*Area;

						*(dummy+((ni-1)*j+i)) = res/(digonal);
					}
					break;
			}
		}
	}
//Added by chandrasekhar --- To caluculate the momentum from the updated velocity once velocity updated with in time loop
//(Except in Refresh surface cells) that's why one more variable masstag is taken.
	if((base->flagMPFV==1)&&(base->masstag==1)&&(base->flagTurb==1)&&((vartag==3100)||(vartag==3200)||(vartag==3300)||(vartag==31)||(vartag==32)||(vartag==33)||(vartag==44)||(vartag==45)||(vartag==46)))
	{
		//printf("\nI am entering update_boundaryTurbDMPFV and vartag=%d",vartag);
		update_boundaryTurbDMPFV(base,z,zBc,vartag);
			
	}
//end--chandrasekhar

	return (sum);
}


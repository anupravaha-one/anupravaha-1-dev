#include "APKboundary_d_t_s.h"

/**
The function boundary_d_t_s() actualy Impliments the pricribed boundary conditions on respective boundary patchs 
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

double boundary_d_t_s(Base* base ,int direction,double *flowsol,int vartag,Zone *z,int size,Dataset *dset,double *dummy)
{
  int ni,nj,nk,size_ind,nobc,numbc,bc_no,flagrho,flagelectcond;
	double p,q,r,Area,qo,res,diffCoef,rho,sigma,vol,Uav=1.0,Dphi[3];
	double nx,ny,nz,bx,by,bz,btx,bty,btz,bn,deltaP,kH,dPx,dPy,dPz,EpX,EpY,EpZ,dEpX,dEpY,dEpZ,Jwall,dJ,DbX,DbY,DbZ;
	double *vols,*beta,*surf,*deltax;
	int i,j,k;
	int *ptrsize,indexP,slabIndex=0;
	double *slabVelocity[6],*slabVolume[6],Fdiff[6];
	int index[3];
	double digonal,sum=0.00,Diff,delta_t,*slabPreviousTimeStep;
	double diagonal,relaxation=1.0;
	double Gamma,Alpha,BetaCoeff,RHS,convcri;
	ZoneBoundary *zBc;


	ni=*(z->Size+0);
	nj=*(z->Size+1);
	nk=*(z->Size+2); 

	//printf("vartag=%d\n",vartag);
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
	
	for(bc_no=0;bc_no<z->Number_of_BC_Patches;bc_no++)
	if((z->ZoneBc+bc_no)->Direction==direction)
	zBc=(z->ZoneBc+bc_no);
	
	if(dset->identity[0]==0)// Constant Value
	{
		if(dset->Bc_Type==10)//NON-ORTHOGONAL GRID Neumann BC
		{
			// Remember Inward flux is positive and outward flux is negative ... OK By Ashish
			qo  = *(dset->Data_Array);
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
			surf= z->ZGeom->SurfaceArea;

			MakeSlabPtrArray(z,slabVelocity,vartag);/*Called for making slabs of u */
			MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/

			convcri = find_convcri(vartag,base);
			switch(direction)
			{
				case 1://east
					i=ni-2;
					do
					{
						slabIndex=0;
						sum=0.0;
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
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 2://west
					i=0;
					do
					{
						slabIndex=0;
						sum=0.0;
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
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 3://north
					j = nj-2;
					do
					{
						slabIndex=0;
						sum=0.0;
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
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 4://south	
					j=0;
					do
					{
						slabIndex=0;
						sum=0.0;
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
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 5://top
					k=nk-2;
					do
					{
						slabIndex=0;
						sum=0.0;
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
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 6://bottom
					k=0;
					do
					{
						slabIndex=0;
						sum=0.0;
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
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
			}
		}
		else if(dset->Bc_Type==0)//ORTHOGONAL GRID Homogeneous Neumann BC ... Put this if condition SECOND FROM TOP which is used often <<< Narsimha
		{
			if(vartag!=16)
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
						}
					break;
					
					case 2:
						i=0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							*(dummy+k*(nj-1)+j)=*(flowsol+indexP);
						}
					break;
					
					case 3:
						j=nj-2;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							*(dummy+k*(ni-1)+i)=*(flowsol+indexP);
						}
					break;
					
					case 4:
						j=0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							*(dummy+k*(ni-1)+i)=*(flowsol+indexP);
						}
					break;
					
					case 5:
						k=nk-2;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							*(dummy+j*(ni-1)+i)=*(flowsol+indexP);
						}
					break;
					
					case 6:
						k=0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							*(dummy+j*(ni-1)+i)=*(flowsol+indexP);
						}
					break;
				}
			
			}
			else
			{
				//printf("\nRadiation solution, symmetry");
				boundary_d_t_sRadsymmetry(base,direction,z,dummy);
			}
		}
		else if(dset->Bc_Type==1)//ORTHOGONAL GRID Non-Homogeneous Neumann BC
		{
			if(vartag!= 16)
			{
				qo=*(dset->Data_Array);
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
			if((vartag==16)&&(base->BaseParameters->flagRadiation==1))
			{
				//printf("\nRadiation solution, insulation");
				boundary_d_t_sRadinsulation(base,direction,z,zBc,dummy);
			}
			if((vartag==16)&&(base->BaseParameters->flagflowRadiation==1))
			{
				//printf("\nRadiation solution, insulation");
				boundary_d_t_sInsulation(base,direction,z,zBc,dset,dummy);
			}
			
		}
		else if(dset->Bc_Type==3)//Homogeneous Dirichlet BC ... Put this if condition on TOP which is used often <<< Narsimha
		{
		  /*	if(base->flagLMMHD==1 && vartag != 4)
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
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							flagelectcond=z->Flowsol->flageleccond; 
							sigma = *(z->PrevFlowsol->Eleconductivity+indexP*flagelectcond);							
							p=*(surf+(indexP)*6*3+3*1+0);
							q=*(surf+(indexP)*6*3+3*1+1);
							r=*(surf+(indexP)*6*3+3*1+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k;index[1]=j;index[2]=i;
							  FindgradientPhiMHD(base,z,index,Dphi);
							
							nx=-p/Area;
							ny=-q/Area;
							nz=-r/Area;
							bx = *(z->Flowsol->MfieldX+indexP);
							by = *(z->Flowsol->MfieldY+indexP);
							bz = *(z->Flowsol->MfieldZ+indexP);
							bn = bx*nx+by*ny+bz*nz;
							vol = *(vols+(indexP));
							Dphi[0] = Dphi[0]/vol;
							Dphi[1] = Dphi[1]/vol;
							Dphi[2] = Dphi[2]/vol;
							dPx = -*(z->Flowsol->Pressure_Grad[0]+indexP)/vol;
							dPy = -*(z->Flowsol->Pressure_Grad[1]+indexP)/vol;
							dPz = -*(z->Flowsol->Pressure_Grad[2]+indexP)/vol;
							dPx = dPx - (dPx*nx+dPy*ny+dPz*nz)*nx;
							dPy = dPy - (dPx*nx+dPy*ny+dPz*nz)*ny;
							dPz = dPz - (dPx*nx+dPy*ny+dPz*nz)*nz;
							DbX = by*nz-bz*ny;
							DbY = bz*nx-bx*nz;
							DbZ = bx*ny-by*nx;
							EpX = sigma*bn*(Dphi[1]*nz-Dphi[2]*ny);
							EpY = sigma*bn*(Dphi[2]*nx-Dphi[0]*nz);
							EpZ = sigma*bn*(Dphi[0]*ny-Dphi[1]*nx);
							deltaP = *(zBc->deltaP+((nj-1)*k+j));
							kH = sqrt(sigma/diffCoef);
							kH = kH*bn;
						
							if(kH*deltaP > 0.01)
							{
							{if(vartag==1)
								qo= kH*(*(flowsol+indexP)-(dPx+EpX)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag==2)
									qo= kH*(*(flowsol+indexP)-(dPy+EpY)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag == 3)
								        qo= kH*(*(flowsol+indexP)-(dPz+EpZ)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
							}
							
						    												
							
							
							digonal=(*(beta+(indexP)*6*3+3*1+0))*diffCoef;
							if(digonal <= 0.0)
								digonal=-digonal;

							*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo*Area)/digonal;
						        }
							else 
							{
							*(dummy+((nj-1)*k+j)) = 0;
							}
						}
					break;
					
					case 2://west
						i=0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							flagelectcond=z->Flowsol->flageleccond; 
							sigma = *(z->PrevFlowsol->Eleconductivity+indexP*flagelectcond);		
							p=*(surf+(indexP)*6*3+3*0+0);
							q=*(surf+(indexP)*6*3+3*0+1);
							r=*(surf+(indexP)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
								index[0]=k;index[1]=j;index[2]=i;
							FindgradientPhiMHD(base,z,index,Dphi);
							nx=-p/Area;
							ny=-q/Area;
							nz=-r/Area;
							bx = *(z->Flowsol->MfieldX+indexP);
							by = *(z->Flowsol->MfieldY+indexP);
							bz = *(z->Flowsol->MfieldZ+indexP);
							bn = bx*nx+by*ny+bz*nz;
							vol = *(vols+(indexP));
							Dphi[0] = Dphi[0]/vol;
							Dphi[1] = Dphi[1]/vol;
							Dphi[2] = Dphi[2]/vol;
							dPx = -*(z->Flowsol->Pressure_Grad[0]+indexP)/vol;
							dPy = -*(z->Flowsol->Pressure_Grad[1]+indexP)/vol;
							dPz = -*(z->Flowsol->Pressure_Grad[2]+indexP)/vol;
							dPx = dPx - (dPx*nx+dPy*ny+dPz*nz)*nx;
							dPy = dPy - (dPx*nx+dPy*ny+dPz*nz)*ny;
							dPz = dPz - (dPx*nx+dPy*ny+dPz*nz)*nz;
							DbX = by*nz-bz*ny;
							DbY = bz*nx-bx*nz;
							DbZ = bx*ny-by*nx;
							EpX = sigma*bn*(Dphi[1]*nz-Dphi[2]*ny);
							EpY = sigma*bn*(Dphi[2]*nx-Dphi[0]*nz);
							EpZ = sigma*bn*(Dphi[0]*ny-Dphi[1]*nx);
							deltaP = *(zBc->deltaP+((nj-1)*k+j));
							kH = sqrt(sigma/diffCoef);
							kH = kH*bn;
						
							if(kH*deltaP > 0.01)
							{
							{if(vartag==1)
								qo= kH*(*(flowsol+indexP)-(dPx+EpX)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag==2)
									qo= kH*(*(flowsol+indexP)-(dPy+EpY)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag == 3)
								        qo= kH*(*(flowsol+indexP)-(dPz+EpZ)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
							}
							
							
							digonal=(*(beta+(indexP)*6*3+3*0+0))*diffCoef;
							
							if(digonal <= 0.0)
								digonal=-digonal;

							*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo*Area)/digonal;
						        }
							else 
							{
							*(dummy+((nj-1)*k+j)) = 0;
							}
						
						}
					break;
					
					case 3://north
						j=nj-2;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							flagelectcond=z->Flowsol->flageleccond; 
							sigma = *(z->PrevFlowsol->Eleconductivity+indexP*flagelectcond);		
							p=*(surf+(indexP)*6*3+3*3+0);
							q=*(surf+(indexP)*6*3+3*3+1);
							r=*(surf+(indexP)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k;index[1]=j;index[2]=i;
							FindgradientPhiMHD(base,z,index,Dphi);
							nx=-p/Area;
							ny=-q/Area;
							nz=-r/Area;
							bx = *(z->Flowsol->MfieldX+indexP);
							by = *(z->Flowsol->MfieldY+indexP);
							bz = *(z->Flowsol->MfieldZ+indexP);
							bn = bx*nx+by*ny+bz*nz;
							vol = *(vols+(indexP));
							Dphi[0] = Dphi[0]/vol;
							Dphi[1] = Dphi[1]/vol;
							Dphi[2] = Dphi[2]/vol;
							dPx = -*(z->Flowsol->Pressure_Grad[0]+indexP)/vol;
							dPy = -*(z->Flowsol->Pressure_Grad[1]+indexP)/vol;
							dPz = -*(z->Flowsol->Pressure_Grad[2]+indexP)/vol;
							dPx = dPx - (dPx*nx+dPy*ny+dPz*nz)*nx;
							dPy = dPy - (dPx*nx+dPy*ny+dPz*nz)*ny;
							dPz = dPz - (dPx*nx+dPy*ny+dPz*nz)*nz;
							DbX = by*nz-bz*ny;
							DbY = bz*nx-bx*nz;
							DbZ = bx*ny-by*nx;
							EpX = sigma*bn*(Dphi[1]*nz-Dphi[2]*ny);
							EpY = sigma*bn*(Dphi[2]*nx-Dphi[0]*nz);
							EpZ = sigma*bn*(Dphi[0]*ny-Dphi[1]*nx);
							deltaP = *(zBc->deltaP+((ni-1)*k+i));
							kH = sqrt(sigma/diffCoef);
							kH = kH*bn;
							if(kH*deltaP > 0.01)
							{
							{if(vartag==1)
								qo= kH*(*(flowsol+indexP)-(dPx+EpX)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag==2)
									qo= kH*(*(flowsol+indexP)-(dPy+EpY)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag == 3)
								        qo= kH*(*(flowsol+indexP)-(dPz+EpZ)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
							}
							
							digonal=*(beta+(indexP)*6*3+3*3+1)*diffCoef;
							
							if(digonal <= 0.0)
								digonal=-digonal;

							*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo*Area)/digonal;
						        }
							else 
							{
							*(dummy+((ni-1)*k+i)) = 0;
							}
													
						}
					break;
					
					case 4://south
						j=0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							flagelectcond=z->Flowsol->flageleccond; 
							sigma = *(z->PrevFlowsol->Eleconductivity+indexP*flagelectcond);		
							p=*(surf+(indexP)*6*3+3*2+0);
							q=*(surf+(indexP)*6*3+3*2+1);
							r=*(surf+(indexP)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k;index[1]=j;index[2]=i;
							FindgradientPhiMHD(base,z,index,Dphi);
							nx=-p/Area;
							ny=-q/Area;
							nz=-r/Area;
							bx = *(z->Flowsol->MfieldX+indexP);
							by = *(z->Flowsol->MfieldY+indexP);
							bz = *(z->Flowsol->MfieldZ+indexP);
							bn = bx*nx+by*ny+bz*nz;
							vol = *(vols+(indexP));
							Dphi[0] = Dphi[0]/vol;
							Dphi[1] = Dphi[1]/vol;
							Dphi[2] = Dphi[2]/vol;
							dPx = -*(z->Flowsol->Pressure_Grad[0]+indexP)/vol;
							dPy = -*(z->Flowsol->Pressure_Grad[1]+indexP)/vol;
							dPz = -*(z->Flowsol->Pressure_Grad[2]+indexP)/vol;
							dPx = dPx - (dPx*nx+dPy*ny+dPz*nz)*nx;
							dPy = dPy - (dPx*nx+dPy*ny+dPz*nz)*ny;
							dPz = dPz - (dPx*nx+dPy*ny+dPz*nz)*nz;
							DbX = by*nz-bz*ny;
							DbY = bz*nx-bx*nz;
							DbZ = bx*ny-by*nx;
							EpX = sigma*bn*(Dphi[1]*nz-Dphi[2]*ny);
							EpY = sigma*bn*(Dphi[2]*nx-Dphi[0]*nz);
							EpZ = sigma*bn*(Dphi[0]*ny-Dphi[1]*nx);
							deltaP = *(zBc->deltaP+((ni-1)*k+i));
							kH = sqrt(sigma/diffCoef);
							kH = kH*bn;
							if(kH*deltaP > 0.01)
							{
							{if(vartag==1)
								qo= kH*(*(flowsol+indexP)-(dPx+EpX)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag==2)
									qo= kH*(*(flowsol+indexP)-(dPy+EpY)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag == 3)
								        qo= kH*(*(flowsol+indexP)-(dPz+EpZ)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
							}
							
							digonal=(*(beta+(indexP)*6*3+3*2+1))*diffCoef;
							
							if(digonal <= 0.0)
								digonal=-digonal;

							*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo*Area)/digonal;
						        }
							else 
							{
							*(dummy+((ni-1)*k+i)) = 0;
							}
						}
					break;
					
					case 5://top
						k=nk-2;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							flagelectcond=z->Flowsol->flageleccond; 
							sigma = *(z->PrevFlowsol->Eleconductivity+indexP*flagelectcond);	
							p=*(surf+(indexP)*6*3+3*5+0);
							q=*(surf+(indexP)*6*3+3*5+1);
							r=*(surf+(indexP)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k;index[1]=j;index[2]=i;
							FindgradientPhiMHD(base,z,index,Dphi);
							nx=-p/Area;
							ny=-q/Area;
							nz=-r/Area;
							bx = *(z->Flowsol->MfieldX+indexP);
							by = *(z->Flowsol->MfieldY+indexP);
							bz = *(z->Flowsol->MfieldZ+indexP);
							bn = bx*nx+by*ny+bz*nz;
							vol = *(vols+(indexP));
							Dphi[0] = Dphi[0]/vol;
							Dphi[1] = Dphi[1]/vol;
							Dphi[2] = Dphi[2]/vol;
							dPx = -*(z->Flowsol->Pressure_Grad[0]+indexP)/vol;
							dPy = -*(z->Flowsol->Pressure_Grad[1]+indexP)/vol;
							dPz = -*(z->Flowsol->Pressure_Grad[2]+indexP)/vol;
							dPx = dPx - (dPx*nx+dPy*ny+dPz*nz)*nx;
							dPy = dPy - (dPx*nx+dPy*ny+dPz*nz)*ny;
							dPz = dPz - (dPx*nx+dPy*ny+dPz*nz)*nz;
							DbX = by*nz-bz*ny;
							DbY = bz*nx-bx*nz;
							DbZ = bx*ny-by*nx;
							EpX = sigma*bn*(Dphi[1]*nz-Dphi[2]*ny);
							EpY = sigma*bn*(Dphi[2]*nx-Dphi[0]*nz);
							EpZ = sigma*bn*(Dphi[0]*ny-Dphi[1]*nx);
							deltaP = *(zBc->deltaP+((ni-1)*j+i));
							kH = sqrt(sigma/diffCoef);
							kH = kH*bn;
							if(kH*deltaP > 0.01)
							{
							{if(vartag==1)
								qo= kH*(*(flowsol+indexP)-(dPx+EpX)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag==2)
									qo= kH*(*(flowsol+indexP)-(dPy+EpY)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag == 3)
								        qo= kH*(*(flowsol+indexP)-(dPz+EpZ)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
							}
							
							digonal=*(beta+(indexP)*6*3+3*5+2)*diffCoef;
							if(digonal <= 0.0)
								digonal=-digonal;

							*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo*Area)/digonal;
						        }
							else 
							{
							*(dummy+((ni-1)*j+i)) = 0;
							}
						}
					break;
						
					case 6://bottom
						k=0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
							flagelectcond=z->Flowsol->flageleccond; 
							sigma = *(z->PrevFlowsol->Eleconductivity+indexP*flagelectcond);	
							p=*(surf+(indexP)*6*3+3*4+0);
							q=*(surf+(indexP)*6*3+3*4+1);
							r=*(surf+(indexP)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							index[0]=k;index[1]=j;index[2]=i;
							FindgradientPhiMHD(base,z,index,Dphi);
							nx=-p/Area;
							ny=-q/Area;
							nz=-r/Area;
							bx = *(z->Flowsol->MfieldX+indexP);
							by = *(z->Flowsol->MfieldY+indexP);
							bz = *(z->Flowsol->MfieldZ+indexP);
							bn = bx*nx+by*ny+bz*nz;
							vol = *(vols+(indexP));
							Dphi[0] = Dphi[0]/vol;
							Dphi[1] = Dphi[1]/vol;
							Dphi[2] = Dphi[2]/vol;
							dPx = -*(z->Flowsol->Pressure_Grad[0]+indexP)/vol;
							dPy = -*(z->Flowsol->Pressure_Grad[1]+indexP)/vol;
							dPz = -*(z->Flowsol->Pressure_Grad[2]+indexP)/vol;
							dPx = dPx - (dPx*nx+dPy*ny+dPz*nz)*nx;
							dPy = dPy - (dPx*nx+dPy*ny+dPz*nz)*ny;
							dPz = dPz - (dPx*nx+dPy*ny+dPz*nz)*nz;
							DbX = by*nz-bz*ny;
							DbY = bz*nx-bx*nz;
							DbZ = bx*ny-by*nx;
							EpX = sigma*bn*(Dphi[1]*nz-Dphi[2]*ny);
							EpY = sigma*bn*(Dphi[2]*nx-Dphi[0]*nz);
							EpZ = sigma*bn*(Dphi[0]*ny-Dphi[1]*nx);
							deltaP = *(zBc->deltaP+((ni-1)*j+i));
							kH = sqrt(sigma/diffCoef);
							kH = kH*bn;
							if(kH*deltaP > 0.01)
							{
							{if(vartag==1)
								qo= kH*(*(flowsol+indexP)-(dPx+EpX)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag==2)
									qo= kH*(*(flowsol+indexP)-(dPy+EpY)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
								else if(vartag == 3)
								        qo= kH*(*(flowsol+indexP)-(dPz+EpZ)/(diffCoef*kH*kH)*(cosh(deltaP*kH)-1))/sinh(deltaP*kH);
							}
							digonal=(*(beta+(indexP)*6*3+3*4+2))*diffCoef;
							if(digonal <= 0.0)
								digonal=-digonal;

							*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo*Area)/digonal;
						        }
							else 
							{
							*(dummy+((ni-1)*j+i)) = 0;
							}
						}
					break;
				}
			 
			 
			}
             else{*/
                 for(size_ind=0;size_ind<size;size_ind++)
			{
				*(dummy+size_ind)=0.0;
			}
		 // }			
		}
		else if(dset->Bc_Type==2)//Dirichlet BC
		{
			if(vartag!=16)
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind)=*(dset->Data_Array);
				}
			}
			else
			{
				boundary_d_t_sRadiation(base,direction,z,zBc,dset,dummy);
			}
		}
		else if(dset->Bc_Type==41)//Robin BC for NON-ORTHOGONAL GRID
		{
			if(base->BaseParameters->flagRadiation==0 && vartag!=16)
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
						printf("Robin BC: variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
						exit(24);
					}
				}while(sum > convcri);
			}
			if(base->BaseParameters->flagflowRadiation==1 && vartag!=16)
				boundary_d_t_sRobin(base,direction,z,zBc,dset,dummy);
		}
		else if(dset->Bc_Type==40)//Robin BC for ORTHOGONAL GRID
		{
			if(base->BaseParameters->flagRadiation==0 && vartag!=16)
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
			if(base->BaseParameters->flagflowRadiation == 1 && vartag!=16)
				boundary_d_t_sRobin(base,direction,z,zBc,dset,dummy);
		}
		else if(dset->Bc_Type==5)//Orlanski BC
                 {
				Uav=*(dset->Data_Array);	
				Alpha=1.0;
				BetaCoeff=Uav*delta_t;
				
				MakeSlabPtrArray(z,slabVelocity,vartag);//Called for making slabs of u 
				MakeSlabPtrArray(z,slabVolume,5);//Called for making slabs of Volumes
			
				vols= z->ZGeom->Volume;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				switch(direction)
				{
					case 1://east

					        for(nobc=0;nobc<numbc;nobc++)
						{ 
						  	if(direction==(z->ZoneBc+nobc)->Direction)
							  slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
						
						}
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
                                                        Gamma=*(slabPreviousTimeStep+((nj-1)*k+j));
							Fdiff[1]=(*(beta+(indexP)*6*3+3*1+0))*(*(flowsol+indexP));
							res=Gamma*Area+BetaCoeff*Fdiff[1];//Diff;
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*1+0))+Alpha*Area;

							*(dummy+((nj-1)*k+j))=res/(digonal);
						}
					break;
					
					case 2://west
					        for(nobc=0;nobc<numbc;nobc++)
						{ 
						  	if(direction==(z->ZoneBc+nobc)->Direction)
							  slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
						
						}

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
							Gamma=*(slabPreviousTimeStep+((nj-1)*k+j));
							Fdiff[0]=(*(beta+(indexP)*6*3+3*0+0))*(*(flowsol+indexP));
							res=Gamma*Area-BetaCoeff*Fdiff[0];
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*0+0))+Alpha*Area;

							*(dummy+((nj-1)*k+j))=res/(digonal);
						}
							
					break;
					
					case 3://north

					        for(nobc=0;nobc<numbc;nobc++)
						{ 
						  	if(direction==(z->ZoneBc+nobc)->Direction)
							  slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
						
						}
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
							Gamma=*(slabPreviousTimeStep+((ni-1)*k+i));
							Fdiff[3]=(*(beta+(indexP)*6*3+3*3+1))*(*(flowsol+indexP));
							res=(Gamma*Area)+BetaCoeff*Fdiff[3];
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*3+1))+Alpha*Area;

							*(dummy+((ni-1)*k+i))=res/(digonal);
						}
					break;
					
					case 4://south	

					        for(nobc=0;nobc<numbc;nobc++)
						{ 
						  	if(direction==(z->ZoneBc+nobc)->Direction)
							  slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
						
						}
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
							Gamma=*(slabPreviousTimeStep+((ni-1)*k+i));
							Fdiff[2]=(*(beta+(indexP)*6*3+3*2+1))*(*(flowsol+indexP));
							res=Gamma*Area-BetaCoeff*Fdiff[2];
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*2+1))+Alpha*Area;

							*(dummy+((ni-1)*k+i))=res/(digonal);
						}
							
					break;
					
					case 5://top
					        for(nobc=0;nobc<numbc;nobc++)
						{ 
						  	if(direction==(z->ZoneBc+nobc)->Direction)
							  slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
						
						}
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
							Gamma=*(slabPreviousTimeStep+((ni-1)*j+i));
							Fdiff[5]=(*(beta+(indexP)*6*3+3*5+2))*(*(flowsol+indexP));
							res=(Gamma*Area)+BetaCoeff*Fdiff[5];
							digonal=BetaCoeff*(*(beta+(indexP)*6*3+3*5+2))+Alpha*Area;

							*(dummy+((ni-1)*j+i))=res/(digonal);
						}
							
					break;
					
					case 6://bottom

					        for(nobc=0;nobc<numbc;nobc++)
						{ 
						  	if(direction==(z->ZoneBc+nobc)->Direction)
							  slabPreviousTimeStep=FindSlabPtrBC(z->ZoneBc+nobc,vartag*10);
						
						}
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
							Gamma=*(slabPreviousTimeStep+((ni-1)*j+i));
							Fdiff[4]=(*(beta+(indexP)*6*3+3*4+2))*(*(flowsol+indexP));
							res=Gamma*Area-BetaCoeff*Fdiff[4];
							digonal=-1.00*BetaCoeff*(*(beta+(indexP)*6*3+3*4+2))+Alpha*Area;

							*(dummy+((ni-1)*j+i))=res/(digonal);
						}
					break;
				}
			 }

		
	}
	if(dset->identity[0]==1)// If Given Function Than it will Activate
	{
		if(dset->Bc_Type==2)//Drichlet BC
		{
			if(vartag!=16)
			{
				for(size_ind=0;size_ind<size;size_ind++)
				{
					*(dummy+size_ind)=*(dset->funct+size_ind);
				}
			}
			else
			{
				boundary_d_t_sRadiation(base,direction,z,zBc,dset,dummy);
			}
		}
		else if(dset->Bc_Type==10)//NON-ORTHOGONAL GRID Neumann BC (function) ... To be verified ... <<<<< Narsimha
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
							
							//Where are the following Alpha, BetaCoeff and Gamma variables put to use?? <<< Narsimha.
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
					printf("variableTag=%d Zone:%s Boundary Direction=%d residual has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,vartag,relaxation);
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
			surf=z->ZGeom->SurfaceArea;
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
					printf("Robin BC (function): variableTag=%d Zone:%s Boundary Direction=%d residual has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,vartag,relaxation);
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
			surf=z->ZGeom->SurfaceArea;

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

						*(dummy+((nj-1)*k+j))=res/(digonal);
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

						*(dummy+((ni-1)*k+i))=res/(digonal);
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
						*(dummy+((ni-1)*k+i))=res/(digonal);
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
						*(dummy+((ni-1)*j+i))=res/(digonal);
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
						*(dummy+((ni-1)*j+i))=res/(digonal);
					}
				break;
			}
		}

	}
        

	  
            if(dset->Bc_Type==9)//ORTHOGONAL GRID zero current efflux BC
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
                                                    
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+1);							
							diffCoef = 1.0;
							digonal=(*(beta+(indexP)*6*3+3*1+0))*diffCoef;
							
							*(dummy+((nj-1)*k+j))=*(flowsol+indexP)+(qo)/digonal;
						}
					break;
					
					case 2://west
                                               
						i=0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						  {                                        
						

							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
								qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+0);
							diffCoef = 1.0;													
							digonal=(*(beta+(indexP)*6*3+3*0+0))*diffCoef;							

							*(dummy+((nj-1)*k+j))=*(flowsol+indexP)-(qo)/digonal;
						}
					break;
					
					case 3://north
                                                
						j=nj-2;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						  {     
                                                    
							
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1.0;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+3);
							digonal=*(beta+(indexP)*6*3+3*3+1)*diffCoef;							

							*(dummy+((ni-1)*k+i))=*(flowsol+indexP)+(qo)/digonal;
						}
					break;
					
					case 4://south
                                                
						j=0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{       
                            
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1.0;
					     	qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+2);
							digonal=(*(beta+(indexP)*6*3+3*2+1))*diffCoef;							

							*(dummy+((ni-1)*k+i))=*(flowsol+indexP)-(qo)/digonal;
						}
					break;
					
					case 5://top
                                                
						k=nk-2;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						  {                                                            

							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1.0;
						qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+5);
							digonal=*(beta+(indexP)*6*3+3*5+2)*diffCoef;
								
							*(dummy+((ni-1)*j+i))=*(flowsol+indexP)+(qo)/digonal;
						}
					break;
						
					case 6://bottom
                                                
						k=0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{                                                      

							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1.0;
								qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+4);
							digonal=(*(beta+(indexP)*6*3+3*4+2))*diffCoef;
							

							*(dummy+((ni-1)*j+i))=*(flowsol+indexP)-(qo)/digonal;
						}
					break;
				}
			
			
			
		}


               if(dset->Bc_Type==90)//Non-ORTHOGONAL GRID zero current efflux BC
               {
			// Remember Inward flux is positive and outward flux is negative ... OK By Ashish
			qo  = *(dset->Data_Array);
			vols= z->ZGeom->Volume;
			beta= z->ZGeom->Beta;
		

			MakeSlabPtrArray(z,slabVelocity,vartag);/*Called for making slabs of u */
			MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/

			convcri = find_convcri(vartag,base);
			switch(direction)
			{
				case 1://east
					i=ni-2;
					do
					{
						slabIndex=0;
						sum=0.0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+1);
							index[0]=k,index[1]=j,index[2]=ni-2;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*1+0);
							
								res=(qo)-diffCoef*Fdiff[1];//Diff;
							
							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 2://west
					i=0;
					do
					{
						slabIndex=0;
						sum=0.0;
						for (k=0;k<nk-1;k++) 
						for (j=0;j<nj-1;j++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+0);
							index[0]=k,index[1]=j,index[2]=0;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*0+0));
						
								res=qo-diffCoef*Fdiff[0];
								
							*(dummy+((nj-1)*k+j))=*(dummy+((nj-1)*k+j)) - relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(nj-1)));
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 3://north
					j = nj-2;
					do
					{
						slabIndex=0;
						sum=0.0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+3);
							index[0]=k,index[1]=nj-2,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*3+1);
							
						
								res=(qo)-diffCoef*Fdiff[3];//Diff;
							
							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 4://south	
					j=0;
					do
					{
						slabIndex=0;
						sum=0.0;
						for (k=0;k<nk-1;k++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+2);
							index[0]=k,index[1]=0,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*2+1));						
						
								res=qo-diffCoef*Fdiff[2];
							
							*(dummy+((ni-1)*k+i))=*(dummy+((ni-1)*k+i)) - relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((nk-1)*(ni-1)));
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 5://top
					k=nk-2;
					do
					{
						slabIndex=0;
						sum=0.0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+5);
							index[0]=nk-2,index[1]=j,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=*(beta+(indexP)*6*3+3*5+2);
						
								res=(qo)-diffCoef*Fdiff[5];//Diff;
								
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) + relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
				
				case 6://bottom
					k=0;
					do
					{
						slabIndex=0;
						sum=0.0;
						for (j=0;j<nj-1;j++) 
						for (i=0;i<ni-1;i++)
						{
							indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
							diffCoef = 1;
							qo=*(z->Flowsol->PotentialMHDpoissonRHS+indexP*6+4);
							index[0]=0,index[1]=j,index[2]=i;
							Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
							digonal=(*(beta+(indexP)*6*3+3*4+2));
						
								res=qo-diffCoef*Fdiff[4];
								
							*(dummy+((ni-1)*j+i))=*(dummy+((ni-1)*j+i)) - relaxation * res/(diffCoef*digonal);
							sum=sum+res*res;
						}
						sum=sqrt(sum/((ni-1)*(nj-1)));
						
						if(sum > 1.0e8)
						{
							printf("variableTag=%d Zone:%s Boundary Direction=%d residual (%e) has diverged beyond 1.0e8 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", vartag,z->zonename,direction,sum,vartag,relaxation);
							exit(24);
						}
					}while(sum > convcri);
					break;
			}
		}

	return (sum);
}

double calcR(int i,int j,int ni,double *Ycor,double *R)
{
	double r1,r2,deltar,*Y,*Y1,*Y2;
	Y=Ycor+j*ni+i;
	Y1=Ycor+(j-1)*ni+i;
	Y2=Ycor+(j+1)*ni+i;
	r2=(*(Y)-*(Y1))/2.0;
	r1=(*(Y2)-*(Y))/2.0;
	deltar=r1-r2;
	*R=r1;
	return(deltar);
}

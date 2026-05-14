#include "APKupdate_boundaryTurbDMPFV.h"
/**
 This function  is used for deriving the boundary condition for phase momentum and global momentum  in case of turbulent flow because user has freedom to give velocity boundary condition which is converted into momentum by multiplying with density and volume fraction of corresponding phases and sum of all phase momentum is being done in order to get global momentum value at the boundary. This function is  useful only when boundary condtion given by user is dirichlet and homogeneous dirichlet. Only Homogeneous  Neumann boundary condition  has been considerd for multiphase flow but that is directly implemented through boundary_dtsTurb() function. 
 */
int update_boundaryTurbDMPFV(Base *base,Zone *CgnsZone,ZoneBoundary *bcptr,int vartag)
{

int flagrhophase=0,nophase=0,nop,k,j,i,dir,NI,NJ,NK,index,VARTAG;

int *ptrRange,*ptrSize,count,che;

int size;//added just to test the boundary values,other than this no use in subroutine-->chandrasekhar
	dir=bcptr->Direction;
	ptrRange = bcptr->Point_Range;
	nop=base->nphase;
	ptrSize = CgnsZone->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
// 	if(vartag==3300)
// 	{
// 	printf("The zone name is %s\n",CgnsZone->zonename);
// 	printf("\tThe boundary name is %s\n",bcptr->BCName);
// 	}
	if (dir==1||dir==2)
		size=(NK-1)*(NJ-1);
	if (dir==3||dir==4)
		size=(NK-1)*(NI-1);
	if (dir==5||dir==6)
		size=(NJ-1)*(NI-1);
	VARTAG=vartag;
	if(vartag%100==0)
		vartag=vartag/100;
	else
		vartag=vartag;
// 	printf("\nVARTAG--%d $$vartag--%d",VARTAG,vartag);
	count=find_nv(base,vartag);
	
	
	if((vartag==31)&&(VARTAG==3100))//called only for updation of mass momentum.(Vx_star) in mass velocity loop
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			if(bcptr->Boundary_Tag == 22)
			{	
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
					break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
					default:
					printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
					exit(1);
					break;
				}
			}
			else
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
								
							}
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumxstar+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}	
			}
		}
	}
	if((vartag==32)&&(VARTAG==3200))//called only for updation of mass momentum.(Vx_star) in mass velocity loop
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			if(bcptr->Boundary_Tag == 22)
			{	
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}
			}
			else
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
							for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumystar+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}	
			}
		}
	}
	if((vartag==33)&&(VARTAG==3300))//called only for updation of mass momentum.(Vx_star) in mass velocity loop
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			if(bcptr->Boundary_Tag == 22)
			{	
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
							for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
							for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;

					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}
			}
			else
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
							for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
		
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumzstar+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}	
			}
		}
	}
	if((vartag==31)&&(VARTAG==31))//called only for updation of momentum(Vx) in Calcvariable loop
	{ 
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			if(bcptr->Boundary_Tag == 22)
			{	
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)*(*((bcptr->ZSlab)->Rho+index));
						}
					
						break;
		
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;

					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)*(*((bcptr->ZSlab)->Rho+index));
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
	
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}
			}
			else
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
								
							}
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
								
							}
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							for(nophase=0;nophase<nop;nophase++)
							{
							*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							}
						}
							break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}	
			}
		}
	}
	if((vartag==32)&&(VARTAG==32))//called only for updation of momentum(Vx) in Calcvariable loop
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			if(bcptr->Boundary_Tag == 22)
			{	
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
				
						break;
		
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;

					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;

					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)*(*((bcptr->ZSlab)->Rho+index));		
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}
			}
			else
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							}	
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}	
			}
		}
	}
	if((vartag==33)&&(VARTAG==33))//called only for updation of momentum(Vx) in Calcvariable loop
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			if(bcptr->Boundary_Tag == 22)
			{	
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;
		
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
							*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;
	
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;

					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;

					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)*(*((bcptr->ZSlab)->Rho+index));	
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}
			}
			else
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
							for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
 							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
 							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{	
							index = k*(NI-1)+i;
							j=(NJ-2);
 							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
		
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
 							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}
						}
						break;
			
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							for(nophase=0;nophase<nop;nophase++)
							{
								*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
							}	
						}
						break;
					default:
						printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
						exit(1);
						break;
				}	
			}
		}
	}
	if((vartag==44)&&(VARTAG==44))
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			
			switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						*((bcptr->ZSlab)->Momentumofphasex[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[1]+index)));
					}
					
					break;
		
				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*((bcptr->ZSlab)->Momentumofphasex[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[1]+index)));
				
					}
				
					break;
		
				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{	
						index = k*(NI-1)+i;
						j=(NJ-2);
						*((bcptr->ZSlab)->Momentumofphasex[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[1]+index)));
					}
					break;
	
				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*((bcptr->ZSlab)->Momentumofphasex[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[1]+index)));
				
					}
					break;
		
				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*((bcptr->ZSlab)->Momentumofphasex[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[1]+index)));
					}
					break;
		
				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						*((bcptr->ZSlab)->Momentumofphasex[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[1]+index)));		
					}
		
					break;
				default:
					printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
					exit(1);
					break;
			}
		}
	}
	if((vartag==45)&&(VARTAG==45))
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			
			switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						*((bcptr->ZSlab)->Momentumofphasey[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[1]+index)));
					}
					
					break;
		
				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*((bcptr->ZSlab)->Momentumofphasey[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[1]+index)));
				
					}
				
					break;
		
				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{	
						index = k*(NI-1)+i;
						j=(NJ-2);
						*((bcptr->ZSlab)->Momentumofphasey[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[1]+index)));
					}
							break;
	
				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*((bcptr->ZSlab)->Momentumofphasey[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[1]+index)));
	
					}
					break;
		
				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*((bcptr->ZSlab)->Momentumofphasey[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[1]+index)));
					}
					break;
		
				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
					index = j*(NI-1)+i;
					k=0;
						*((bcptr->ZSlab)->Momentumofphasey[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[1]+index)));		
					}
		
					break;
				default:
					printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
					exit(1);
					break;
			}
		}
	}
	if((vartag==46)&&(VARTAG==46))
	{
		if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
		{
			
			
			switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						*((bcptr->ZSlab)->Momentumofphasez[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[1]+index)));
					}
					
					break;
		
				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*((bcptr->ZSlab)->Momentumofphasez[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[1]+index)));
			
					}
			
					break;
		
				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{	
						index = k*(NI-1)+i;
						j=(NJ-2);
						*((bcptr->ZSlab)->Momentumofphasez[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[1]+index)));
					}
					break;
	
				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*((bcptr->ZSlab)->Momentumofphasez[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[1]+index)));

					}
					break;
		
				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*((bcptr->ZSlab)->Momentumofphasez[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[1]+index)));
					}
					break;
		
				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						*((bcptr->ZSlab)->Momentumofphasez[1]+index)=(*((bcptr->ZSlab)->VarthetaPhase[1]+index)*(*((bcptr->ZSlab)->RhoofPhase[1]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[1]+index)));		
					}

					break;
				default:
					printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
					exit(1);
					break;
			}
		}
	}
// 	if((vartag==31)&&(VARTAG==3100))
// 	{
// 	for(che=0;che<1;che++)
// 	printf("The value of the converted global momentumX is %lf\n",*((bcptr->ZSlab)->Momentumxstar+che));
// 	}
// 	if((vartag==32)&&(VARTAG==3200))
// 	{
// 	for(che=0;che<1;che++)
// 	printf("The value of the converted global momentumY is %lf\n",*((bcptr->ZSlab)->Momentumystar+che));
// 	}
// 	if((vartag==33)&&(VARTAG==3300))
// 	{
// 		for(che=0;che<1;che++)
// 			printf("The value of the converted global momentumYZ is %lf\n",*((bcptr->ZSlab)->Momentumzstar+che));
// 	}
// 	if(((vartag==31)&&(VARTAG==31)))
// 	{
// 	for(che=0;che<1;che++)
// 	printf("The value of the converted global momentumX is %lf\n",*((bcptr->ZSlab)->Momentumx+che));
// 	}
// 	if(((vartag==32)&&(VARTAG==32)))
// 	{
// 		for(che=0;che<1;che++)
// 			printf("The value of the converted global momentumY is %lf\n",*((bcptr->ZSlab)->Momentumy+che));
// 	}
// 	if(((vartag==33)&&(VARTAG==33)))
// 	{
// 		for(che=0;che<1;che++)
// 			printf("The value of the converted global momentumZ is %lf\n",*((bcptr->ZSlab)->Momentumz+che));
// 	}
// 	if((vartag==44)&&(VARTAG==44))
// 	{
// 		for(che=0;che<1;che++)
// 			printf("The momentumX of phase 1 is %lf\n",*((bcptr->ZSlab)->Momentumofphasex[1]+che));
// 	
// 	}
// 	if((vartag==45)&&(VARTAG==45))
// 	{
// 		for(che=0;che<1;che++)
// 			printf("The momentumY of phase 1 is %lf\n",*((bcptr->ZSlab)->Momentumofphasey[1]+che));
// 	
// 	}
// 	if((vartag==46)&&(VARTAG==46))
// 	{
// 		for(che=0;che<1;che++)
// 			printf("The momentumY of phase 1 is %lf\n",*((bcptr->ZSlab)->Momentumofphasez[1]+che));
// 	
// 	}
	return(0);
	
}	

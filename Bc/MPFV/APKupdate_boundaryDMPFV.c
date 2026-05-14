#include "APKupdate_boundaryDMPFV.h"
/**
 This function  is used for deriving the boundary condition for phase momentum and global momentum because user has freedom to give velocity boundary condition which is converted into momentum by multiplying with density and volume fraction of corresponding phases and sum of all phase momentum is being done in order to get global momentum value at the boundary. This function is  useful only when boundary condtion given by user is dirichlet and homogeneous dirichlet. Only Homogeneous  Neumann boundary condition  has been considerd for multiphase flow but that is directly implemented through boundary_dts() function. 
 */
int update_boundaryDMPFV(Base *base,Zone *CgnsZone,ZoneBoundary *bcptr)
{

int flagrhophase=0,nophase=0,nop,k,j,i,dir,NI,NJ,NK,index,che;

int *ptrRange,*ptrSize,count;
	dir=bcptr->Direction;
	ptrRange = bcptr->Point_Range;
	nop=base->nphase;
	ptrSize = CgnsZone->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	count=find_nv(base,41);
	
//  	printf("The zone name is %s\n",CgnsZone->zonename);
// 	printf("The boundary name is %s\n\n",bcptr->BCName);
	
	
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
					
					
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						
	
						if(flagrhophase==1)
						{	
							
						
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
							
						
						}
						else
						{	
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
					
						}
							
						
					}
					
/************************overall momentum calculation from individual momentum***************************************/ 		
					*((bcptr->ZSlab)->Momentumx+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						
					
					}
					*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumx+index);
					
					
				}
				
			break;
	
			case 2: /*Dir 2 ->Imin,West*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						
						
						if(flagrhophase==1)
						{
							
						
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
							
							
						
						}
					
						else
						{	
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
						
							
						}	
					
					}
					
/************************overall momentum calculation from individual momentum***************************************/ 	
					*((bcptr->ZSlab)->Momentumx+index)=0.0;
					
						
					for(nophase=0;nophase<nop;nophase++)
					{
						*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						
						
					}
					*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumx+index);
				
				}
				
			break;
	
			case 3:/*Dir 3 ->Jmax,North*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{	
					index = k*(NI-1)+i;
					j=(NJ-2);
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
							
							
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
							
							
	
								
								
							
							
						}
						else
						{	
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
							
							
						}	
					
						
					}		
/************************overall momentum calculation from individual momentum***************************************/ 	
					*((bcptr->ZSlab)->Momentumx+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						
					}
					*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumx+index);
				
				}
			break;
	
			case 4:/*Dir 4 ->Jmin,South*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = k*(NI-1)+i;
					j=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
	// 						*((bcptr->ZSlab)->Momentumofphasex[nophase])=0;
							
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
	
	// 						
						}	
						else
						{	
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
						
	 					}	
					
						
					}
					
						
/************************overall momentum calculation from individual momentum***************************************/ 		
					*((bcptr->ZSlab)->Momentumx+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							
					}
					*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumx+index);
				
				}
			break;
	
			case 5:/*Dir 5 ->Kmax,Top*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=(NK-2);
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
	// 						
								
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
								
	// 						
						}
					
						else
						{	
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
						
						}	
					}	
/************************overall momentum calculation from individual momentum***************************************/ 
					*((bcptr->ZSlab)->Momentumx+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						
						
						}
						
						*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumx+index);
				}
			break;
	
			case 6:/*Dir 6 ->Kmin,Bottom*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
							
							
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
							
							
						}	
					
						else
						{	
							*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index)));
						
						
						}	
					}		
/************************overall momentum calculation from individual momentum***************************************/ 
					*((bcptr->ZSlab)->Momentumx+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						*((bcptr->ZSlab)->Momentumx+index)=*((bcptr->ZSlab)->Momentumx+index)+(*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							
							
						
					}
				
					*((bcptr->ZSlab)->Momentumxstar+index)=*((bcptr->ZSlab)->Momentumx+index);		
				}
				
			break;
			default:
			printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
			exit(1);
			break;
		}
	}
	
	count=find_nv(base,42);
	
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
					
					
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						
	// 					
						if(flagrhophase==1)
						{	
// 						printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
						
							
					
							
						}
						else
						{	
					
							
					
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
						
						}	
						
					}
					
/************************overall momentum calculation from individual momentum***************************************/ 		
					
					*((bcptr->ZSlab)->Momentumy+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
												
						*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							
							
						
					}
					*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumy+index);
				}
				
			break;
	
			case 2: /*Dir 2 ->Imin,West*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						
						
						if(flagrhophase==1)
						{
						
							
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
						
						}
					
						else
						{	
						
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
							
							
						}	
					
					}
					
/************************overall momentum calculation from individual momentum***************************************/ 	
					
					*((bcptr->ZSlab)->Momentumy+index)=0.0;
					
						
					for(nophase=0;nophase<nop;nophase++)
					{
												
						*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							
							
						
					}
					*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumy+index);
				}
				
			break;
	
			case 3:/*Dir 3 ->Jmax,North*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{	
					index = k*(NI-1)+i;
					j=(NJ-2);
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
						
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
						
						}
						else
						{	
							
						
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
							
							
						}	
					
						
					}		
/************************overall momentum calculation from individual momentum***************************************/ 	
					
					*((bcptr->ZSlab)->Momentumy+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
													
						*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							
							
						
					}
					*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumy+index);
				}
			break;
	
			case 4:/*Dir 4 ->Jmin,South*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = k*(NI-1)+i;
					j=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
	// 						
	// 						printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
							
	// 						*((bcptr->ZSlab)->Momentumofphasey[nophase])=0;	
							
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
							
						}	
						else
						{	
							
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
							
						
						}	
					
						
					}
					
						
/************************overall momentum calculation from individual momentum***************************************/ 		
					
					*((bcptr->ZSlab)->Momentumy+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						
							
						*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							
							
						
					}
				
					*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumy+index);
				}
			break;
	
			case 5:/*Dir 5 ->Kmax,Top*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=(NK-2);
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
							
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
						
						}
					
						else
						{	
						
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
						
						}	
					}	
/************************overall momentum calculation from individual momentum***************************************/ 
					
					*((bcptr->ZSlab)->Momentumy+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						
							
						*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
							
						
					}
						
					*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumy+index);
				}
			break;
	
			case 6:/*Dir 6 ->Kmin,Bottom*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
							
							
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
								
	 					}	
					
						else
						{	
							
						
							*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index)));
							
							
						}	
					}		
/************************overall momentum calculation from individual momentum***************************************/ 
					
					*((bcptr->ZSlab)->Momentumy+index)=0.0;
					
					for(nophase=0;nophase<nop;nophase++)
					{
						
						
						*((bcptr->ZSlab)->Momentumy+index)=*((bcptr->ZSlab)->Momentumy+index)+(*((bcptr->ZSlab)->Momentumofphasey[nophase]+index));
					}	
					
					*((bcptr->ZSlab)->Momentumystar+index)=*((bcptr->ZSlab)->Momentumy+index);
					
				}
				
			break;
			default:
			printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
			exit(1);
			break;
		}
	}

	count=find_nv(base,43);
	
	if(((bcptr->Data_Set+(count))->Bc_Type==2)||((bcptr->Data_Set+(count))->Bc_Type==3))
	{
	
		//printf("dir=%d\n",dir);
		switch(dir)
		{
			case 1: /*Dir 1 ->Imax,East*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=(NI-2);
					
					
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						
	 					
						if(flagrhophase==1)
						{	
							
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}
						else
						{	
							
							
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
						
					}
					
/************************overall momentum calculation from individual momentum***************************************/ 		
					
					
					*((bcptr->ZSlab)->Momentumz+index)=0.0;
					for(nophase=0;nophase<nop;nophase++)
					{
						
													
						*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
					}
					*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumz+index);
				}
				
			break;
	
			case 2: /*Dir 2 ->Imin,West*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						
						
						if(flagrhophase==1)
						{
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}
					
						else
						{	
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
					
					}
					
/************************overall momentum calculation from individual momentum***************************************/ 	
					
					
					*((bcptr->ZSlab)->Momentumz+index)=0.0;
						
					for(nophase=0;nophase<nop;nophase++)
					{
						
							
						*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
					}
					*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumz+index);
				}
				
			break;
	
			case 3:/*Dir 3 ->Jmax,North*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{	//printf("3hello");
					index = k*(NI-1)+i;
					j=(NJ-2);
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
	 						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}
						else
						{	
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
					
						
					}		
/************************overall momentum calculation from individual momentum***************************************/ 	
					
					*((bcptr->ZSlab)->Momentumz+index)=0.0;
					for(nophase=0;nophase<nop;nophase++)
					{
						
						*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
					}
					*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumz+index);
				}
			break;
	
			case 4:/*Dir 4 ->Jmin,South*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = k*(NI-1)+i;
					j=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
						else
						{
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
					
						
					}
					
						
/************************overall momentum calculation from individual momentum***************************************/ 		
					
					*((bcptr->ZSlab)->Momentumz+index)=0.0;
					for(nophase=0;nophase<nop;nophase++)
					{
						
						*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
					}
				
					*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumz+index);
				}
			break;
	
			case 5:/*Dir 5 ->Kmax,Top*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=(NK-2);
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}
					
						else
						{
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						
						}	
					}	
/************************overall momentum calculation from individual momentum***************************************/ 
					
					*((bcptr->ZSlab)->Momentumz+index)=0.0;
					for(nophase=0;nophase<nop;nophase++)
					{	
						*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
					}
						
					*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumz+index);
				}
			break;
	
			case 6:/*Dir 6 ->Kmin,Bottom*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=0;
					for(nophase=0;nophase<nop;nophase++)
					{	flagrhophase=0;
						flagrhophase=(base->BaseParameters->flagRhoPhase[nophase]);
						if(flagrhophase==1)
						{
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((bcptr->ZSlab)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
					
						else
						{	
						
							*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)=(*((bcptr->ZSlab)->VarthetaPhase[nophase]+index)*(*((CgnsZone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index)));
						}	
					}		
/************************overall momentum calculation from individual momentum***************************************/ 
					
					*((bcptr->ZSlab)->Momentumz+index)=0.0;
					for(nophase=0;nophase<nop;nophase++)
					{
						
						*((bcptr->ZSlab)->Momentumz+index)=*((bcptr->ZSlab)->Momentumz+index)+(*((bcptr->ZSlab)->Momentumofphasez[nophase]+index));
					}
				
					*((bcptr->ZSlab)->Momentumzstar+index)=*((bcptr->ZSlab)->Momentumz+index);		
				}
				
			break;
			default:
			printf("Case not found,Terminating the Program from update_boundaryDMPFV\n");
			exit(1);
			break;
		}
	}
/*	for(che=0;che<1;che++)
		printf("The value of the converted momentumX for continuous phase  is %lf\n",*((bcptr->ZSlab)->Momentumofphasex[0]+che));
	for(che=0;che<1;che++)
		printf("The value of the converted momentumX for dispersed phase is %lf\n",*((bcptr->ZSlab)->Momentumofphasex[1]+che));
	for(che=0;che<1;che++)
		printf("The value of the converted global momentumX is %lf\n",*((bcptr->ZSlab)->Momentumx+che));
	for(che=0;che<1;che++)
		printf("The value of the converted momentumY for continuous phase  is %lf\n",*((bcptr->ZSlab)->Momentumofphasey[0]+che));
	for(che=0;che<1;che++)
		printf("The value of the converted momentumY for dispersed phase is %lf\n",*((bcptr->ZSlab)->Momentumofphasey[1]+che));
	for(che=0;che<1;che++)
		printf("The value of the converted global momentumY is %lf\n",*((bcptr->ZSlab)->Momentumy+che));
	for(che=0;che<1;che++)
		printf("The value of the converted momentumZ for continuous phase  is %lf\n",*((bcptr->ZSlab)->Momentumofphasez[0]+che));
	for(che=0;che<1;che++)
		printf("The value of the converted momentumZ for dispersed phase is %lf\n",*((bcptr->ZSlab)->Momentumofphasez[1]+che));
	for(che=0;che<1;che++)
		printf("The value of the converted global momentumZ is %lf\n",*((bcptr->ZSlab)->Momentumz+che));*/	
return(0);
	
}	

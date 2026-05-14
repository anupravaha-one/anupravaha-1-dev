#include "APKupdate_boundaryDisperse.h"

int update_boundaryDisperse(Base *base,Zone *CgnsZone,ZoneBoundary *bcptr)
{

int flagrhoDisperse=0,k,j,i,dir,NI,NJ,NK,index;

int *ptrRange,*ptrSize,count;
	dir=bcptr->Direction;
	ptrRange = bcptr->Point_Range;
	ptrSize = CgnsZone->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
//---------------------------------------------------
int num_zone;
//int a,b,c,numiface,numBc,iface,zBc;
Zone *z;
//------------------------------------------------------
z=CgnsZone; /* z is the pointer to first Zone*/
num_zone=base->Number_of_zones;

	//printf("\n ZoneName-%s",CgnsZone->zonename);
	//printf("\n Boundary-%s",bcptr->BCName);
//----------------------------------------------	
	count=find_nv(base,62);
	if((bcptr->Data_Set+(count))->Bc_Type==2)
	{
		switch(dir)
		{
			case 1: /*Dir 1 ->Imax,East*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=(NI-2);
					
					
						flagrhoDisperse=0;
						flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	
						if(flagrhoDisperse==1)
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}
						else
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}	
				}
			break;
			case 2: /*Dir 2 ->Imin,West*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
					       {	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
					 	}
					else
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}	
				}
			break;
			case 3:/*Dir 3 ->Jmax,North*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{	
					index = k*(NI-1)+i;
					j=(NJ-2);
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}
					else
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}		
				}		
			break;
	
			case 4:/*Dir 4 ->Jmin,South*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = k*(NI-1)+i;
					j=0;
						flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}
					else
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}		
				}
			break;
	
			case 5:/*Dir 5 ->Kmax,Top*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=(NK-2);
						flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}
					else
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}		
				}
			break;
	
			case 6:/*Dir 6 ->Kmin,Bottom*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=0;
						flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}
					else
						{	
							*((bcptr->ZSlab)->MomentumDispersex+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersex+index)));
							*((bcptr->ZSlab)->MomentumDispersexOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersexOld+index)));
						}		
				}
				
			break;
			default:
			printf("Case not found,Terminating the Program from update_boundaryDisperse\n");
			exit(1);
			break;
		}
	}
	//printf("\nMomentum-X=%lf",*((bcptr->ZSlab)->MomentumDispersex));
	count=find_nv(base,63);
	
	if((bcptr->Data_Set+(count))->Bc_Type==2)
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
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
						//printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}	
						
				}
			break;
	
			case 2: /*Dir 2 ->Imin,West*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
						//printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
				}
				
			break;
	
			case 3:/*Dir 3 ->Jmax,North*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{	
					index = k*(NI-1)+i;
					j=(NJ-2);
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
						//printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
				}
			break;
	
			case 4:/*Dir 4 ->Jmin,South*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = k*(NI-1)+i;
					j=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
						//printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
				}
			break;
	
			case 5:/*Dir 5 ->Kmax,Top*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=(NK-2);
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
						//printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
				}
			break;
	
			case 6:/*Dir 6 ->Kmin,Bottom*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
					if(flagrhoDisperse==1)
						{	
						//printf("Value updated = %e\n",*((bcptr->ZSlab)->Momentumofphasex[nophase]+index));
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersey+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersey+index)));
						*((bcptr->ZSlab)->MomentumDisperseyOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDisperseyOld+index)));
						}
				}
				
			break;
			default:
			printf("Case not found,Terminating the Program from update_boundaryDisperse\n");
			exit(1);
			break;
		}
	}
	//printf("\nMomentum-Y=%lf",*((bcptr->ZSlab)->MomentumDispersey));
	if(base->BaseParameters->flag2d == 0)
	{
	count=find_nv(base,64);
	
	if((bcptr->Data_Set+(count))->Bc_Type==2)
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
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	 					
						if(flagrhoDisperse==1)
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}	
						
				}
			break;
	
			case 2: /*Dir 2 ->Imin,West*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				{
					index = k*(NJ-1)+j;
					i=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	 					
						if(flagrhoDisperse==1)
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}	
				}	
			break;
	
			case 3:/*Dir 3 ->Jmax,North*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{	//printf("3hello");
					index = k*(NI-1)+i;
					j=(NJ-2);
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	 					
						if(flagrhoDisperse==1)
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}	
					
				}
			break;
	
			case 4:/*Dir 4 ->Jmin,South*/
				for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = k*(NI-1)+i;
					j=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	 					
						if(flagrhoDisperse==1)
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}	
					
				}
			break;
	
			case 5:/*Dir 5 ->Kmax,Top*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=(NK-2);
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	 					
						if(flagrhoDisperse==1)
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}	
					
				}
			break;
	
			case 6:/*Dir 6 ->Kmin,Bottom*/
				for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
				for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
				{
					index = j*(NI-1)+i;
					k=0;
					flagrhoDisperse=0;
					flagrhoDisperse=(base->BaseParameters->flagRhoDisperse);
						
	 					
						if(flagrhoDisperse==1)
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((bcptr->ZSlab)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}
						else
						{	
						*((bcptr->ZSlab)->MomentumDispersez+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersez+index)));
						*((bcptr->ZSlab)->MomentumDispersezOld+index)=(*((bcptr->ZSlab)->VarthetaDisperse+index)*(*((CgnsZone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((bcptr->ZSlab)->MomentumDispersezOld+index)));
						}	
							
				}
				
			break;
			default:
			printf("Case not found,Terminating the Program from update_boundaryDisperse\n");
			exit(1);
			break;
		}
	}
	}
//printf("_+_+_+_+_+_+_+++++++++_+_+_+_____+_+_+_+_+_+_+_+_+_+_+_+_+_+_\n\n\n\n");
// for(j=0;j<num_zone;j++)
// 	{
// 		numiface=(z+j)->Number_of_Interfaces;
// 		numBc=(z+j)->Number_of_BC_Patches;
// 		ptrSize=(z+j)->Size;
// 
// 		NI=*(ptrSize+0);
// 		NJ=*(ptrSize+1);
// 		NK=*(ptrSize+2);
// 		//printf("\nNi %d NJ %d Nk %d\n",NI,NJ,NK);
// 		for(dir=1;dir<=6;dir++)
// 		{
// 			for(iface=0;iface<numiface;iface++)
// 			{
// 				if(((z+j)->Iface+iface)->Direction==dir)
// 				{
// 					for(c=0;c<NK-1;c++)
// 					for(b=0;b<NJ-1;b++)
// 					for(a=0;a<NI-1;a++)
// 					{
// 						if(dir==1||dir==2)
// 							index= c*(NJ-1)+b;
// 						if(dir==3||dir==4)
// 							index= c*(NI-1)+a;
// 						if(dir==5||dir==6)
// 							index= b*(NI-1)+a;
// 						//printf("\n Udisperse %lf",*((((zone)->Iface+iface)->ZSlab)->slabUDisperse+index));  
// 						//*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);
// 					
// 					}
// 				
// 				}
// 			}
// 			for(zBc=0;zBc<numBc;zBc++)
// 			{
// 				if(((z+j)->ZoneBc+zBc)->Direction==dir)
// 				{
// 				//printf("\n----------------------------------------\n");
// 				//printf("\n----------------------------------------\n");
// 				//printf("\nNi %d NJ %d Nk %d\n",NI,NJ,NK);	
// 					for(c=0;c<NK-1;c++)
// 					for(b=0;b<NJ-1;b++)
// 					for(a=0;a<NI-1;a++)
// 					{
// 						if(dir==1||dir==2)
// 							index= c*(NJ-1)+b;
// 						if(dir==3||dir==4)
// 							index= c*(NI-1)+a;
// 						if(dir==5||dir==6)
// 							index= b*(NI-1)+a;
// 						//printf("\n %d %lf\t %lf\t %lf\t %lf\t  \n",dir,*(((z+j)->ZoneBc+zBc)->ZSlab->MomentumDispersex+index),*(((z+j)->ZoneBc+zBc)->ZSlab->MomentumDispersey+index),*(((z+j)->ZoneBc+zBc)->ZSlab->RhoDisperse+index*(base->BaseParameters->flagRhoDisperse)),*(((z+j)->ZoneBc+zBc)->ZSlab->VarthetaDisperse+index));
// 					//printf("\n###################################\n");
// 					}
// 				}
// 			}
// 		}
// 		
// 	}	
return(0);
	
}	

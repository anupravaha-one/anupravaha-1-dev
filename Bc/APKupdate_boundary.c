#include "APKupdate_boundary.h"
/**
This module inplements the boundary conditions according to the the BC type entered by the user
for the particular variable at any boundary patch.
*/
int update_boundary(Base *base,int nv,Zone *z,ZoneBoundary *ZoneBc)
{	
	int bc_type,che;
	
	if(nv!=5)
	{
		int ni,nj,nk,direct_index,size,vartag;
		double *domain,*slab;
		vartag = nv;

		ni=*(z->Size+0);
		nj=*(z->Size+1);
		nk=*(z->Size+2); 
		
		nv = find_nv(base,nv);
 	
		direct_index=ZoneBc->Direction;
		ZoneBc->NueRms=0.0;
 	
		
		bc_type=((ZoneBc)->Data_Set+nv)->Bc_Type;

		if (direct_index==1||direct_index==2)
			size=(nk-1)*(nj-1);
		if (direct_index==3||direct_index==4)
			size=(nk-1)*(ni-1);
		if (direct_index==5||direct_index==6)
			size=(nj-1)*(ni-1);
		
		domain = MakeDomainArray(z,vartag);
		
		slab   = FindSlabPtrBC(ZoneBc,vartag);
		
		
		if((base->flagTurbModel == 1 || base->flagTurbModel == 3) && z->flagSolidFluid == 0)
		{
			if(!(((vartag==62)&&(bc_type==2||bc_type==3))||((vartag==63)&&(bc_type==2||bc_type==3))||((vartag==64)&&(bc_type==2||bc_type==3))))
			{
				//printf("\nI am entering Vartag=%d$$Bc_Type=%d$zone-%s$Bc-%s",vartag,bc_type,z->zonename,ZoneBc->BCName);
				ZoneBc->NueRms=boundary_d_t_sTurb(base,ZoneBc,direct_index,domain,vartag,z,size,(ZoneBc->Data_Set+nv),slab);
			}
		}
		else if((base->flagMPFV == 0)&&(base->flagParticle == 0)&&(base->flagDroplet == 0))
		{
			ZoneBc->NueRms=boundary_d_t_s(base,direct_index,domain,vartag,z,size,(ZoneBc->Data_Set+nv),slab);
// 			if(vartag == 4)
// 				printf("Zone:%s BCNueRms=%e vartag=%d, Direction=%d\n",z->zonename,ZoneBc->NueRms,vartag,direct_index);
		}
		else if(mulphase(vartag,bc_type)==0)
		{
			ZoneBc->NueRms=boundary_d_t_s(base,direct_index,domain,vartag,z,size,(ZoneBc->Data_Set+nv),slab);
		}
		//added for checking purpose-->chandrasekhar
// 		for(che=0;che<1;che++)
//  		{
// 			printf("\nBoundary name is %s\n",ZoneBc->BCName);
// 			if(vartag==4)
//  				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->P+che));
// 			else if(vartag==1)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->u+che));
// 			else if(vartag==2)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->v+che));
// 			else if(vartag==65)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->VarthetaDisperse+che));
// 			else if(vartag==62)
// 				printf("\nThe value of %d at the boundary %s is %lf\n",vartag,ZoneBc->BCName,*(ZoneBc->ZSlab->MomentumDispersex+che));
// 			else if(vartag==63)
// 				printf("\nThe value of %d at the boundary %s is %lf\n",vartag,ZoneBc->BCName,*(ZoneBc->ZSlab->MomentumDispersey+che));
// 			else if(vartag==31)
// 				printf("\nThe value of %d at the boundary is %lf and masstag=%d\n",vartag,*(ZoneBc->ZSlab->Momentumx+che),base->masstag);
// 			else if(vartag==32)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumy+che));
// 			else if(vartag==33)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumz+che));
// 			else if(vartag==34)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->VarthetaPhase[0]+che));
// 			else if(vartag==35)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->VarthetaPhase[1]+che));
// 			else if(vartag==41)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumofphasex[0]+che));
// 			else if(vartag==42)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumofphasey[0]+che));
// 			else if(vartag==43)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumofphasez[0]+che));
// 			else if(vartag==44)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumofphasex[1]+che));
// 			else if(vartag==45)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumofphasey[1]+che));
// 			else if(vartag==46)
// 				printf("\nThe value of %d at the boundary is %lf\n",vartag,*(ZoneBc->ZSlab->Momentumofphasez[1]+che));
//  			else if(vartag==600)
// 				printf("\nThe value of %d at the boundary %s is %lf\n",vartag,ZoneBc->BCName,*(ZoneBc->ZSlab->Scal[0]+che));
// 			else if(vartag==601)
// 				printf("\nThe value of %d at the boundary %s is %lf\n",vartag,ZoneBc->BCName,*(ZoneBc->ZSlab->Scal[1]+che));
//  		}
		//added for checking purpose ends-->chandrasekhar
	}
	return 0;
}

int mulphase(int vartag,int type)
{
	int c=0;
	if(vartag>=31 && vartag<=33 && (type==2||type==3))
		c=1;
	if(vartag>=41 && vartag<=55 && (type==2||type==3))
	c=1;
//The following lines are commented as the disperse-phase momentum equations are being solved for velocity components and not momentum components
/*	if(vartag>=62 && vartag<=64 && (type==2||type==3))//adhiraj
		c=1;
*/	
	return (c);
} //by sushil	

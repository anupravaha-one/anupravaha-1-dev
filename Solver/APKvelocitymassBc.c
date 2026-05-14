#include "APKvelocitymassBc.h"
/**
The function velocitymassBc() impliments the boundary conditions for massvelocities. It applies the same boundary
conditions entered by the user for velocity components, to their respective mass velocity components. 
 */

int velocitymassBc(Base *base,int vartag,Zone *z,ZoneBoundary *ZoneBc)
{
	int ni,nj,nk,direct_index,size,nv;
	double *domain,*slab,bc_type;
	
	if(nv!=5)
	{
		ni=*(z->Size+0);
		nj=*(z->Size+1);
		nk=*(z->Size+2); 
		
		nv = find_nv(base,vartag);
		
		direct_index=ZoneBc->Direction;
		ZoneBc->NueRms=0.0;
 	
		
		bc_type=((ZoneBc)->Data_Set+nv)->Bc_Type;
		
		if (direct_index==1||direct_index==2)
			size=(nk-1)*(nj-1);
		if (direct_index==3||direct_index==4)
			size=(nk-1)*(ni-1);
		if (direct_index==5||direct_index==6)
			size=(nj-1)*(ni-1);
		
		domain = MakeDomainArray(z,vartag*100);
		
		slab   = FindSlabPtrBC(ZoneBc,vartag*100);
		
		if((base->flagTurbModel == 1 || base->flagTurbModel == 3) && z->flagSolidFluid == 0)
		{
			ZoneBc->NueRms=boundary_d_t_sTurb(base,ZoneBc,direct_index,domain,vartag*100,z,size,(ZoneBc->Data_Set+nv),slab);
		}
		else if(base->flagMPFV==0)
		{
			ZoneBc->NueRms=boundary_d_t_s(base,direct_index,domain,vartag*100,z,size,(ZoneBc->Data_Set+nv),slab);
		}
		else if(base->flagMPFV==1)
		{
			if(mulphase(vartag,bc_type)==0)
			{
				ZoneBc->NueRms=boundary_d_t_s(base,direct_index,domain,vartag*100,z,size,(ZoneBc->Data_Set+nv),slab);
			}
		}
	}
	return 0;
}



// //This is same as  above////
// /**
// The function velocitymassBc() impliments the boundary conditions for massvelocities. It applies the same boundary
// conditions entered by the user for velocity components, to their respective mass velocity components. 
//  */
// 
// int velocitymassBc(Base *base,int nv,Zone *z,ZoneBoundary *ZoneBc)
// {
// 	char variable_name;
// 	int ni,nj,nk,direct_index,size;
// 	int var_index;
// 	
// 	ZoneBc->NueRms=0.0;
// 	ni=*(z->Size+0);
// 	nj=*(z->Size+1);
// 	nk=*(z->Size+2); 
// 	
// 	Find_Direction(&(ZoneBc->Direction),ZoneBc->Point_Range);
// 	direct_index=ZoneBc->Direction;
// 	
// 	var_index=nv-1;
// 		
// 	/*looping for all directions starts here*/	
// 	if (direct_index==1||direct_index==2)
// 	{
// 		size=(nk-1)*(nj-1);/*Size of the boundary patch*/
// 	}
// 	else if (direct_index==3||direct_index==4)
// 	{
// 		size=(nk-1)*(ni-1);
// 	}
// 	else if (direct_index==5||direct_index==6)
// 	{
// 		size=(nj-1)*(ni-1);
// 	}
// 
// 	variable_name=(ZoneBc->Data_Set+var_index)->Name_of_variable[0];
// 	switch(variable_name)
// 	{
// 		/*test for all variables starts here*/
// 		case 'u':
// 			if((base->flagTurbModel == 1 || base->flagTurbModel == 3) && z->flagSolidFluid == 0)
// 				ZoneBc->NueRms = boundary_d_t_sTurb(base,ZoneBc,direct_index,z->Flowsol->umass,nv*100,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->umass));
// 			else
// 				ZoneBc->NueRms = boundary_d_t_s(base,direct_index,z->Flowsol->umass,nv*100,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->umass));
// 			break;
// 
// 		case 'v':
// 			if((base->flagTurbModel == 1 || base->flagTurbModel == 3) && z->flagSolidFluid == 0)
// 				ZoneBc->NueRms = boundary_d_t_sTurb(base,ZoneBc,direct_index,z->Flowsol->vmass,nv*100,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->vmass));
// 			else
// 				ZoneBc->NueRms = boundary_d_t_s(base,direct_index,z->Flowsol->vmass,nv*100,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->vmass));
// 			break;
// 
// 		case 'w':
// 			if((base->flagTurbModel == 1 || base->flagTurbModel == 3) && z->flagSolidFluid == 0)
// 				ZoneBc->NueRms = boundary_d_t_sTurb(base,ZoneBc,direct_index,z->Flowsol->wmass,nv*100,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->wmass));
// 			else
// 				ZoneBc->NueRms = boundary_d_t_s(base,direct_index,z->Flowsol->wmass,nv*100,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->wmass));
// 			break;
// 
// 		case 'P':
// 			ZoneBc->NueRms=boundary_d_t_s(base,direct_index,z->Flowsol->P,nv,z,size,(ZoneBc->Data_Set+var_index),((ZoneBc->ZSlab)->P));
// 			break;
// 
// 	}
// 
// 	return(0);
// }
// 

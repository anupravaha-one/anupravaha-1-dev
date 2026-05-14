#include "APKUpdate_Boundary_Charge.h"
/**
This module inplements the boundary condition for charge density on wire.
This function in turn calls ImplementChargeWireBC() for implementing/updating
wire boundary value of charge density after every 'N' timesteps.
*/

int Update_Boundary_Charge(Base *base,int nv,Zone *z,ZoneBoundary *ZoneBc)
{
	FILE *sav,*Tstep;

	char variable_name;
	double Omega, Gama;
	int maxn,nki,njk,nij,ni,nj,nk,direct_index,size;
	int ret,Timestep, Criteria,cntTime; 

	ni=*(z->Size+0);
	nj=*(z->Size+1);
	nk=*(z->Size+2); 
		
	
	if((base->flagEfield==1)||(base->flagEfield==2))
	{

		Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","r");
		fscanf(Tstep,"%d",&Timestep);
		fclose(Tstep);
		
		cntTime = Timestep;
		
		sav=fopen("Settings/ESPInput/WireChargeUpdateData.dat","r");
		fscanf(sav,"%lf",&Omega);
		fscanf(sav,"%lf",&Gama);
		fscanf(sav,"%d",&Criteria);
		fclose(sav);

	}
	
	nv = find_nv(base,nv);

//End of updation of variable tags
		nij=ni*nj;
		nki=nk*ni;
		njk=nj*nk;
		if(nij>njk)
		{
		if(nij>nki)
		maxn=nij;
		else
		maxn=nki;
		}
		else
		{	if(njk>nki)
			maxn=njk;

			else 
			maxn=nki;
		}

		direct_index=ZoneBc->Direction;
		ZoneBc->NueRms=0.0; 
	/*looping for all directions starts here*/
		if (direct_index==1||direct_index==2)//east or west
		{
			size=(nk-1)*(nj-1);/*Size of the boundary patch*/
			
			variable_name=(ZoneBc->Data_Set+nv)->Name_of_variable[0];

			switch(variable_name)
			{	
				case 'I':
					Timestep = (Timestep%Criteria);
					if (Timestep==0)
					{
						ret = ImplementChargeWireBC(base,direct_index,z->Flowsol->RhoIon,9,z,size,(ZoneBc->Data_Set+nv),((ZoneBc->ZSlab)->RhoIon));
						WriteWireCharge(direct_index,z,size,((ZoneBc->ZSlab)->RhoIon),cntTime);
					}
				break;

			}
			
		}
		else if (direct_index==3||direct_index==4)
		{
			
			size=(nk-1)*(ni-1);
			variable_name=(ZoneBc->Data_Set+nv)->Name_of_variable[0];

			switch(variable_name)
			{	
				case 'I':
					Timestep = (Timestep%Criteria);
					if (Timestep==0)
					{
						ret = ImplementChargeWireBC(base,direct_index,z->Flowsol->RhoIon,9,z,size,(ZoneBc->Data_Set+nv),((ZoneBc->ZSlab)->RhoIon));
						WriteWireCharge(direct_index,z,size,((ZoneBc->ZSlab)->RhoIon),cntTime);
					}
				break;

			}
		}
		else if (direct_index==5||direct_index==6)
		{
			
			size=(nj-1)*(ni-1);
			variable_name=(ZoneBc->Data_Set+nv)->Name_of_variable[0];

			switch(variable_name)
			{	
				case 'I':
					Timestep = (Timestep%Criteria);
					if (Timestep==0)
					{
						ret = ImplementChargeWireBC(base,direct_index,z->Flowsol->RhoIon,9,z,size,(ZoneBc->Data_Set+nv),((ZoneBc->ZSlab)->RhoIon));
						WriteWireCharge(direct_index,z,size,((ZoneBc->ZSlab)->RhoIon),cntTime);
					}
				break;

			}
		}

	
		return 0;
}


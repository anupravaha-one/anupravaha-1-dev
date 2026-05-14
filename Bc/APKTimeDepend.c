#include "APKTimeDepend.h"
/**
This function TimeDepend() is called for implementation of Time dependent Boundary conditions for any 
variables. It takes time, variable tag and Cgnsbase pointers.  Various Switch Cases are removed to make the function simpler.
*/
int TimeDepend(double time,int vartag,Base *base)
{
	
	int numZ,numbcP,numZL,numbcPL,nv;
	
	numZ=base->Number_of_zones;
	nv=find_nv(base,vartag);
	for(numZL=0;numZL<numZ;numZL++)
	{
		
		numbcP=(base->CgnsZone+numZL)->Number_of_BC_Patches;
		for(numbcPL=0;numbcPL<numbcP;numbcPL++)
		{
			if((((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->flagTime==1)
			{
			   if ( ((((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->FunctionIndex) != NULL ) {
				Call_user_func(*((((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->FunctionIndex),  (((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->funct,     (base->CgnsZone+numZL),((base->CgnsZone+numZL)->ZoneBc+numbcPL),time);	
				Call_user_func(*((((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->FunctionIndex+1),(((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->funct_beta,(base->CgnsZone+numZL),((base->CgnsZone+numZL)->ZoneBc+numbcPL),time);
				Call_user_func(*((((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->FunctionIndex+2),(((base->CgnsZone+numZL)->ZoneBc+numbcPL)->Data_Set+nv)->funct_gama,(base->CgnsZone+numZL),((base->CgnsZone+numZL)->ZoneBc+numbcPL),time);
			   }
			}
		}
	}
	return(0);
}

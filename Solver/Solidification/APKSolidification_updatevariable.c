#include "APKSolidification_updatevariable.h"
int Solidification_updatevariable(Base *base)
{
	int numzones,i,total,NI,NJ,NK,no,numBc,numiface,index;
	int number,direction;
	Zone *zone;
	int *size,*Point_Range;
	double E_sm,E_fm,Melt_Temp,CP_s,CP_f,K_s,K_f,Rho_s,Rho_f,Mu_f,Mu_s,f;
	numzones = base->Number_of_zones;
	E_sm = base->BaseParameters->E_sm;
	E_fm = base->BaseParameters->E_fm;
	Melt_Temp = base->BaseParameters->Melt_Temp;
	CP_s = base->BaseParameters->CP_s;
	CP_f = base->BaseParameters->CP_f;
	K_s = base->BaseParameters->K_s;
	K_f = base->BaseParameters->K_f;
	Rho_s = base->BaseParameters->Rho_s;
	Rho_f = base->BaseParameters->Rho_f;
	Mu_s = base->BaseParameters->Mu_s;
	Mu_f = base->BaseParameters->Mu_f;
	
	for (no=0;no<numzones;no++)
	{
		zone = (base->CgnsZone+no);
		size=(base->CgnsZone+no)->Size;
		NI = *(size);
		NJ = *(size+1);
		NK = *(size+2);
		total = (NI-1)*(NJ-1)*(NK-1);
		for (i=0;i<total;i++)
		{
			*((zone->Flowsol)->f+i) = *((zone->Flowsol)->Latent_Heat_Content+i)/(E_fm - E_sm);
			f = *((zone->Flowsol)->f+i);
			*((zone->Flowsol)->Diffcoeff+i) = (K_s + (K_f - K_s) *f) / (CP_s + (CP_f - CP_s) *f);
			*((zone->Flowsol)->SpecificHeat+i) = CP_s + (CP_f - CP_s) *f;
			//*((zone->Flowsol)->Rho+i) = Rho_s + (Rho_f - Rho_s) *f;
		}
		numBc=zone->Number_of_BC_Patches;
		for (index=0;index<numBc;index++)
		{
			direction= (zone->ZoneBc+index)->Direction;
			Point_Range=((zone)->ZoneBc+index)->Point_Range;
			if(direction==1||direction==2)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
			if(direction==3||direction==4)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
			if(direction==5||direction==6)
			number=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

			
			for (i=0;i<number;i++)
			{
				*(((zone->ZoneBc+index)->ZSlab)->f+i) = *(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i)/(E_fm - E_sm);
				f = *(((zone->ZoneBc+index)->ZSlab)->f+i);
				*(((zone->ZoneBc+index)->ZSlab)->Diffcoeff+i) = (K_s + (K_f - K_s) *f) / (CP_s + (CP_f - CP_s) *f);
				*(((zone->ZoneBc+index)->ZSlab)->SpecificHeat+i) = CP_s + (CP_f - CP_s) *f;
				//*(((zone->ZoneBc+index)->ZSlab)->Rho+i) = Rho_s + (Rho_f - Rho_s) *f;
			}
		}
		
		numiface=zone->Number_of_Interfaces;
		for (index=0;index<numiface;index++)
		{
			direction= (zone->Iface+index)->Direction;
			Point_Range=((zone)->Iface+index)->Point_Range;
			if(direction==1||direction==2)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
			if(direction==3||direction==4)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
			if(direction==5||direction==6)
			number=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));
			for (i=0;i<number;i++)
			{
				*(((zone->Iface+index)->ZSlab)->f+i) = *(((zone->Iface+index)->ZSlab)->Latent_Heat_Content+i)/(E_fm - E_sm);
				f = *(((zone->Iface+index)->ZSlab)->f+i);
				*(((zone->Iface+index)->ZSlab)->Diffcoeff+i) =(K_s + (K_f - K_s) *f) / (CP_s + (CP_f - CP_s) *f);
				*(((zone->Iface+index)->ZSlab)->SpecificHeat+i) = CP_s + (CP_f - CP_s) *f;
				//*(((zone->Iface+index)->ZSlab)->Rho+i) = Rho_s + (Rho_f - Rho_s) *f;
				
			}
		}
		
	}
	return 0;	
}

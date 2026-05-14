#include "APKupdate_dele_iteration.h"
#define RELAX 1.0
/** This is the main subroutine for solidification and melting problems.
This function updates the latent heat content of each cell.
These function makes sure that solid only becomes liquid when it's latent heat content becomes
equal to latent heat of phase chage.
For details about this updation procedure refer Thesis of bhanesh on "NUMERICAL SOLIDIFICATION AND MELTING PROBLEMS"
*/

int update_dele_iteration(Base *base,int zone_no)
{
	int numzones,i,total,NI,NJ,NK,numBc,numiface,index;
	int number,direction;
	Zone *zone;
	int *size,*Point_Range;
	double E_sm,E_fm,update_dele_p;
	double f_inv,solidus,liquidus,melt,C_p;
	char tag,name[33];
	void *ptr;
	Interface *facePtr;
	Zone *ptrZoneNeighbor;
	
	numzones = base->Number_of_zones;
	E_sm = base->BaseParameters->E_sm;
	E_fm = base->BaseParameters->E_fm;
	if((base->BaseParameters->flag_mushy)==0)
	melt = base->BaseParameters->Melt_Temp;
	else if((base->BaseParameters->flag_mushy)==1)
	{
		solidus = base->BaseParameters->Solidus_Temp;
		liquidus= base->BaseParameters->Liquidus_Temp;
	}
	else
	{
		printf("wrong falg for solidification problem\nExiting in dele_iteration.\n");
		exit(0);
	}
	//for (no=0;no<numzones;no++)
	{
		zone = (base->CgnsZone+zone_no);
		size=(base->CgnsZone+zone_no)->Size;
		NI = *(size);
		NJ = *(size+1);
		NK = *(size+2);
		total = (NI-1)*(NJ-1)*(NK-1);
		for (i=0;i<total;i++)
		{
			//Fluid Fraction 'f' updation
			if((base->BaseParameters->flag_mushy)==0)
			f_inv = melt;
			else if((base->BaseParameters->flag_mushy)==1)
			f_inv = solidus + (*((zone->Flowsol)->Latent_Heat_Content+i)/(E_fm-E_sm)) *(liquidus-solidus); 
				
			C_p = *((zone->Flowsol)->SpecificHeat+i);
			update_dele_p = *((zone->Flowsol)->Latent_Heat_Content+i) + 
					RELAX*(*((zone->Flowsol)->enthalpy+i) - C_p * f_inv);
			//printf("cp =%e,melt=%e,E_sm =%e\n",C_p,f_inv,C_p * f_inv);
			//If any problem comes in the solution of the problem
			//check the above statement............... IMPORTANT
			if (update_dele_p > (E_fm-E_sm))
			update_dele_p = E_fm - E_sm;
			else if (update_dele_p < 0.0)
			update_dele_p = 0.0;
			//printf("%e\t",update_dele_p);
			*((zone->Flowsol)->Latent_Heat_Content+i) = update_dele_p;
			
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
				C_p = *(((zone->ZoneBc+index)->ZSlab)->SpecificHeat+i);
				if((base->BaseParameters->flag_mushy)==0)
				f_inv = melt;
				else if((base->BaseParameters->flag_mushy)==1)
				f_inv = solidus + (*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i)/(E_fm-E_sm)) *(liquidus-solidus); 
			
				update_dele_p = *(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i) + 
					RELAX*(*(((zone->ZoneBc+index)->ZSlab)->enthalpy+i) - C_p * f_inv);
				
				if (update_dele_p > (E_fm-E_sm))
				update_dele_p = E_fm - E_sm;
				else if (update_dele_p < 0.0)
				update_dele_p = 0.0;
				*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i)=update_dele_p;
				/*if(*(((zone->ZoneBc+index)->ZSlab)->enthalpy+i) <=E_sm)
				*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i) = 0.0;
				else if(*(((zone->ZoneBc+index)->ZSlab)->enthalpy+i)>=E_fm) 
				*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i) = E_fm - E_sm;
				else
				*(((zone->ZoneBc+index)->ZSlab)->Latent_Heat_Content+i) =*(((zone->ZoneBc+index)->ZSlab)->enthalpy+i) - E_sm;*/
			}

		}
		
		numiface=zone->Number_of_Interfaces;
		for (index=0;index<numiface;index++)
		{
			
			direction= (zone->Iface+index)->Direction;
			ptr = (void *)FindInterfaceBCPtrMod(zone,direction,name,&tag);
			if (tag!='i')
			{
				printf("Some thing wrong in update_dele_iteration\n");
				exit(1);
			}
			facePtr=(Interface *)ptr;
			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);
			FillerInterfaceVariables(21,facePtr,ptrZoneNeighbor,direction-1);
			
			/*Point_Range=((zone)->Iface+index)->Point_Range;
			if(direction==1||direction==2)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
			if(direction==3||direction==4)
			number=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
			if(direction==5||direction==6)
			number=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));
			for (i=0;i<number;i++)
			{
				
				if(*(((zone->Iface+index)->ZSlab)->enthalpy+i)<=E_sm)
				*(((zone->Iface+index)->ZSlab)->Latent_Heat_Content+i)=0.0;
				else if(*(((zone->Iface+index)->ZSlab)->enthalpy+i)>=E_fm)
				*(((zone->Iface+index)->ZSlab)->Latent_Heat_Content+i)=E_fm-E_sm;
				else
				*(((zone->Iface+index)->ZSlab)->Latent_Heat_Content+i)=*(((zone->Iface+index)->ZSlab)->enthalpy+i) - E_sm;           
			}*/
		}
	}
	return 0;	
}

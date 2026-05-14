#include "APKWriteCgnsFileDisperse.h"
/**
This function WriteCgnsFileDisperse() is called from the function WriteCgnsFile() if disperse phase flow is choosen. 
This function writes the disperse phase flow variables in CGNS file. 
*/
int WriteCgnsFileDisperse(char* CgnsFilename,Root CgnsRoot,int index_file,int sol_index ,int index_base,int index_zone,int solarr_index)
{	
	 DataType_t soldtype=4;
	int ier;

	//Disperse phase flow problem*****************************************************************************/
	
	// For Disperse phase momentum
	
	char *temp;
	//char *tempx1[10],*tempy1[10],*tempz1[10];
	//char *tempx1;
/*	char tempx[]="Dispersephasemomentum-X";
	char tempy[]="Dispersephasemomentum-Y";
	char tempz[] ="Dispersephasemomentum-Z";*/
	char tempvisc[]="Disperseviscosity";
	char tempden[] ="Dispersedensity";
 	char tempVolFrac[]="VolumeFraction-Disp";
	//char tempNumDen[]="Number-Density";
	char tempU[]="DispersephaseVelocity-X";
	char tempV[]="DispersephaseVelocity-Y";
	char tempW[]="DispersephaseVelocity-Z";
	char tempDia[]="Particle-Diameter";
	char tempFrac[]="FuelMassFraction";
	char tempTdrop[]="DropletTemperature";
	char tempCPdrop[]="DropletSpecHeat";
	temp=(char *)malloc(40*sizeof(char));
	
	/*for(i=0;i<10;i++)
	{	
		printf("\n\n\n i=%d\n\n",i);
		sprintf(temp,"%d",i);
		tempx1[i] =  strcat(tempx,temp);
		tempy1[i]= strcat(tempy,temp);
		tempz1[i]= strcat(tempz,temp);
	}*/
 
	printf("\nDisperse\n");
/*  
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MomentumDispersex!=NULL)	
		{
			
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MomentumDispersex,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		
		//tempx1=  strcat(tempy,temp);
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MomentumDispersey!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempy,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MomentumDispersey,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		//tempx1=  strcat(tempz,temp);
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MomentumDispersez!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->MomentumDispersez,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
*/		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->U_Disperse!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempU,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->U_Disperse,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->V_Disperse!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempV,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->V_Disperse,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->W_Disperse!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempW,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->W_Disperse,&solarr_index);	
			CGNSErrorCheck(ier);		
		}	

		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->VarthetaDisperse!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempVolFrac,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->VarthetaDisperse,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		//For Number Density
/*		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->NumDensity!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempNumDen,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->NumDensity,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
*/
		//for fuel mass fraction
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->DropFuelMassFrac!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempFrac,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->DropFuelMassFrac,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tdrop!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempTdrop,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Tdrop,&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagCPdrop==1)
			{
				if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->CPdrop!=NULL)
				{
					ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempCPdrop,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->CPdrop,&solarr_index);
					CGNSErrorCheck(ier);
				}
			}
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ParticleDia!=NULL)
		{
		if((CgnsRoot.CgnsBase+(index_base-1))->BaseParameters->flagVariableDia==1)
			{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempDia,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->PrevFlowsol->ParticleDia,&solarr_index);	
			CGNSErrorCheck(ier);		
			}
		}
				
	if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagRhoDisperse==1)
	{	
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoDisperse!=NULL)
		{
			
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempden,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoDisperse,&solarr_index);	
			CGNSErrorCheck(ier);
			
		}
		
	}	
	if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagViscoeffDisperse==1)	
	{
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ViscosityDisperse!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,soldtype,tempvisc,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ViscosityDisperse,&solarr_index);	
			CGNSErrorCheck(ier);
		}
	}	
		
		

	return(0);
}



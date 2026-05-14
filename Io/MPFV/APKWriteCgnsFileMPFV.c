#include "APKWriteCgnsFileMPFV.h"
/**
 This function is used to write the  cgns  solution file for multiphase flow(bubbly).Here all the variables which has been solved for mutiphase flow has been considered.
 */
int WriteCgnsFileMPFV(char* CgnsFilename,Root CgnsRoot,int index_file,int sol_index ,int index_base,int index_zone,int solarr_index)
{	
	 //DataType_t soldtype=4;
	int nophase,ier,nop;

	//Multiphase flow problem*****************************************************************************/
	
	// For individual phase momentum
	
	char *temp;
	char *tempx,*tempy,*tempz;
	char *P1,*P2,*P3,*tempx1;
	temp=(char *)malloc(40*sizeof(char));
	tempx=(char *)malloc(40*sizeof(char));
	tempy=(char *)malloc(40*sizeof(char));
	tempz=(char *)malloc(40*sizeof(char));
	sprintf(tempx,"Multiphasemomentum-X");
	sprintf(tempy,"Multiphasemomentum-Y");
	sprintf(tempz,"Multiphasemomentum-Z");
	
	//char tempx[]="multiphasemomentum-X";
	//char tempy[]="multiphasemomentum-Y";
	//char tempz[] ="multiphasemomentum-Z";
	char tempvisc[]="Phaseviscosity";
	char tempden[] ="Phasedensity";
	char tempV[]="Vartheta";
	//temp=(char *)malloc(40*sizeof(char));
	/*for(i=0;i<10;i++)
	{	
		printf("\n\n\n i=%d\n\n",i);
		sprintf(temp,"%d",i);
		tempx1[i] =  strcat(tempx,temp);
		tempy1[i]= strcat(tempy,temp);
		tempz1[i]= strcat(tempz,temp);
	}*/
 
 
 //For overall momentum
 char ovmx[]="momentum-X";
 char ovmy[]="momentum-Y";
 char ovmz[]="momentum-Z";
 
 
 char ovmxs[]="momentumstar-X";
 char ovmys[]="momentumstar-Y";
 char ovmzs[]="momentumstar-Z";
	nop=CgnsRoot.CgnsBase->nphase;
//  	printf("\nnophase=%d\n",nop);
	P1 = tempx;
	P2 =tempy;
	P3 =tempz;		

	for(nophase=0;nophase<nop;nophase++)
	{	
// 		printf("\n\n\nnophase=%d\n\n\n",nophase);
// 		printf("\n\n\ninside  phase loop\n\n");
		
		sprintf(temp,"%d",nophase);

		strcat(tempx,temp);
		
		
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasex[nophase]!=NULL)	
		{
			
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,P1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasex[nophase],&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		
		strcat(P2,temp);
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasey[nophase]!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,P2,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasey[nophase],&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		strcat(P3,temp);
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasez[nophase]!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,P3,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumofphasez[nophase],&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		
		tempx1=  strcat(tempV,temp);
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->VarthetaPhase[nophase]!=NULL)
		{
			ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,tempx1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->VarthetaPhase[nophase],&solarr_index);	
			CGNSErrorCheck(ier);		
		}
		tempx1=  strcat(tempden,temp);
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagRhoPhase[nophase]==1)
		{	
			
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoofPhase[nophase]!=NULL)
			{
				
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,tempx1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->RhoofPhase[nophase],&solarr_index);	
				CGNSErrorCheck(ier);
				
			}
			
		}	
		
		tempx1=  strcat(tempvisc,temp);
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->flagViscoeffPhase[nophase]==1)	
		{
			if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ViscosityofPhase[nophase]!=NULL)
			{
				ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,tempx1,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->ViscosityofPhase[nophase],&solarr_index);	
				CGNSErrorCheck(ier);
			}
		}	
	}	
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumx!=NULL)
		{
		ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmx,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumx,&solarr_index);	
			CGNSErrorCheck(ier);
		ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmxs,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumxstar,&solarr_index);	
		CGNSErrorCheck(ier);
		}
	
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumy!=NULL)
		{
		ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmy,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumy,&solarr_index);	
			CGNSErrorCheck(ier);
		ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmys,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumystar,&solarr_index);	
		CGNSErrorCheck(ier);
		}
		
		if(((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumz!=NULL)
		{
		ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmz,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumz,&solarr_index);	
			CGNSErrorCheck(ier);
		ier = cg_field_write(index_file,index_base,index_zone,sol_index,4,ovmzs,((CgnsRoot.CgnsBase+(index_base-1))->CgnsZone+(index_zone-1))->Flowsol->Momentumzstar,&solarr_index);	
		CGNSErrorCheck(ier);
		}
return(0);
}




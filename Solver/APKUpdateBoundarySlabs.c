#include "APKUpdateBoundarySlabs.h"
/**
This function UpdateBoundarySlabs() is called to fill the previous timestep slab values (uOld ,umassOld, vOld, vmassOld etc)
with current timestep slab values(u,umass,v,vmass etc). These values are used in the implementation of convective outflow
boundary conditions and also in Crank Nicolson time stepping scheme.
*/
int UpdateBoundarySlabs(Zone *zone)
{
	int k,j,i,NK,NJ,NI,iface,*ptrSize;//,iScal;

	char tag,name[33];
	void *ptr;
	ZoneBoundary *bcPtr;
		
	ptrSize = zone->Size;
		
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	for(iface=0;iface<6;iface++)
	{
		ptr = FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		
		if(tag=='b')
		{
			bcPtr=(ZoneBoundary *)ptr;
			
			if((bcPtr->Data_Set->Bc_Type)==5)
			{
			 
			 
				
				switch(iface+1)
				{
					case 1:
					  //	  printf("\nI am in case 1\n");
						for(k=0;k<(NK-1);k++)
						for(j=0;j<(NJ-1);j++)
						{
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
								*(bcPtr->ZSlab->uOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->u+k*(NJ-1)+j);

							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
								*(bcPtr->ZSlab->vOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->v+k*(NJ-1)+j);

							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
								*(bcPtr->ZSlab->wOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->w+k*(NJ-1)+j);
							
							if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
								*(bcPtr->ZSlab->TOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->T+k*(NJ-1)+j);

							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							  	*(bcPtr->ZSlab->umassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->uOld+k*(NJ-1)+j);
								
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							  	*(bcPtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->vOld+k*(NJ-1)+j);	
								
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							        *(bcPtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->wOld+k*(NJ-1)+j);		

							//Added by  mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->U_Disperse+k*(NJ-1)+j);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->V_Disperse+k*(NJ-1)+j);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->W_Disperse+k*(NJ-1)+j);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NJ-1)+j);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NJ-1)+j);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Tdrop+k*(NJ-1)+j);
							/*for(iScal=0;iScal<10;iScal++)
							{
								if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
									*(bcPtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
							}*/
						}
					break;
					case 2:
					  //  printf("\nI am in case 2\n");
						for(k=0;k<(NK-1);k++)
						for(j=0;j<(NJ-1);j++)
						{
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
							{
								*(bcPtr->ZSlab->uOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->u+k*(NJ-1)+j);
							}
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
							{
								*(bcPtr->ZSlab->vOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->v+k*(NJ-1)+j);
							}
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
							{
								*(bcPtr->ZSlab->wOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->w+k*(NJ-1)+j);
							}
							
							if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
								*(bcPtr->ZSlab->TOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->T+k*(NJ-1)+j);
								
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							  	*(bcPtr->ZSlab->umassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->uOld+k*(NJ-1)+j);
								
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							  	*(bcPtr->ZSlab->vmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->vOld+k*(NJ-1)+j);	
								
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
								*(bcPtr->ZSlab->wmassOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->wOld+k*(NJ-1)+j);	

							//Added by mrunalini
	
						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->U_Disperse+k*(NJ-1)+j);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->V_Disperse+k*(NJ-1)+j);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->W_Disperse+k*(NJ-1)+j);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NJ-1)+j);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NJ-1)+j);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NJ-1)+j)=*(bcPtr->ZSlab->Tdrop+k*(NJ-1)+j);	

						}
					break;
					case 3:
					  //	  printf("\nI am in case 3\n");
						for(k=0;k<(NK-1);k++)
						for(i=0;i<(NI-1);i++)
						{
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
								*(bcPtr->ZSlab->uOld+k*(NI-1)+i)=*(bcPtr->ZSlab->u+k*(NI-1)+i);
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
								*(bcPtr->ZSlab->vOld+k*(NI-1)+i)=*(bcPtr->ZSlab->v+k*(NI-1)+i);
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
								*(bcPtr->ZSlab->wOld+k*(NI-1)+i)=*(bcPtr->ZSlab->w+k*(NI-1)+i);
							if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
								*(bcPtr->ZSlab->TOld+k*(NI-1)+i)=*(bcPtr->ZSlab->T+k*(NI-1)+i);
							
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
								*(bcPtr->ZSlab->umassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->uOld+k*(NI-1)+i);
								
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							  	*(bcPtr->ZSlab->vmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->vOld+k*(NI-1)+i);	
								
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
								*(bcPtr->ZSlab->wmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->wOld+k*(NI-1)+i);		

							if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+k*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+k*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+k*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+k*(NI-1)+i);

							/*for(iScal=0;iScal<10;iScal++)
							{
								if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
									*(bcPtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
									}*/
						}
					break;
					case 4:
					  //  printf("\nI am in case 4\n");
						for(k=0;k<(NK-1);k++)
						for(i=0;i<(NI-1);i++)
						{
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
								*(bcPtr->ZSlab->uOld+k*(NI-1)+i)=*(bcPtr->ZSlab->u+k*(NI-1)+i);
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
								*(bcPtr->ZSlab->vOld+k*(NI-1)+i)=*(bcPtr->ZSlab->v+k*(NI-1)+i);
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
								*(bcPtr->ZSlab->wOld+k*(NI-1)+i)=*(bcPtr->ZSlab->w+k*(NI-1)+i);
							
							if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
								*(bcPtr->ZSlab->TOld+k*(NI-1)+i)=*(bcPtr->ZSlab->T+k*(NI-1)+i);
							
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
								*(bcPtr->ZSlab->umassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->uOld+k*(NI-1)+i);
								
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							  	*(bcPtr->ZSlab->vmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->vOld+k*(NI-1)+i);	
								
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
								*(bcPtr->ZSlab->wmassOld+k*(NI-1)+i)=*(bcPtr->ZSlab->wOld+k*(NI-1)+i);	

							if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+k*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+k*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+k*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+k*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+k*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+k*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+k*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+k*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+k*(NI-1)+i);

							
							/*for(iScal=0;iScal<10;iScal++)
							{
								if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
									*(bcPtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
									}*/
						}
					break;
					case 5:
					  //  printf("\nI am in case 5\n");
						for(j=0;j<(NJ-1);j++)
						for(i=0;i<(NI-1);i++)
						{
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
								*(bcPtr->ZSlab->uOld+j*(NI-1)+i)=*(bcPtr->ZSlab->u+j*(NI-1)+i);
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
								*(bcPtr->ZSlab->vOld+j*(NI-1)+i)=*(bcPtr->ZSlab->v+j*(NI-1)+i);
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
								*(bcPtr->ZSlab->wOld+j*(NI-1)+i)=*(bcPtr->ZSlab->w+j*(NI-1)+i);
							
							if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
							       *(bcPtr->ZSlab->TOld+j*(NI-1)+i)=*(bcPtr->ZSlab->T+j*(NI-1)+i);
							
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							       *(bcPtr->ZSlab->umassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->uOld+j*(NI-1)+i);
								
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							       *(bcPtr->ZSlab->vmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->vOld+j*(NI-1)+i);	
								
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							       *(bcPtr->ZSlab->wmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->wOld+j*(NI-1)+i);		

							//Added by adhi+shash for disperse phase flow	//mrunalini

						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+j*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+j*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+j*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+j*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+j*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+j*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+j*(NI-1)+i);



							/*for(iScal=0;iScal<10;iScal++)
							{
								if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
									*(bcPtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
							}*/
						}
					break;
					case 6:

					  //	  printf("\nI am in case 6\n");	
			                        for(j=0;j<(NJ-1);j++)
						for(i=0;i<(NI-1);i++)
						{
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->uOld)!=NULL)
								*(bcPtr->ZSlab->uOld+j*(NI-1)+i)=*(bcPtr->ZSlab->u+j*(NI-1)+i);
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vOld)!=NULL)
								*(bcPtr->ZSlab->vOld+j*(NI-1)+i)=*(bcPtr->ZSlab->v+j*(NI-1)+i);
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wOld)!=NULL)
								*(bcPtr->ZSlab->wOld+j*(NI-1)+i)=*(bcPtr->ZSlab->w+j*(NI-1)+i);
							
							if((zone->Flowsol->T)!=NULL && (bcPtr->ZSlab->TOld)!=NULL)
								*(bcPtr->ZSlab->TOld+j*(NI-1)+i)=*(bcPtr->ZSlab->T+j*(NI-1)+i);
							
							if((zone->Flowsol->u)!=NULL && (bcPtr->ZSlab->umassOld)!=NULL)
							       *(bcPtr->ZSlab->umassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->uOld+j*(NI-1)+i);
								
							if((zone->Flowsol->v)!=NULL && (bcPtr->ZSlab->vmassOld)!=NULL)
							       *(bcPtr->ZSlab->vmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->vOld+j*(NI-1)+i);	
								
							if((zone->Flowsol->w)!=NULL && (bcPtr->ZSlab->wmassOld)!=NULL)
							       *(bcPtr->ZSlab->wmassOld+j*(NI-1)+i)=*(bcPtr->ZSlab->wOld+j*(NI-1)+i);		

							//Added by adhi+shash for disperse phase flow	//mrunalini

						if((zone->Flowsol->U_Disperse)!=NULL && (bcPtr->ZSlab->U_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->U_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->U_Disperse+j*(NI-1)+i);
							
						if((zone->Flowsol->V_Disperse)!=NULL && (bcPtr->ZSlab->V_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->V_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->V_Disperse+j*(NI-1)+i);	
						
						if((zone->Flowsol->W_Disperse)!=NULL && (bcPtr->ZSlab->W_DisperseOld)!=NULL)
							*(bcPtr->ZSlab->W_DisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->W_Disperse+j*(NI-1)+i);
														
						if((zone->Flowsol->VarthetaDisperse)!=NULL && (bcPtr->ZSlab->VarthetaDisperseOld)!=NULL)
							*(bcPtr->ZSlab->VarthetaDisperseOld+j*(NI-1)+i)=*(bcPtr->ZSlab->VarthetaDisperse+j*(NI-1)+i);
						
						if((zone->Flowsol->DropFuelMassFrac)!=NULL && (bcPtr->ZSlab->DropFuelMassFracOld)!=NULL)
							*(bcPtr->ZSlab->DropFuelMassFracOld+j*(NI-1)+i)=*(bcPtr->ZSlab->DropFuelMassFrac+j*(NI-1)+i);

						if((zone->Flowsol->Tdrop)!=NULL && (bcPtr->ZSlab->TdropOld)!=NULL)
							*(bcPtr->ZSlab->TdropOld+j*(NI-1)+i)=*(bcPtr->ZSlab->Tdrop+j*(NI-1)+i);


							/*for(iScal=0;iScal<10;iScal++)
							{
								if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
									*(bcPtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
							}*/
						}
					break;
				}
			}
		}
	}
	
	return 0;
}

#include "APKTURBUpdateOldZSlabsKEpsilon.h"
/**
This function is used to update the Oldtime step values of k,epsilon while solving for turbulent fluid flow
in the interface and boundary slabs.  This function is called in UpdateOldZSlabs().
*/
int UpdateOldZSlabsKEpsilon(Zone *zone)
{
	int k,j,i,NK,NJ,NI,iface,*ptrSize,iScal;
	char tag,name[33];
	void *ptr;
	ZoneBoundary *bcPtr;
	Interface *IfacePtr;
		
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
			
			switch(iface+1)
			{
				case 1:

					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 2:
					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(bcPtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 3:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 4:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 5:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
				case 6:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (bcPtr->ZSlab->ScalOld[iScal])!=NULL)
								*(bcPtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(bcPtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
			}
		}
		else
		if(tag=='i')
		{
			IfacePtr=(Interface *)ptr;
			
			switch(iface+1)
			{
				case 1:

					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 2:
					for(k=0;k<(NK-1);k++)
					for(j=0;j<(NJ-1);j++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NJ-1)+j)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NJ-1)+j);
						}
					}
				break;
				case 3:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 4:
					for(k=0;k<(NK-1);k++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+k*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+k*(NI-1)+i);
						}
					}
				break;
				case 5:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
				case 6:
					for(j=0;j<(NJ-1);j++)
					for(i=0;i<(NI-1);i++)
					{
						for(iScal=0;iScal<2;iScal++)
						{
							if((zone->Flowsol->Scal[iScal])!=NULL && (IfacePtr->ZSlab->ScalOld[iScal])!=NULL)
								*(IfacePtr->ZSlab->ScalOld[iScal]+j*(NI-1)+i)=*(IfacePtr->ZSlab->Scal[iScal]+j*(NI-1)+i);
						}
					}
				break;
			}
		}
	
	}
	
	return 0;
}




#include "APKUpdateInterfaceStruct.h"

/**
This function is used to calculate and store the Application Programm Inerface (API) parameters of the interface structure.
The function takes the Base poiter as an input argument. For all the interfaces of all the Zones it first calculates 
the Transformation Matrix (for index tranformations from present Zone to neighboring Zone) using Derive_Transformation_Matrix(). 
The assignment of direction to the interfaces based on their span of indices is done by Find_Direction(). The translation vector 
for cell index transformations is caclculated using FindCellTransformRule().
After termination of the interface loop the function loops over all the boundaries in a Zone and assigns them a direction based on 
the spanof their indices.
*/
int UpdateInterfaceStruct(Base *CgnsBase)
{
	FILE *k;
	k=fopen("Settings/transform_matrix.dat","w");
	/*printf("i'm in update func \n");*/
	int numZones,iFace;
	int nZones,numIFace,numBC;
		
	nZones=CgnsBase->Number_of_zones;
	fprintf(k,"Number of Zones:%d\n",nZones);

	for(numZones=0;numZones<nZones;numZones++)
	{
		numIFace=(CgnsBase->CgnsZone+numZones)->Number_of_Interfaces;
		/*printf("\n\nzone no =%d , no of interfaces =%d\n",numZones+1,numIFace);*/
		fprintf(k,"##################");
		fprintf(k,"\nZONE No =%d , No of Interfaces =%d\n",numZones+1,numIFace);
		for(iFace=0;iFace<numIFace;iFace++)
		{
			/*printf(" iFace=%d\n",iFace+1);*/
			Derive_Transformation_Matrix((CgnsBase->CgnsZone+numZones)->Iface+iFace); 
			/*printf("i am goin to find dir func.\n");*/
			Find_Direction(&(((CgnsBase->CgnsZone+numZones)->Iface+iFace)->Direction),((CgnsBase->CgnsZone+numZones)->Iface+iFace)->Point_Range);
			//printf("%d\n",(((CgnsBase->CgnsZone+numZones)->Iface+iFace)->Direction));
			/*printf("i am out of Finddir func.\n");*/
			
			FindCellTransformRule((CgnsBase->CgnsZone+numZones)->Iface+iFace);
			/*printf("i am out of FindcellTransMatr\n");*/
			      
		}
		
		numBC=(CgnsBase->CgnsZone+numZones)->Number_of_BC_Patches;
		for(iFace=0;iFace<numBC;iFace++)
		{
			Find_Direction(&(((CgnsBase->CgnsZone+numZones)->ZoneBc+iFace)->Direction),((CgnsBase->CgnsZone+numZones)->ZoneBc+iFace)->Point_Range);
			//printf("%d\n",(((CgnsBase->CgnsZone+numZones)->ZoneBc+iFace)->Direction));
		}
	}
	fclose(k);
	return 0;
}


#include "APKMakeSlabPtrArray.h"
/** 
This function is modified completely from version 3.0.  The reasons are explained in Profiling and Optimization manual.  Refer documentation of APKMakeSlabPtrArray_com.c for more details.
**/


int MakeSlabPtrArray(Zone *zone,double *ptrArray[6],int vartag)
{
	int cntDir;

	switch (vartag)
	{
		case 1:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabU[cntDir];
			break;
		
		case 10:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabU_old[cntDir];
			break;
			
		case 100:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabUmass[cntDir];
			break;
		
		case 2:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabV[cntDir];
			break;
		
		case 20:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabV_old[cntDir];
			break;
			
		case 200:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVmass[cntDir];
			break;
			
		case 3:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabW[cntDir];
			break;
			
		case 30:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabW_old[cntDir];
			break;	
		case 300:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabWmass[cntDir];
			break;
		
		case 4:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabPress[cntDir];
			break;
		
		case 5:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVol[cntDir];
			break;
		case 6:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabT[cntDir];
			break;
		case 7:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Enthalpy[cntDir];
			break;
		case 8:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Potential[cntDir];
			break;
		case 9:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_RhoIon[cntDir];
			break;
		case 14:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabMixFrac[cntDir];
			break;
		case 17:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Potential_MHD[cntDir];
			break;
		case 19:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_H[cntDir];
			break;	
			
		case 56:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Total_MfieldX[cntDir];
			break;
		case 57:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Total_MfieldY[cntDir];
			break;
		case 58:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Total_MfieldZ[cntDir];
			break;
		case 59:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Psi_MHD[cntDir];
			break;
		case 990:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Flow_Potential[cntDir];
			break;
	//Changed by Mrunalini	
		case 62:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxDisperse[cntDir];
			break;
		case 620:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabVxDisperse[cntDir];
			break;
		case 6200:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabUDmass[cntDir];
			break;
		case 63:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVyDisperse[cntDir];
			break;	
		case 630:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabVyDisperse[cntDir];
			break;	
		case 6300:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVDmass[cntDir];
			break;	
		case 64:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzDisperse[cntDir];
			break;
		case 640:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabVzDisperse[cntDir];
			break;
		case 6400:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabWDmass[cntDir];
			break;
		case 65:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVarthetaDisperse[cntDir];
			break;
		case 650:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabVarthetaDisperse[cntDir];
			break;
		case 66:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabRhoDisperse[cntDir];
			break;
		case 660:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabRhoDisperse[cntDir];
			break;	
/*		case 67:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabNumDensity[cntDir];
			break;	
*/
		/*	case 670:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabNumDensity[cntDir];
			break;	
*/
		case 68:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir]=zone->Flowsol->SlabDropFuelMassFrac[cntDir];
			break;
		case 680:
			for(cntDir=0;cntDir<6;cntDir++)
			     ptrArray[cntDir]=zone->PrevFlowsol->SlabDropFuelMassFrac[cntDir];
			break;
		case 69:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir]=zone->Flowsol->SlabTdrop[cntDir];
			break;
		case 690:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir]=zone->PrevFlowsol->SlabTdrop[cntDir];



			//Added by sushil for multiphase
		case 31:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVx[cntDir];
			break;	
		case 310://Added by chandrasekhar , for 320 and 330 also
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVx_Old[cntDir];
			break;
		case 3100:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxstar[cntDir];
			break;
		case 32:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVy[cntDir];
			break;
		case 320:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVy_Old[cntDir];
			break;
		case 3200:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVystar[cntDir];
			break;
		case 33:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVz[cntDir];
			break;
		case 330:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVz_Old[cntDir];
			break;
		case 3300:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzstar[cntDir];
			break;	
		case 34:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVartheta[0][cntDir];
			break;
		case 35:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVartheta[1][cntDir];
			break;	
		case 36:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVartheta[2][cntDir];
			break;	
		case 37:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVartheta[3][cntDir];
			break;	
		case 38:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVartheta[4][cntDir];
			break;	
		case 41:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxPhase[0][cntDir];
			break;	
		case 42:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVyPhase[0][cntDir];
			break;	
		case 43:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzPhase[0][cntDir];
			break;	
		case 44:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxPhase[1][cntDir];
			break;	
		case 45:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVyPhase[1][cntDir];
			break;	
		case 46:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzPhase[1][cntDir];
			break;
		case 47:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxPhase[2][cntDir];
			break;	
		case 48:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVyPhase[2][cntDir];
			break;	
		case 49:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzPhase[2][cntDir];
			break;	
		case 50:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxPhase[3][cntDir];
			break;	
		case 51:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVyPhase[3][cntDir];
			break;
		case 52:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzPhase[3][cntDir];
			break;
				
		case 53:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVxPhase[4][cntDir];
			break;
		case 54:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVyPhase[4][cntDir];
			break;
		case 55:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVzPhase[4][cntDir];
			break;		
		case 410:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVcx_Old[cntDir];
			break;
		case 420:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVcy_Old[cntDir];
			break;
		case 430:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabVcz_Old[cntDir];
			break;
			//End.........
		case 600:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[0][cntDir];
			break;		
		case 601:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[1][cntDir];
			break;		
		case 602:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[2][cntDir];
			break;		
		case 603:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[3][cntDir];
			break;		
		case 604:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[4][cntDir];
			break;		
		case 605:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[5][cntDir];
			break;		
		case 606:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[6][cntDir];
			break;		
		case 607:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[7][cntDir];
			break;		
		case 608:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[8][cntDir];
			break;		
		case 609:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar[9][cntDir];
			break;		
			
		case 6000:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[0][cntDir];
			break;		
		case 6010:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[1][cntDir];
			break;		
		case 6020:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[2][cntDir];
			break;		
		case 6030:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[3][cntDir];
			break;		
		case 6040:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[4][cntDir];
			break;		
		case 6050:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[5][cntDir];
			break;		
		case 6060:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[6][cntDir];
			break;		
		case 6070:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[7][cntDir];
			break;		
		case 6080:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[8][cntDir];
			break;		
		case 6090:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabScalar_old[9][cntDir];
			break;		
			
		case 710:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDensity[cntDir];
			break;	
			
		case 711:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Visc[cntDir];
			break;	
		case 712:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_K[cntDir];
			break;		
		case 713:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_CP[cntDir];
			break;	
		case 714:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Eff_Visc[cntDir];
			break;	
		case 715:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->Slab_Tvisc[cntDir];
			break;
		case 726:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDummyViscosity_disperse[cntDir];
			break;
		case 727:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDiffVolfrac_disperse[cntDir];
			break;
		//Added by chandrasekhar --- for Turbulent MPFV
		case 730:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDummyeffViscoeffofphase[0][cntDir];
			break;
		case 731:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDummyeffViscoeffofphase[1][cntDir];
			break;
		case 732:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDummyeffViscoeffofphase[2][cntDir];
			break;
		case 733:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDummyeffViscoeffofphase[3][cntDir];
			break;
		case 734:
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->Flowsol->SlabDummyeffViscoeffofphase[3][cntDir];
			break;
		//End --- chandrasekhar
		
			break;
		/*case 7100: //Is it required.
			for(cntDir=0;cntDir<6;cntDir++)
				ptrArray[cntDir] = zone->PrevFlowsol->SlabCarrierDensity[cntDir];
			break;*/	
		default:
			printf("Till now it's not been implemented in MakeSlabPtrArray new.c for vartag=%d\n\n",vartag);
			printf("Exiting from this function\n\n");
			exit(12);
		break;
	}
	return 0;
}
		













// this one has been done for groundwater contaminant transport... 

int MakeSlabPtrArrayDeltaInverse(Zone *zone,double *ptrArray[6])
{
	int cntDir;
	void *ptr;
	char tag,zoneName[33];
	Interface *facePtr;
	ZoneBoundary *BCPtr;
	
	for(cntDir=0;cntDir<6;cntDir++)
	{
		ptr=(void *)FindInterfaceBCPtrMod(zone,cntDir+1,zoneName,&tag);
		if(tag=='i')
		{
			facePtr=(Interface *)ptr;
			ptrArray[cntDir]=facePtr->ZSlab->Delta;//FindSlabPtrIface(facePtr,vartag);
			//facePtr->ZSlab->Delta;
		}
		else
		if(tag=='b')
		{
			BCPtr=(ZoneBoundary *)ptr;
			ptrArray[cntDir]=BCPtr->ZSlab->Delta;// BCPtr;FindSlabPtrBC(BCPtr,vartag);
		}
		else
			exit(1);
	}
	return 0;
}

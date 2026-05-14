#include "APKMakeDomainArray.h"
/**
This function, MakeDomainArray() is to get the variable pointer from the 'Flowsol' data structure. The pointer of a variable
is required in the solver when it is solving for that variable.
*/

double *MakeDomainArray(Zone *zone,int vartag)
{
	double *domain;
	
	switch(vartag)
	{
		case 1:/**/
			domain=zone->Flowsol->u;
		break;
		case 10:/**/
		case 1000:
			domain=zone->PrevFlowsol->u;
		break;
		case 100:/**/
			domain=zone->Flowsol->umass;
		break;
		
		case 2:/**/
			domain=zone->Flowsol->v;
		break;
		case 20:/**/
		case 2000:
			domain=zone->PrevFlowsol->v;
		break;
		case 200:/**/
			domain=zone->Flowsol->vmass;
		break;
		case 21:/**/
			domain=zone->PrevFlowsol->Latent_Heat_Content;
		break;
		case 3:/**/
			domain=zone->Flowsol->w;
		break;
		case 30:/**/
		case 3000:
			domain=zone->PrevFlowsol->w;
		break;
		
		case 300:/**/
			domain=zone->Flowsol->wmass;
		break;
			
		case 4:/**/
			domain=zone->Flowsol->P;
		break;
		case 40:/**/
			domain=zone->PrevFlowsol->P;
		break;
		
		case 5:/**/
			domain=zone->ZGeom->Volume;
		break;
		
		case 6:
			domain=zone->Flowsol->T;
		break;
		
		case 60:
			domain=zone->PrevFlowsol->T;
		break;
		
		case 7:
			domain=zone->Flowsol->enthalpy;
		break;

		case 70:
			domain=zone->PrevFlowsol->enthalpy;
		break;
		
		//Electrical Parameters
		case 8:
			domain=zone->Flowsol->Potential;
		break;
		
		case 80:
			domain=zone->PrevFlowsol->Potential;
		break;
		
		case 9:
			domain=zone->Flowsol->RhoIon;
		break;
		
		case 90:
			domain=zone->PrevFlowsol->RhoIon;
		break;

		//end of Electrical parameters
		
                case 14:
			domain=zone->Flowsol->MixFrac;
		break;

		case 140:
			domain=zone->PrevFlowsol->MixFrac;
		break;
                case 15:
			domain=zone->Flowsol->FuelMassFrac;
		break;

		case 150:
			domain=zone->PrevFlowsol->FuelMassFrac;
		break;
		
		case 16:
			domain=zone->Flowsol->RadIntensity;
		break;
		
		case 160:
			domain=zone->PrevFlowsol->RadIntensity;
		break;	
		case 17://naveen
			domain=zone->Flowsol->PotentialMHD;
		break;
		
		case 170://naveen
			domain=zone->PrevFlowsol->PotentialMHD;
		break;
		
		case 19://groundwater
			domain=zone->Flowsol->H;
		break;
		
		case 190://groundwater
			domain=zone->PrevFlowsol->H;
		break;
		
		case 22://naveen
			domain=zone->Flowsol->UcrossBX;
		break;
		
		case 23://naveen
			domain=zone->Flowsol->UcrossBY;
		break;
		
		case 24://naveen
			domain=zone->Flowsol->UcrossBZ;
		break;
		
		case 25://naveen
			domain=zone->Flowsol->MfieldX;
		break;
		
		case 26://naveen
			domain=zone->Flowsol->MfieldY;
		break;
		
		case 27://naveen
			domain=zone->Flowsol->MfieldZ;
		break;
		case 56://rajendra
			domain=zone->Flowsol->TotalMfieldX;
		break;
		case 560://rajendra
			domain=zone->PrevFlowsol->TotalMfieldX;
		break;
		case 57://rajendra
			domain=zone->Flowsol->TotalMfieldY;
		break;
		case 570://rajendra
			domain=zone->PrevFlowsol->TotalMfieldY;
		break;
		case 58://rajendra
			domain=zone->Flowsol->TotalMfieldZ;
		break;
		case 580://rajendra
			domain=zone->PrevFlowsol->TotalMfieldZ;
		break;
		case 59://rajendra
			domain=zone->Flowsol->Psi;
		break;
		case 590://rajendra
			domain=zone->PrevFlowsol->Psi;
		break;
		//Added by sushil:For multiphase flow varibles......
		
		case 31:
			domain=zone->Flowsol->Momentumx;
		break;
		
		case 32:
			domain=zone->Flowsol->Momentumy;
		break;
		case 33:
			domain=zone->Flowsol->Momentumz;
		break;
		
				
		case 310:
			domain=zone->PrevFlowsol->Momentumx;
		break;
				
		case 320:
			domain=zone->PrevFlowsol->Momentumy;
		break;
		
		case 330:
			domain=zone->PrevFlowsol->Momentumz;
		break;
		
		case 3100:
	         	domain=zone->Flowsol->Momentumxstar;
		break;
		case 3200:
			domain=zone->Flowsol->Momentumystar;
		break;
		case 3300:
			domain=zone->Flowsol->Momentumzstar;
		break;
		
		
		case 41:
			domain=zone->Flowsol->Momentumofphasex[0];
		break;
		case 42:
			domain=zone->Flowsol->Momentumofphasey[0];
		break;
		case 43:
			domain=zone->Flowsol->Momentumofphasez[0];
		break;
		case 44:
			domain=zone->Flowsol->Momentumofphasex[1];
		break;
		case 45:
			domain=zone->Flowsol->Momentumofphasey[1];
		break;
		case 46:
			domain=zone->Flowsol->Momentumofphasez[1];
		break;
		case 47:
			domain=zone->Flowsol->Momentumofphasex[2];
		break;
		case 48:
			domain=zone->Flowsol->Momentumofphasey[2];
		break;
		case 49:
			domain=zone->Flowsol->Momentumofphasez[2];
		break;
		case 50:
			domain=zone->Flowsol->Momentumofphasex[3];
		break;
		case 51:
			domain=zone->Flowsol->Momentumofphasey[3];
		break;
		case 52:
			domain=zone->Flowsol->Momentumofphasez[3];
		break;
		case 53:
			domain=zone->Flowsol->Momentumofphasex[4];
		break;
		case 54:
			domain=zone->Flowsol->Momentumofphasey[4];
		break;
		case 55:
			domain=zone->Flowsol->Momentumofphasez[4];
		break;
		case 410:
			domain=zone->PrevFlowsol->Momentumofphasex[0];
		break;
		case 420:
			domain=zone->PrevFlowsol->Momentumofphasey[0];
		break;
		case 430:
			domain=zone->PrevFlowsol->Momentumofphasez[0];
		break;
		case 440:
			domain=zone->PrevFlowsol->Momentumofphasex[1];
		break;
		case 450:
			domain=zone->PrevFlowsol->Momentumofphasey[1];
		break;
		case 460:
			domain=zone->PrevFlowsol->Momentumofphasez[1];
		break;
		case 470:
			domain=zone->PrevFlowsol->Momentumofphasex[2];
		break;
		case 480:
			domain=zone->PrevFlowsol->Momentumofphasey[2];
		break;
		case 490:
			domain=zone->PrevFlowsol->Momentumofphasez[2];
		break;
		case 500:
			domain=zone->PrevFlowsol->Momentumofphasex[3];
		break;
		case 510:
			domain=zone->PrevFlowsol->Momentumofphasey[3];
		break;
		case 520:
			domain=zone->PrevFlowsol->Momentumofphasez[3];
		break;
		case 530:
			domain=zone->PrevFlowsol->Momentumofphasex[4];
		break;
		case 540:
			domain=zone->PrevFlowsol->Momentumofphasey[4];
		break;
		case 550:
			domain=zone->PrevFlowsol->Momentumofphasez[4];
		break;
		case 34:
			domain=zone->Flowsol->VarthetaPhase[0];
		break;
		
		case 35:
			domain=zone->Flowsol->VarthetaPhase[1];
		break;
		case 36:
			domain=zone->Flowsol->VarthetaPhase[2];
		break;
		case 37:
			domain=zone->Flowsol->VarthetaPhase[3];
		break;
		case 38:
			domain=zone->Flowsol->VarthetaPhase[4];
		break;
		case 340:
			domain=zone->PrevFlowsol->VarthetaPhase[0];
		break;
		
		case 350:
			domain=zone->PrevFlowsol->VarthetaPhase[1];
		break;
		case 360:
			domain=zone->PrevFlowsol->VarthetaPhase[2];
		break;
		case 370:
			domain=zone->PrevFlowsol->VarthetaPhase[3];
		break;
		case 380:
			domain=zone->PrevFlowsol->VarthetaPhase[4];
		break;
		
		case 720:
			domain=zone->Flowsol->ViscoeffDummy[0];
			break;
		case 721:
			domain=zone->Flowsol->ViscoeffDummy[1];
			break;
		case 722:
			domain=zone->Flowsol->ViscoeffDummy[2];
			break;
		case 723:
			domain=zone->Flowsol->ViscoeffDummy[3];
			break;
		case 724:
			domain=zone->Flowsol->ViscoeffDummy[4];
			break;
		
		//End of multiphase flow variables.......
		
		case 600:
			domain=zone->Flowsol->Scal[0];
		break;
		
		case 6000:
			domain=zone->PrevFlowsol->Scal[0];
		break;	
		
		case 601:
			domain=zone->Flowsol->Scal[1];
		break;
		
		case 6010:
			domain=zone->PrevFlowsol->Scal[1];
		break;	
		
		case 602:
			domain=zone->Flowsol->Scal[2];
		break;
		
		case 6020:
			domain=zone->PrevFlowsol->Scal[2];
		break;	
			
		case 603:
			domain=zone->Flowsol->Scal[3];
		break;
		
		case 6030:
			domain=zone->PrevFlowsol->Scal[3];
		break;	
		
		case 604:
			domain=zone->Flowsol->Scal[4];
		break;
		
		case 6040:
			domain=zone->PrevFlowsol->Scal[4];
		break;	
		
		case 605:
			domain=zone->Flowsol->Scal[5];
		break;
		
		case 6050:
			domain=zone->PrevFlowsol->Scal[5];
		break;	
		
		case 606:
			domain=zone->Flowsol->Scal[6];
		break;
		
		case 6060:
			domain=zone->PrevFlowsol->Scal[6];
		break;	
		
		case 607:
			domain=zone->Flowsol->Scal[7];
		break;
		
		case 6070:
			domain=zone->PrevFlowsol->Scal[7];
		break;	
		
		case 608:
			domain=zone->Flowsol->Scal[8];
		break;
		
		case 6080:
			domain=zone->PrevFlowsol->Scal[8];
		break;	
		
		case 609:
			domain=zone->Flowsol->Scal[9];
		break;
		
		case 6090:
			domain=zone->PrevFlowsol->Scal[9];
		break;
		
		case 700:
			domain=zone->Flowsol->Scalcoeff[0];
		break;
		
		case 701:
			domain=zone->Flowsol->Scalcoeff[1];
		break;
		
		case 702:
			domain=zone->Flowsol->Scalcoeff[2];
		break;
		
		case 703:
			domain=zone->Flowsol->Scalcoeff[3];
		break;
		
		case 704:
			domain=zone->Flowsol->Scalcoeff[4];
		break;
		
		case 705:
			domain=zone->Flowsol->Scalcoeff[5];
		break;
		
		case 706:
			domain=zone->Flowsol->Scalcoeff[6];
		break;
		
		case 707:
			domain=zone->Flowsol->Scalcoeff[7];
		break;
		
		case 708:
			domain=zone->Flowsol->Scalcoeff[8];
		break;
		
		case 709:
			domain=zone->Flowsol->Scalcoeff[9];
		break;
		
		case 710:
			domain=zone->Flowsol->Rho;
		break;
		
		case 711:
			domain=zone->Flowsol->Viscoeff;
		break;
		
		case 712:
			domain=zone->Flowsol->Diffcoeff;
		break;
		
		case 713:
			domain=zone->Flowsol->SpecificHeat;
		break;
		case 714:
			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;
		case 715:
			domain=zone->Flowsol->Tcoeff->TViscoeff;
		break;
                case 716:
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
		break;
		case 717:
			domain=zone->Flowsol->ChargeDiffcoeff;
		break;
		case 718:
			domain=zone->Flowsol->Tcoeff->TurbConductivity;
		break;
		case 719:
			domain=zone->Flowsol->UnSatPermeability	;
		break;
		
		//$$$$$$$$$$$$$$$$$$ for scaler in groundwater  Grad in X direction$$$$$$$$$$$$$$$$$$$$//
		
		case 60000:
			domain=zone->Flowsol->GradScalX[0];
		break;
		case 60001:
			domain=zone->Flowsol->GradScalX[1];
		break;
		
		case 60002:
			domain=zone->Flowsol->GradScalX[2];
		break;
		case 60003:
			domain=zone->Flowsol->GradScalX[3];
		break;
		
		case 60004:
			domain=zone->Flowsol->GradScalX[4];
		break;
		
		case 60005:
			domain=zone->Flowsol->GradScalX[5];
		break;
		
		case 60006:
			domain=zone->Flowsol->GradScalX[6];
		break;
		
		case 60007:
			domain=zone->Flowsol->GradScalX[7];
		break;
		
		case 60008:
			domain=zone->Flowsol->GradScalX[8];
		break;
		
		case 60009:
			domain=zone->Flowsol->GradScalX[9];
		break;
		
		//$$$$$$$$$$$$$$$$$$ for scaler in groundwater  Grad in Y direction$$$$$$$$$$$$$$$$$$$$//
				
		case 60100:
			domain=zone->Flowsol->GradScalY[0];
		break;
			
		case 60101:
			domain=zone->Flowsol->GradScalY[1];
		break;
		
		case 60102:
			domain=zone->Flowsol->GradScalY[2];
		break;
		case 60103:
			domain=zone->Flowsol->GradScalY[3];
		break;
		
		case 60104:
			domain=zone->Flowsol->GradScalY[4];
		break;
		
		case 60105:
			domain=zone->Flowsol->GradScalY[5];
		break;
		
		case 60106:
			domain=zone->Flowsol->GradScalY[6];
		break;
		
		case 60107:
			domain=zone->Flowsol->GradScalY[7];
		break;
		
		case 60108:
			domain=zone->Flowsol->GradScalY[8];
		break;
		
		case 60109:
			domain=zone->Flowsol->GradScalY[9];
		break;
		
		//$$$$$$$$$$$$$$$$$$ for scaler in groundwater  Grad in Z direction$$$$$$$$$$$$$$$$$$$$//
		case 60200:
			domain=zone->Flowsol->GradScalZ[0];
		break;
			
		case 60201:
			domain=zone->Flowsol->GradScalZ[1];
		break;
		
		case 60202:
			domain=zone->Flowsol->GradScalZ[2];
		break;
		case 60203:
			domain=zone->Flowsol->GradScalZ[3];
		break;
		
		case 60204:
			domain=zone->Flowsol->GradScalZ[4];
		break;
		
		case 60205:
			domain=zone->Flowsol->GradScalZ[5];
		break;
		
		case 60206:
			domain=zone->Flowsol->GradScalZ[6];
		break;
		
		case 60207:
			domain=zone->Flowsol->GradScalZ[7];
		break;
		
		case 60208:
			domain=zone->Flowsol->GradScalZ[8];
		break;
		
		case 60209:
			domain=zone->Flowsol->GradScalZ[9];
		break;

		case 990:
			domain=zone->Flowsol->Flow_Potential;
		break;

		case 9900:
			domain=zone->PrevFlowsol->Flow_Potential;
		break;
		
		
		//$$$$$$$$$$$$$$$$$$ end for sclar in groundwater$$$$$$$$$$$$$$$$$$$$//
		//Added by Adhi+shash :For Disperse phase flow varibles......
		
		case 62:
			domain=zone->Flowsol->U_Disperse;
		break;
		
		case 63:
			domain=zone->Flowsol->V_Disperse;
		break;
		case 64:
			domain=zone->Flowsol->W_Disperse;
		break;
		
				
		case 620:
		case 62000:
			domain=zone->PrevFlowsol->U_Disperse;
		break;		
		case 630:
		case 63000:
			domain=zone->PrevFlowsol->V_Disperse;
		break;
		
		case 640:
		case 64000:
			domain=zone->PrevFlowsol->W_Disperse;
		break;
		
		case 6200:
			domain=zone->Flowsol->U_Dmass;
		break;
		case 6300:
			domain=zone->Flowsol->V_Dmass;
		break;
		case 6400:
			domain=zone->Flowsol->W_Dmass;
		break;

		case 65:
			domain=zone->Flowsol->VarthetaDisperse;
		break;
		
		case 650:
			domain=zone->PrevFlowsol->VarthetaDisperse;
		break;
		case 66:
			domain=zone->Flowsol->RhoDisperse;
		break;
		case 660:
			domain=zone->PrevFlowsol->RhoDisperse;
		break;
/*		case 67:
			domain=zone->Flowsol->NumDensity;//shashwat
		break;
		case 670:
			domain=zone->PrevFlowsol->NumDensity;
		break;
*/
		case 68:
			domain=zone->Flowsol->DropFuelMassFrac;//shashwat
		break;
		case 680:
			domain=zone->PrevFlowsol->DropFuelMassFrac;
		break;
		case 69:
			domain=zone->Flowsol->Tdrop;//shashwat
		break;
		case 690:
			domain=zone->PrevFlowsol->Tdrop;
		break;
		case 7100:
			domain=zone->PrevFlowsol->Rho;
		break;
		case 725:
			domain=zone->Flowsol->CPdrop;
		break;
		case 726:
			domain=zone->Flowsol->Tcoeff->DummyViscoeff_disperse;
		break;
		case 727:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac_disperse;
		break;
		//end for disperse phase by adhiraj			
		default:
			printf("Tag(%d) is not present in Make domain Array.c\nexiting....\n",vartag);
			exit(0);
		break;	
	}
	return (domain);
}

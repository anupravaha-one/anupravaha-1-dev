#include "APKFindSlabPtrBC.h"
/**
This function, FindSlabPtrBC() gives boundary slab pointer value of a variable, when its 'vartag' is supplied.
*/

double *FindSlabPtrBC(ZoneBoundary *facePtr,int vartag)
{
	double *slab;
	
	switch(vartag)
	{
		case 1:/**/
			slab=facePtr->ZSlab->u;
		break;
		case 10:/**/
			slab=facePtr->ZSlab->uOld;
		break;
		case 100:/**/
			slab=facePtr->ZSlab->umass;
		break;
		case 1000:/**/
		  	slab=facePtr->ZSlab->umassOld;
                     
		break;

		case 2:/**/
			slab=facePtr->ZSlab->v;
		break;
		case 20:/**/
			slab=facePtr->ZSlab->vOld;
		break;
		
		case 200:/**/
			slab=facePtr->ZSlab->vmass;
		break;
		case 2000:/**/
		  	slab=facePtr->ZSlab->vmassOld;
                       
		break;

		case 3:/**/
			slab=facePtr->ZSlab->w;
		break;
		case 30:/**/
			slab=facePtr->ZSlab->wOld;
		break;
		
		case 300:/**/
			slab=facePtr->ZSlab->wmass;
		break;
		case 3000:/**/
		  	slab=facePtr->ZSlab->wmassOld;
                      
		break;
		
		case 4:/**/
			slab=facePtr->ZSlab->P;
		break;
		case 5:/**/
			slab=facePtr->ZSlab->Volume;
		break;
		case 6:/**/
			slab=facePtr->ZSlab->T;
		break;
		case 60:/**/
			slab=facePtr->ZSlab->TOld;
		break;
		
		case 7:/**/
			slab=facePtr->ZSlab->enthalpy;
		break;
		case 70:/**/
			slab=facePtr->ZSlab->enthalpyOld;
		break;
		
		//Electrical Parameters
		case 8:/**/
			slab=facePtr->ZSlab->Potential;
		break;
		case 80:/**/
			slab=facePtr->ZSlab->PotentialOld;
		break;
		case 9:/**/
			slab=facePtr->ZSlab->RhoIon;
		break;
		case 90:/**/
			slab=facePtr->ZSlab->RhoIonOld;
		break;
		//End of Electrical Parameters
		
                case 14:/**/
			slab=facePtr->ZSlab->MixFrac;
		break;
		case 140:/**/
			slab=facePtr->ZSlab->MixFracOld;
		break;
                case 15:/**/
			slab=facePtr->ZSlab->FuelMassFrac;
		break;
		case 150:/**/
			slab=facePtr->ZSlab->FuelMassFracOld;
		break;
		
		case 16:
			slab=facePtr->ZSlab->RadIntensity;
		break;
		
		case 17://naveen
			slab=facePtr->ZSlab->PotentialMHD;
		break;
		case 170://naveen
			slab=facePtr->ZSlab->PotentialMHDold;
		break;
		case 19:
			slab=facePtr->ZSlab->H;
		break;
		case 22://naveen
			slab=facePtr->ZSlab->UcrossBX;
		break;
		case 23://naveen
			slab=facePtr->ZSlab->UcrossBY;
		break;
		case 24://naveen
			slab=facePtr->ZSlab->UcrossBZ;
		break;
		case 25://naveen
			slab=facePtr->ZSlab->MfieldX;
		break;
		case 26://naveen
			slab=facePtr->ZSlab->MfieldY;
		break;
		case 27://naveen
			slab=facePtr->ZSlab->MfieldZ;
		break;
		case 56://rajendra
			slab=facePtr->ZSlab->TotalMfieldX;
		break;
		case 560://rajendra
			slab=facePtr->ZSlab->TotalMfieldXold;
		break;
		case 57://rajendra
			slab=facePtr->ZSlab->TotalMfieldY;
		break;
		case 570://rajendra
			slab=facePtr->ZSlab->TotalMfieldYold;
		break;
		case 58://rajendra
			slab=facePtr->ZSlab->TotalMfieldZ;
		break;
		case 580://rajendra
			slab=facePtr->ZSlab->TotalMfieldZold;
		break;
		case 59://rajendra
			slab=facePtr->ZSlab->Psi;
		break;
		case 590://rajendra
			slab=facePtr->ZSlab->Psiold;
		break;
		//By SUSHIL:For multiphase flow variables.............
		
		
		case 31:
			slab=facePtr->ZSlab->Momentumx;
		break;
		
		case 310:
			slab=facePtr->ZSlab->MomentumxOld;
		break;
		case 32:
			slab=facePtr->ZSlab->Momentumy;
		break;
		
		case 320:
			slab=facePtr->ZSlab->MomentumyOld;
		break;
		case 33:
			slab=facePtr->ZSlab->Momentumz;
		break;
		case 330:
			slab=facePtr->ZSlab->MomentumzOld;
		break;
		
		case 3100:
	         	slab=facePtr->ZSlab->Momentumxstar;
		break;
		case 3200:
			slab=facePtr->ZSlab->Momentumystar;
		break;
		case 3300:
			slab=facePtr->ZSlab->Momentumzstar;
		break;
		
		case 31000:
	         	slab=facePtr->ZSlab->MomentumxstarOld;
		break;
		case 32000:
			slab=facePtr->ZSlab->MomentumystarOld;
		break;
		case 33000:
			slab=facePtr->ZSlab->MomentumzstarOld;
		break;
		case 41:
			slab=facePtr->ZSlab->Momentumofphasex[0];
		break;
		case 42:
			slab=facePtr->ZSlab->Momentumofphasey[0];
		break;
		case 43:
			slab=facePtr->ZSlab->Momentumofphasez[0];
		break;
		case 44:
			slab=facePtr->ZSlab->Momentumofphasex[1];
		break;
		case 45:
			slab=facePtr->ZSlab->Momentumofphasey[1];
		break;
		case 46:
			slab=facePtr->ZSlab->Momentumofphasez[1];
		break;
		case 47:
			slab=facePtr->ZSlab->Momentumofphasex[2];
		break;
		case 48:
			slab=facePtr->ZSlab->Momentumofphasey[2];
		break;
		case 49:
			slab=facePtr->ZSlab->Momentumofphasez[2];
		break;
		case 50:
			slab=facePtr->ZSlab->Momentumofphasex[3];
		break;
		case 51:
			slab=facePtr->ZSlab->Momentumofphasey[3];
		break;
		case 52:
			slab=facePtr->ZSlab->Momentumofphasez[3];
		break;
		case 53:
			slab=facePtr->ZSlab->Momentumofphasex[4];
		break;
		case 54:
			slab=facePtr->ZSlab->Momentumofphasey[4];
		break;
		case 55:
			slab=facePtr->ZSlab->Momentumofphasez[4];
		break;
		case 410:
			slab=facePtr->ZSlab->MomentumofphasexOld[0];
		break;
		case 420:
			slab=facePtr->ZSlab->MomentumofphaseyOld[0];
		break;
		case 430:
			slab=facePtr->ZSlab->MomentumofphasezOld[0];
		break;
		case 440:
			slab=facePtr->ZSlab->MomentumofphasexOld[1];
		break;
		case 450:
			slab=facePtr->ZSlab->MomentumofphaseyOld[1];
		break;
		case 460:
			slab=facePtr->ZSlab->MomentumofphasezOld[1];
		break;
		case 470:
			slab=facePtr->ZSlab->MomentumofphasexOld[2];
		break;
		case 480:
			slab=facePtr->ZSlab->MomentumofphaseyOld[2];
		break;
		case 490:
			slab=facePtr->ZSlab->MomentumofphasezOld[2];
		break;
		case 500:
			slab=facePtr->ZSlab->MomentumofphasexOld[3];
		break;
		case 510:
			slab=facePtr->ZSlab->MomentumofphaseyOld[3];
		break;
		case 520:
			slab=facePtr->ZSlab->MomentumofphasezOld[3];
		break;
		case 530:
			slab=facePtr->ZSlab->MomentumofphasexOld[4];
		break;
		case 540:
			slab=facePtr->ZSlab->MomentumofphaseyOld[4];
		break;
		case 550:
			slab=facePtr->ZSlab->MomentumofphasezOld[4];
		break;
		case 34:
			slab=facePtr->ZSlab->VarthetaPhase[0];
		break;
		
		case 35:
			slab=facePtr->ZSlab->VarthetaPhase[1];
		break;
		case 36:
			slab=facePtr->ZSlab->VarthetaPhase[2];
		break;
		case 37:
			slab=facePtr->ZSlab->VarthetaPhase[3];
		break;
		case 38:
			slab=facePtr->ZSlab->VarthetaPhase[4];
		break;
		case 340:
			slab=facePtr->ZSlab->VarthetaPhaseOld[0];
		break;
		
		case 350:
			slab=facePtr->ZSlab->VarthetaPhaseOld[1];
		break;
		case 360:
			slab=facePtr->ZSlab->VarthetaPhaseOld[2];
		break;
		case 370:
			slab=facePtr->ZSlab->VarthetaPhaseOld[3];
		break;
		case 380:
			slab=facePtr->ZSlab->VarthetaPhaseOld[4];
		break;
		
		case 720:
			slab=facePtr->ZSlab->ViscoeffDummy[0];
			break;
		case 721:
			slab=facePtr->ZSlab->ViscoeffDummy[1];
			break;
		case 722:
			slab=facePtr->ZSlab->ViscoeffDummy[2];
			break;
		case 723:
			slab=facePtr->ZSlab->ViscoeffDummy[3];
			break;
		case 724:
			slab=facePtr->ZSlab->ViscoeffDummy[4];
			break;
		
			
		//Chandrasekhar --- for Turbulent MPFV 730,731,732,733,734 for effective dummy viscosities of phases
		case 730:
			slab=facePtr->ZSlab->DummyeffViscoeffofphase[0];
			break;
		case 731:
			slab=facePtr->ZSlab->DummyeffViscoeffofphase[1];
			break;
		case 732:
			slab=facePtr->ZSlab->DummyeffViscoeffofphase[2];
			break;
		case 733:
			slab=facePtr->ZSlab->DummyeffViscoeffofphase[3];
			break;
		case 734:
			slab=facePtr->ZSlab->DummyeffViscoeffofphase[4];
			break;
		//End of multiphase flow variables.......
		
		
		
		case 21:/**/
			slab=facePtr->ZSlab->Latent_Heat_Content;
		break;
		case 600:
			slab=facePtr->ZSlab->Scal[0];
		break;
		
		case 601:
			slab=facePtr->ZSlab->Scal[1];
		break;
		
		case 602:
			slab=facePtr->ZSlab->Scal[2];
		break;
			
		case 603:
			slab=facePtr->ZSlab->Scal[3];
		break;
		
		case 604:
			slab=facePtr->ZSlab->Scal[4];
		break;
		
		case 605:
			slab=facePtr->ZSlab->Scal[5];
		break;
		
		case 606:
			slab=facePtr->ZSlab->Scal[6];
		break;
		
		case 607:
			slab=facePtr->ZSlab->Scal[7];
		break;
		
		case 608:
			slab=facePtr->ZSlab->Scal[8];
		break;
		
		case 609:
			slab=facePtr->ZSlab->Scal[9];
		break;
		
		case 6000:
			slab=facePtr->ZSlab->ScalOld[0];
		break;
		
		case 6010:
			slab=facePtr->ZSlab->ScalOld[1];
		break;
		
		case 6020:
			slab=facePtr->ZSlab->ScalOld[2];
		break;
		
		case 6030:
			slab=facePtr->ZSlab->ScalOld[3];
		break;
		
		case 6040:
			slab=facePtr->ZSlab->ScalOld[4];
		break;
		
		case 6050:
			slab=facePtr->ZSlab->ScalOld[5];
		break;
		
		case 6060:
			slab=facePtr->ZSlab->ScalOld[6];
		break;
		
		case 6070:
			slab=facePtr->ZSlab->ScalOld[7];
		break;
		
		case 6080:
			slab=facePtr->ZSlab->ScalOld[8];
		break;
		
		case 6090:
			slab=facePtr->ZSlab->ScalOld[9];
		break;
		
		case 700:
			slab=facePtr->ZSlab->Scalcoeff[0];
		break;
		
		case 701:
			slab=facePtr->ZSlab->Scalcoeff[1];
		break;
		
		case 702:
			slab=facePtr->ZSlab->Scalcoeff[2];
		break;
		
		case 703:
			slab=facePtr->ZSlab->Scalcoeff[3];
		break;
		
		case 704:
			slab=facePtr->ZSlab->Scalcoeff[4];
		break;
		
		case 705:
			slab=facePtr->ZSlab->Scalcoeff[5];
		break;
		
		case 706:
			slab=facePtr->ZSlab->Scalcoeff[6];
		break;
		
		case 707:
			slab=facePtr->ZSlab->Scalcoeff[7];
		break;
		
		case 708:
			slab=facePtr->ZSlab->Scalcoeff[8];
		break;
		
		case 709:
			slab=facePtr->ZSlab->Scalcoeff[9];
		break;
		
		case 710:
			slab=facePtr->ZSlab->Rho;
		break;
		
		case 711:
			slab=facePtr->ZSlab->Viscoeff;
		break;
		
		case 712:
			slab=facePtr->ZSlab->Diffcoeff;
		break;
		
		case 713:
			slab=facePtr->ZSlab->SpecificHeat;
		break;
		case 714:
			slab=facePtr->ZSlab->EffViscoeff;
		break;
		case 715:
			slab=facePtr->ZSlab->TViscoeff;
		break;
                case 716:
			slab=facePtr->ZSlab->EffScalcoeff;
		break;
		case 717:
			slab=facePtr->ZSlab->ChargeDiffcoeff;
		break;
		case 718:
			slab=facePtr->ZSlab->TurbConductivity;
		break;
		case 719:
			slab=facePtr->ZSlab->UnSatPermeability;
		break;
		
		//$$$$$$$$$$$$$$$$$$ for scaler in groundwater  Grad in X direction$$$$$$$$$$$$$$$$$$$$//
		
		case 60000:
			slab=facePtr->ZSlab->GradScalX[0];
		break;
		case 60001:
			slab=facePtr->ZSlab->GradScalX[1];
		break;
		
		case 60002:
			slab=facePtr->ZSlab->GradScalX[2];
		break;
		case 60003:
			slab=facePtr->ZSlab->GradScalX[3];
		break;
		
		case 60004:
			slab=facePtr->ZSlab->GradScalX[4];
		break;
		
		case 60005:
			slab=facePtr->ZSlab->GradScalX[5];
		break;
		
		case 60006:
			slab=facePtr->ZSlab->GradScalX[6];
		break;
		
		case 60007:
			slab=facePtr->ZSlab->GradScalX[7];
		break;
		
		case 60008:
			slab=facePtr->ZSlab->GradScalX[8];
		break;
		
		case 60009:
			slab=facePtr->ZSlab->GradScalX[9];
		break;
		
		//$$$$$$$$$$$$$$$$$$ for scaler in groundwater  Grad in Y direction$$$$$$$$$$$$$$$$$$$$//
				
		case 60100:
			slab=facePtr->ZSlab->GradScalY[0];
		break;
			
		case 60101:
			slab=facePtr->ZSlab->GradScalY[1];
		break;
		
		case 60102:
			slab=facePtr->ZSlab->GradScalY[2];
		break;
		case 60103:
			slab=facePtr->ZSlab->GradScalY[3];
		break;
		
		case 60104:
			slab=facePtr->ZSlab->GradScalY[4];
		break;
		
		case 60105:
			slab=facePtr->ZSlab->GradScalY[5];
		break;
		
		case 60106:
			slab=facePtr->ZSlab->GradScalY[6];
		break;
		
		case 60107:
			slab=facePtr->ZSlab->GradScalY[7];
		break;
		
		case 60108:
			slab=facePtr->ZSlab->GradScalY[8];
		break;
		
		case 60109:
			slab=facePtr->ZSlab->GradScalY[9];
		break;
		
		//$$$$$$$$$$$$$$$$$$ for scaler in groundwater  Grad in Z direction$$$$$$$$$$$$$$$$$$$$//
		case 60200:
			slab=facePtr->ZSlab->GradScalZ[0];
		break;
			
		case 60201:
			slab=facePtr->ZSlab->GradScalZ[1];
		break;
		
		case 60202:
			slab=facePtr->ZSlab->GradScalZ[2];
		break;
		case 60203:
			slab=facePtr->ZSlab->GradScalZ[3];
		break;
		
		case 60204:
			slab=facePtr->ZSlab->GradScalZ[4];
		break;
		
		case 60205:
			slab=facePtr->ZSlab->GradScalZ[5];
		break;
		
		case 60206:
			slab=facePtr->ZSlab->GradScalZ[6];
		break;
		
		case 60207:
			slab=facePtr->ZSlab->GradScalZ[7];
		break;
		
		case 60208:
			slab=facePtr->ZSlab->GradScalZ[8];
		break;
		
		case 60209:
			slab=facePtr->ZSlab->GradScalZ[9];
		break;
		
	//$$$$$$$$$$$$$$$$$$ END FOR GROUNDWATER $$$$$$$$$$$$$$$$$$$$//
		case 990:
			slab=facePtr->ZSlab->Flow_Potential;
		break;
		//By adhi+shash :For Disperse phase flow variables.............
		//changed by Mrunalini
		
		case 62:
			slab=facePtr->ZSlab->U_Disperse;
		break;
		
		case 620:
			slab=facePtr->ZSlab->U_DisperseOld;
		break;

		case 6200:
			 slab=facePtr->ZSlab->U_DisperseOld;
		break;

		case 62000:
			slab=facePtr->ZSlab->U_DmassOld;
		break;

		case 63:
			slab=facePtr->ZSlab->V_Disperse;
		break;
		
		case 630:
			slab=facePtr->ZSlab->V_DisperseOld;
		break;

		case 6300:
			 slab=facePtr->ZSlab->V_Dmass;
		break;

		case 63000:
			slab=facePtr->ZSlab->V_DmassOld;
		break;

		case 64:
			slab=facePtr->ZSlab->W_Disperse;
		break;

		case 640:
			slab=facePtr->ZSlab->W_DisperseOld;
		break;

		case 6400:
			  slab=facePtr->ZSlab->W_Dmass;
		break;	

		case 64000:
			slab=facePtr->ZSlab->W_DmassOld;
		break;

		case 65:
			slab=facePtr->ZSlab->VarthetaDisperse;
		break;
		
		case 650:
			slab=facePtr->ZSlab->VarthetaDisperseOld;
		break;
		
		case 66:
			slab=facePtr->ZSlab->RhoDisperse;
		break;
		/*case 67:
			slab=facePtr->ZSlab->NumDensity;//shashwat
		break;
Commented as we no longer solve for number density*/
		case 68:
			slab=facePtr->ZSlab->DropFuelMassFrac;
		break;
		case 680:
			slab=facePtr->ZSlab->DropFuelMassFracOld;
		break;
		case 69:
			slab=facePtr->ZSlab->Tdrop;
		break;
		case 690:
			slab=facePtr->ZSlab->TdropOld;
		break;
		case 725:
			slab=facePtr->ZSlab->CPdrop;
		break;
		case 726:
			slab=facePtr->ZSlab->DummyViscoeff_disperse;//chandrasekhar--Turbulent gas-particle flow(Dummy viscosity of the disperse phase)
		break;
		case 727:
			slab=facePtr->ZSlab->DiffVolfrac_disperse;//chandrasekhar--Turbulent gas-particle flow(Dummy viscosity of the disperse phase)
		break;
		//End of Disperse phase  flow variables.......
		default:
			printf("Vartag is = %d Proper variable tag is not supplied to FindSlabPtrBC and exiting\n",vartag);
				exit(1);
		break;
		
	}
	return(slab);
}

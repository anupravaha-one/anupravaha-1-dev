#include "APKMakeDiffusionCoeffDomain.h"
/**
The function MakeDiffusionCoeffDomain() takes the variable tag and pointer corresponding to a zone boundary
as input and returns the corresponding Diffusion Coefficient pointer in the domain.
*/
double *MakeDiffusionCoeffDomain(Zone *zone,int vartag)
{
	double *domain;

	switch(vartag)
	{
		case 1: 
		case 2: 
		case 3: 
		case 100: 
		case 200: 
		case 300: 
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->Viscoeff;
			else
				domain=zone->Flowsol->Tcoeff->EffViscoeff;
			break;
		
		case 10: 
		case 20: 
		case 30: 
			domain=zone->PrevFlowsol->Viscoeff;
			break;

		case 4:
			if((zone->FlagMPFV==1))
				domain=zone->Flowsol->Rho;
			else
			domain=zone->Flowsol->PressureCoeff;
		break;

		/*break;

		case 5:
			domain=zone->ZGeom->Volume;
		break;*/

		case 6:
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->Diffcoeff;
			else
				domain=zone->Flowsol->Tcoeff->EffTurbConductivity;
		break;

		case 60:
			domain=zone->PrevFlowsol->Diffcoeff;
		break;

		case 7:
                	if(zone->Flowsol->flagTurbComb==1)
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
                        else
			domain=zone->Flowsol->Diffcoeff;
		break;

		case 70:
                	if(zone->Flowsol->flagTurbComb==1)
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
                        else
			domain=zone->PrevFlowsol->Diffcoeff;
		break;
		
		//Electric Parameters
		case 8:
			
			domain=(double *)malloc(sizeof(double));
			*(domain)=1.0;
		break;
		case 80:
			
			domain=(double *)malloc(sizeof(double));
			*(domain)=1.0;
		break;
		
		case 9:
			domain=zone->Flowsol->ChargeDiffcoeff;
		break;
		
		case 90:
			domain=zone->PrevFlowsol->ChargeDiffcoeff;
		break;
		
		//End of Electrical Parameters

                case 14:
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
		break;

                case 140:
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
		break;

                case 15:
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
		break;

                case 150:
			domain=zone->Flowsol->Tcoeff->EffScalcoeff;
		break;
		
		case 17:
			domain=zone->Flowsol->PotentialMHDCoeff;
		break;
		case 170://added by rajendra (not by naveen)
			domain=zone->Flowsol->PotentialMHDCoeff;
		break;
		case 56: //rajendra
			//if(zone->Flowsol->flagTurb==0)
			domain=zone->Flowsol->Magdiff;
// 			else
// 			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;

		case 560: //rajendra
			//if(zone->Flowsol->flagTurb==0)
			domain=zone->Flowsol->Magdiff;
// 			else
// 			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;

		case 57: //rajendra
			//if(zone->Flowsol->flagTurb==0)
			domain=zone->Flowsol->Magdiff;
// 			else
// 			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;

		case 570: //rajendra
			//if(zone->Flowsol->flagTurb==0)
			domain=zone->Flowsol->Magdiff;
// 			else
// 			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;

		case 58: //rajendra
			//if(zone->Flowsol->flagTurb==0)
			domain=zone->Flowsol->Magdiff;
// 			else
// 			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;

		case 580: //rajendra
			//if(zone->Flowsol->flagTurb==0)
			domain=zone->Flowsol->Magdiff;
// 			else
// 			domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;
		
		case 59://rajendra
			domain=zone->Flowsol->PsiCoeff;
		break;
		case 590://rajendra
			domain=zone->Flowsol->PsiCoeff;
		break;
		case 19:	//Groundwater Flow
		 	if(zone->Flowsol-> flagTestSaturation==1)
		 		domain=zone->Flowsol->UnSatPermeability;
		 	else
				domain=zone->Flowsol->HeadCoeff;
		break;
		
		case 600:
			domain=zone->Flowsol->Scalcoeff[0];
		break;
		
		case 6000:
			domain=zone->PrevFlowsol->Scalcoeff[0];
		break;
		
		case 601:
			domain=zone->Flowsol->Scalcoeff[1];
		break;
		
		case 6010:
			domain=zone->PrevFlowsol->Scalcoeff[1];
		break;
		
		case 602:
			domain=zone->Flowsol->Scalcoeff[2];
		break;
		
		case 6020:
			domain=zone->PrevFlowsol->Scalcoeff[2];
		break;	
			
		case 603:
			domain=zone->Flowsol->Scalcoeff[3];
		break;
		
		case 6030:
			domain=zone->PrevFlowsol->Scalcoeff[3];
		break;	
		
		case 604:
			domain=zone->Flowsol->Scalcoeff[4];
		break;
		
		case 6040:
			domain=zone->PrevFlowsol->Scalcoeff[4];
		break;	
		
		case 605:
			domain=zone->Flowsol->Scalcoeff[5];
		break;
		
		case 6050:
			domain=zone->PrevFlowsol->Scalcoeff[5];
		break;	
		
		case 606:
			domain=zone->Flowsol->Scalcoeff[6];
		break;
		
		case 6060:
			domain=zone->PrevFlowsol->Scalcoeff[6];
		break;	
		
		case 607:
			domain=zone->Flowsol->Scalcoeff[7];
		break;
		
		case 6070:
			domain=zone->PrevFlowsol->Scalcoeff[7];
		break;	
		
		case 608:
			domain=zone->Flowsol->Scalcoeff[8];
		break;
		
		case 6080:
			domain=zone->PrevFlowsol->Scalcoeff[8];
		break;	
		
		case 609:
			domain=zone->Flowsol->Scalcoeff[9];
		break;
		
		case 6090:
			domain=zone->PrevFlowsol->Scalcoeff[9];
		break;
		
		
		//Added by sushil for multiphase flow
		case 31:
		case 32:
		case 33: 
		case 3100:
		case 3200:
		case 3300:	
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->Viscoeff;
			
			else
				domain=zone->Flowsol->Tcoeff->EffViscoeff;
		break;
		case 310:
		case 320:
		case 330:
		 
			domain=zone->PrevFlowsol->Viscoeff;
		break;
		
		case 41:
		case 42:
		case 43: 
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->ViscoeffDummy[0];
			else
				domain=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[0];
		break;
		
		case 410:
		case 420:
		case 430: 
			domain=zone->PrevFlowsol->ViscoeffDummy[0];
		break;
		
		case 44:
		case 45:
		case 46: 
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->ViscoeffDummy[1];
			else
				domain=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[1];
		break;
		
		case 440:
		case 450:
		case 460: 
			domain=zone->PrevFlowsol->ViscoeffDummy[1];
		break;
		
		case 47:
		case 48:
		case 49: 
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->ViscoeffDummy[2];
			else
				domain=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[2];
		break;
		
		case 470:
		case 480:
		case 490: 
			domain=zone->PrevFlowsol->ViscoeffDummy[2];
		break;
		
		case 50:
		case 51:
		case 52: 
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->ViscoeffDummy[3];
			else
				domain=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[3];
		break;
		
		case 500:
		case 510:
		case 520: 
			domain=zone->PrevFlowsol->ViscoeffDummy[3];
		break;
		
		case 53:
		case 54:
		case 55: 
			if(zone->Flowsol->flagTurb==0)
				domain=zone->Flowsol->ViscoeffDummy[4];
			else
				domain=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[4];
		break;
		
		case 530:
		case 540:
		case 550: 
			domain=zone->PrevFlowsol->ViscoeffDummy[4];
		break;
		
		case 34:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac[0];
		break;
		case 35:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac[1];
		break;
		case 36:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac[2];
		break;
		case 37:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac[3];
		break;
		case 38:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac[4];
		break;
		case 68://dropfuelmassfrac for evaporating droplet//shash 
		case 680: //this is temporarily equal to viscoeff assuming lewis number unity.
			domain=zone->Flowsol->DropDiffcoeff;
		break;
		case 62:
		case 63:
		case 64:
		case 6200:
		case 6300:
		case 6400:		
			domain=zone->Flowsol->Tcoeff->DummyViscoeff_disperse;
		break;

		case 65:
			domain=zone->Flowsol->Tcoeff->DiffVolfrac_disperse;
		break;
		//end
		default:
			printf("\nDiffusion Coeff for given vartag(%d) not found,exiting from MakeDiffusionCoeffDomain.c\n",vartag);
			exit(1);
		break;
	}
	return (domain);
}

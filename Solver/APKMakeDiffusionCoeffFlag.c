#include "APKMakeDiffusionCoeffFlag.h"
/**
The function MakeDiffusionCoeffFlag() takes the variable tag and pointer corresponding to a zone boundary
as input and returns the flag indicating whether the diffusion coefficient is a constant
or a variable.
*/
int MakeDiffusionCoeffFlag(Zone *zone,int vartag)
{
	int flag;

	switch(vartag)
	{
		case 1:/**/
		case 2:/**/
		case 3:/**/
		case 10:/**/
		case 20:/**/
		case 30:/**/
		case 100:/**/
		case 200:/**/
		case 300:/**/
			if(zone->Flowsol->flagTurb==0)
				flag=zone->Flowsol->flagViscoeff;
			else
				flag = 1;
		break;
			
		case 4:
			if(zone->Flowsol->Momentumx!=NULL)
				flag = 1;
			else
			flag=0;
		break;
		case 40:
			flag=0;
		break;
		
		/*case 5:
			flag=zone->ZGeom->Volume;
		break;*/
		
		case 6:
		case 60:
			if(zone->Flowsol->flagTurb==0)
				flag=zone->Flowsol->flagDiffcoeff;
			else
				flag = 1;
		break;
		
		case 7://Shyam
                	if(zone->Flowsol->flagTurbComb==0)
			flag=zone->Flowsol->flagDiffcoeff;
                        else
                        flag = 1;
		break;

		case 70:
                 	if(zone->Flowsol->flagTurbComb==0)
			flag=zone->Flowsol->flagDiffcoeff;
                        else
                        flag = 1;
		break;
		
		//Electric Parameters
		case 8:
			flag=0;
		break;
		
		case 80:
			flag=0;
		break;
		
		case 9:
			flag=zone->Flowsol->flagChargeDiffCoeff;
		break;
		
		case 90:
			flag=zone->Flowsol->flagChargeDiffCoeff;
		break;
		
		//End of Electrical Parameters

                case 14:
                	flag = 1; //In  Tur. Combustion, effective trub. diffusivity is a varying one.

                break;

                case 140:
                	flag = 1;
                break;

                case 15:
                	flag = 1;
                break;

                case 150:
                	flag = 1;
                break;
		case 56://rajendra
			//if(zone->Flowsol->flagTurb==0)
			flag=zone->Flowsol->flagMagdiff;
// 			else
// 			flag = 1;
		break;
		case 560://rajendra
			//if(zone->Flowsol->flagTurb==0)
			flag=zone->Flowsol->flagMagdiff;
// 			else
// 			flag = 1;
		break;
		case 57://rajendra
			//if(zone->Flowsol->flagTurb==0)
			flag=zone->Flowsol->flagMagdiff;
// 			else
// 			flag = 1;
		break;
		case 570://rajendra
			//if(zone->Flowsol->flagTurb==0)
			flag=zone->Flowsol->flagMagdiff;
// 			else
// 			flag = 1;
		break;
		case 58://rajendra
			//if(zone->Flowsol->flagTurb==0)
			flag=zone->Flowsol->flagMagdiff;
// 			else
// 			flag = 1;
		break;
		case 580://rajendra
			//if(zone->Flowsol->flagTurb==0)
			flag=zone->Flowsol->flagMagdiff;
// 			else
// 			flag = 1;
		break;
		case 600:
			flag=zone->Flowsol->flagScalcoeff[0];
		break;
		
		case 6000:
			flag=zone->Flowsol->flagScalcoeff[0];
		break;
		
		case 601:
			flag=zone->Flowsol->flagScalcoeff[1];
		break;
		
		case 6010:
			flag=zone->Flowsol->flagScalcoeff[1];
		break;
		
		case 602:
			flag=zone->Flowsol->flagScalcoeff[2];
		break;
		
		case 6020:
			flag=zone->Flowsol->flagScalcoeff[2];
		break;	
			
		case 603:
			flag=zone->Flowsol->flagScalcoeff[3];
		break;
		
		case 6030:
			flag=zone->Flowsol->flagScalcoeff[3];
		break;	
		
		case 604:
			flag=zone->Flowsol->flagScalcoeff[4];
		break;
		
		case 6040:
			flag=zone->Flowsol->flagScalcoeff[4];
		break;	
		
		case 605:
			flag=zone->Flowsol->flagScalcoeff[5];
		break;
		
		case 6050:
			flag=zone->Flowsol->flagScalcoeff[5];
		break;	
		
		case 606:
			flag=zone->Flowsol->flagScalcoeff[6];
		break;
		
		case 6060:
			flag=zone->Flowsol->flagScalcoeff[6];
		break;	
		
		case 607:
			flag=zone->Flowsol->flagScalcoeff[7];
		break;
		
		case 6070:
			flag=zone->Flowsol->flagScalcoeff[7];
		break;	
		
		case 608:
			flag=zone->Flowsol->flagScalcoeff[8];
		break;
		
		case 6080:
			flag=zone->Flowsol->flagScalcoeff[8];
		break;	
		
		case 609:
			flag=zone->Flowsol->flagScalcoeff[9];
		break;
		
		case 6090:
			flag=zone->Flowsol->flagScalcoeff[9];
		break;
		
		//Added by sushil for multiphase flow.
		
		case 31:/**/
		case 32:
		case 33:
		case 3100:
		case 3200:
		case 3300:	
			if(zone->Flowsol->flagTurb==0)
				flag=zone->Flowsol->flagViscoeff;
			else
				flag = 1;
		break;
		case 310:/**/
		case 320:
		case 330:
			flag=zone->Flowsol->flagViscoeff;
			break;
		
		case 34:
		case 35:
		case 36:
		case 37:
		case 38:
		case 340:
		case 350:
		case 360:
		case 370:
		case 380:
			flag = 1;
		break;
			
		case 41:/**/
		case 44:
		case 47:
		case 50:
		case 53:
		
			if(zone->Flowsol->flagTurb==0)
				flag=zone->Flowsol->flagViscoeff;
			else
				flag = 1;
		break;
		case 410:/**/
		case 440:
		case 470:
		case 500:
		case 530:
		
			if(zone->Flowsol->flagTurb==0)
				flag=zone->Flowsol->flagViscoeff;
			else
				flag = 1;
		break;
		//End of multiphase flow.
		case 68://dropfuelmassfrac for evaporating droplet//shash 
		case 680: //this is temporarily equal to viscoeff assuming lewis number unity.
			flag=zone->Flowsol->flagDropDiffcoeff;
		break;
		case 62:
		case 63:
		case 64:
		case 6200:
		case 6300:
		case 6400:
			if(zone->Flowsol->flagTurb==0)
				flag=0;
			else
				flag = 1;
		break;
		case 65:
			if(zone->Flowsol->flagTurb==0)
				flag=0;
			else
				flag = 1;
		break;
		default:
			printf("\nFlag of Diffusion Coeff for given vartag not found,exiting from MakeDiffusionCoeffFlag.c\n");
			exit(1);
		break;
	}
	return (flag);
}

#include "APKfind_name.h"
/***************************************************IMPORTANT*********************************

		This function returns the name of the variable according to its 
		vartag (just to print the name while asking for the boundary conditions)
				  
*********************************************************************************************/
char *find_name(Base *base,int count)
{
	switch(base->VARTAG[count-1])
	{
		case 1:
			{
			if((base->flagParticle==1)||(base->flagDroplet==1))
			return("u velocity Carrier phase");
			else
			return("u velocity ");
			}
		break;
		case 2:
			if((base->flagParticle==1)||(base->flagDroplet==1))
			return("v velocity Carrier phase");
			else
			return("v velocity ");
		break;
		case 3:
			if((base->flagParticle==1)||(base->flagDroplet==1))
			return("w velocity Carrier phase");
			else
			return("w velocity ");
		break;
		case 4:
			return("Pressure");
		break;
		case 6:
			return("Temperature");
		break;
		case 7:
			return("enthalpy");
		break;
		
		case 8:
			return("Electric Potential");
		break;
		case 9:
			return("Ion Charge Density");
		break;
		case 11:
			return("X-dirn Efield");
		break;
		case 12:
			return("Y-dirn Efield");
		break;
		case 13:
			return("Z-dirn Efield");
		break;
		
		case 14:
			return("Mixture Fraction");
		break;
                case 15:
			return("Fuel Mass Fraction");
		break;
		case 16:
			return("Radiative Intensity");
		break;
		case 17: //naveen
			return("Electric potential MHD");
		break;
		case 56://rajendra
			return("TotalMfieldX");
		break;
		case 57://rajendra
			return("TotalMfieldY");
		break;
		case 58://rajendra
			return("TotalMfieldZ");
		break;
		case 59://rajendra
			return("Psi");
		break;
		case 19: //Groundwater
			return("Head Pressure(pressure head)");
		break;
		case 600:
			if(base->BaseParameters->Turb->numTurbScalars>=1)
			{
				if(base->flagTurbModel==8)
					return("Modified Kinematic Eddy Viscosity-nu_t^tilda");
				
				else
					return("Turbulent Kinetic Energy (k)");
			}
			return("1 Scalar");
		break;
		case 601:
			if(base->BaseParameters->Turb->numTurbScalars>=2)
			{
				if(base->flagTurbModel==2 || base->flagTurbModel==6)
					return("Specific Dissipation Rate of 'k' (Omega)");
				
				else
					return("Dissipation Rate of 'k' (Epsilon)");
			}
			return("2 Scalar");
		break;
		case 602:
			if(base->BaseParameters->Turb->numTurbScalars>=3)
			{
				return("v2");
			}
			return("3 Scalar");
		break;
		case 603:
			if(base->BaseParameters->Turb->numTurbScalars>=4)
			{
				return("Damping Scalar (f)");
			}
			return("4 Scalar");
		break;
		case 604:
			return("5 Scalar");
		break;
		case 605:
			return("6 Scalar");
		break;
		case 606:
			return("7 Scalar");
		break;
		case 607:
			return("8 Scalar");
		break;
		case 608:
			return("9 Scalar");
		break;
		case 609:
			return("10 Scalar");
		break;
		case 990:
			return("Flow_Potential");
		break;
		//Added by sushil
		
		
		case 34:
			return("0th phase_volumefraction");
		break;
		case 35:
			return("1st phase_volumefraction");
		break;
		case 36:
			return("2nd phase_volumefraction");
		break;
		case 37:
			return("3rd phase_volumefraction");
		break;
		case 38:
			return("4th phase_volumefraction");
		break;
		case 41:
			return("0th phase_Individual_Momentum_X");
		break;
		case 42:
			return("0th phase_Individual_Momentum_Y");
		break;
		case 43:
			return("0th phase_Individual_Momentum_Z");
		break;
		case 44:
			return("1st phase_Individual_Momentum_X");
		break;
		case 45:
			return("1st phase_Individual_Momentum_Y");
		break;
		case 46:
			return("1st phase_Individual_Momentum_Z");
		break;
		case 47:
			return("2nd phase_Individual_Momentum_X");
		break;
		case 48:
			return("2nd phase_Individual_Momentum_Y");
		break;
		case 49:
			return("2nd phase_Individual_Momentum_Z");
		break;
		case 50:
			return("3rd phase_Individual_Momentum_X");
		break;
		case 51:
			return("3rd phase_Individual_Momentum_Y");
		break;
		case 52:
			return("3rd phase_Individual_Momentum_Z");
		break;
		case 53:
			return("4th phase_Individual_Momentum_X");
		break;
		case 54:
			return("4th phase_Individual_Momentum_Y");
		break;
		case 55:
			return("4th phase_Individual_Momentum_Z");
		break;
		//added by adhi+shash for disperse phase flow
		case 62:
			return("Disperse Phase U velocity");
		break;
		case 63:
			return("Disperse Phase V velocity");
		break;
		case 64:
			return("Disperse Phase W velocity");
		break;
		case 65:
			return("Normalized Disperse phase volume fraction");
		break;
/*		case 67:
			return("Number Density");//shashwat
		break;
*/
		case 68:
			return("Drop Fuel Mass Fraction");//shashwat
		break;
		case 69:
			return("Drop Temperature");//shashwat
		break;
		case 710:
			return("Carrier phase Density");//shashwat
		break;
		
		default :
		{
			printf("\nCase not found, exiting from find_name.c\n");
			exit(1);
		}
		break;
	}
}



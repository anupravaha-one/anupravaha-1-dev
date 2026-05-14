#include "APKfind_convcri.h"

double find_convcri(int vartag,Base *base)
{

	switch(vartag)
	{
		case 1:
		case 2:
		case 3:
		case 100:
		case 200:
		case 300:
			return (*(base->BaseParameters->convergenceCriterion+0));
		break;
		case 4:
			return (base->BaseParameters->PressLoopConv);
		break;
		case 6:
		case 7:
                	return (*(base->BaseParameters->convergenceCriterion+1));
		break;
		case 17:
                        return (base->BaseParameters->PotentialMHDLoopConv);
		break;
		case 56://raj
		case 57:
		case 58:
                        return (base->BaseParameters->InducedMHDLoopConv);
		break;
            	case 59:
                        return (base->BaseParameters->PsiMHDLoopConv);
		break;
		case 990:
			return (*(base->BaseParameters->convergenceCriterion+1));
		break;
		
		case 8:
		case 9:
			return(1e-08);
		break;
		
		case 14:
		case 15:
			return(1e-08);
		break;
			
		case 600:
			return (*(base->BaseParameters->convergenceCriterion+2));
		break;
		case 601:
			return (*(base->BaseParameters->convergenceCriterion+3));
		break;
		case 602:
			return (*(base->BaseParameters->convergenceCriterion+4));
		break;
		case 603:
			return (*(base->BaseParameters->convergenceCriterion+5));
		break;
		case 604:
			return (*(base->BaseParameters->convergenceCriterion+6));
		break;
		case 605:
			return (*(base->BaseParameters->convergenceCriterion+7));
		break;

		case 606:
			return (*(base->BaseParameters->convergenceCriterion+8));
		break;
		case 607:
			return (*(base->BaseParameters->convergenceCriterion+9));
		break;
		case 608:
			return (*(base->BaseParameters->convergenceCriterion+10));
		break;
		case 609:
			return (*(base->BaseParameters->convergenceCriterion+11));
		break;
		case 715:
			return (1.0e-8);
		break;
		
		//Added for multiphase flow
		case 31:
		case 32:
		case 33:
		case 3100:
		case 3200:
		case 3300:
		case 41:
		case 42:
		case 43:
		case 44:
		case 45:
		case 46:
			return (*(base->BaseParameters->convergenceCriterion+0));
			
			break;
		case 34:
		case 35:
		case 36:
		case 37:
		case 38:
			return(1e-8);
			break;
		//end
		//added by adhi+shash
		case 62:
		case 63:
		case 64:
		case 6200:
		case 6300:
		case 6400:
			return (*(base->BaseParameters->convergenceCriterion+0));
			break;
		case 65:
			return (1.0e-8);
			break;
	/*	case 67:
			return (1.0e-5);
			break;	*/
		case 68:
			return (1.0e-7);
		break;
		case 69:
			return (*(base->BaseParameters->convergenceCriterion+1));
		break;
		case 710:
			return (1.0e-8);
		break;
		default:
			printf("The variableTag [%d] is not yet added in find_convcri.c\n\n", vartag);
			exit(2);
			return (1e-08);
		break;
	}
	return (1e-08);
}


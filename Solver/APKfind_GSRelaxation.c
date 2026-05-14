#include "APKfind_GSRelaxation.h"
/**
This function returns the Gauss-Siedel (over/under) relaxation parameter for a given variable
*/

double find_GSRelaxation(int vartag,Base *base)
{

	switch(vartag)
	{
		case 1:
		case 2:
		case 3:
		case 100:
		case 200:
		case 300:
			return (*(base->BaseParameters->OMEGA+0));
			break;
		case 4:
			return (*(base->BaseParameters->OMEGA+1));
			break;
		case 6:
			return (*(base->BaseParameters->OMEGA+2));
			break;
		case 990:
			return (*(base->BaseParameters->OMEGA+0));
			break;
		
		case 600:
			return (*(base->BaseParameters->OMEGA+4));
			break;
		case 601:
			return (*(base->BaseParameters->OMEGA+5));
			break;
		case 602:
			return (*(base->BaseParameters->OMEGA+6));
			break;
		case 603:
			return (*(base->BaseParameters->OMEGA+7));
			break;
		case 604:
			return (*(base->BaseParameters->OMEGA+8));
			break;
		case 605:
			return (*(base->BaseParameters->OMEGA+9));
			break;
		case 606:
			return (*(base->BaseParameters->OMEGA+10));
			break;
		case 607:
			return (*(base->BaseParameters->OMEGA+11));
			break;
		case 608:
			return (*(base->BaseParameters->OMEGA+12));
			break;
		case 609:
			return (*(base->BaseParameters->OMEGA+13));
			break;
		default:
// 			printf("The variableTag [%d] is not yet added in find_GSRelaxation.c\n\n", vartag);
// 			exit(2);
			return (*(base->BaseParameters->OMEGA+0));
			break;
	}
	return (*(base->BaseParameters->OMEGA+0));
}


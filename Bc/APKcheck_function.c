#include "APKcheck_function.h"

/**
This function checks the choice of the function index it should be an integer between 1->15 other wise this 
function will give an erroer message. 
*/
int check_function(int c)
{ 
	int test;
 	
	if(c==0||c==1)
	{
	test=1;
	}
	else
	{
	test=0;
	printf("\nINCORRECT INPUT\n");
	printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
	}

return(test);
}

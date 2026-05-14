#include "APKTerminate.h"

/* -------------------------------------------------------------- */
/*   Function for Terminating the Main if some error is found     */
/* -------------------------------------------------------------- */

void Terminate(char *msg)
{
	printf("\n Terminating the program ....\n");
        printf("%s\n",msg);
        exit(1);
}

#include"APKCGNSErrorCheck.h"


/* ------------------------------------------------------------- */
/*  Function for checking & printing any CGNS error encountered  */
/* ------------------------------------------------------------- */
/* err = 0 corresponds to NO CGNS ERROR                          */
/* ------------------------------------------------------------- */


/**
While reading or writing CGNS file if some error occurs than the following message will be displayed in 
the screen ``Following CGNS Error Encountered". Present function will not tell about the error type and line 
number in file, It will only give an error messange. 
*/
void CGNSErrorCheck(int err)
{
   if(err!=0)
   {
      printf("Following CGNS Error Encountered\n");
      printf("--------------------------------\n");
      cg_error_exit();
   }
}


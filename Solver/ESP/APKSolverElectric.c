#include "APKSolverElectric.h"

/**
This function SolverElectric() is the main function where the program execution for electric field starts and sequentially carries out the 
solution of electric field depending upon the type of solver selected for solving the electric field equations.
At present, only SemiCoupled method is implemented for solving the two non-linear
electric field equations namely potential and charge density.

Variable tags for electric potential and charge density are 8 and 9 respectively.
*/

int SolverElectric(Root CgnsRoot)
{
	char Potential[] = "SOLUTIONS/Output_Electric.cgns";
	
	
	printf("Started solvin for electric field\n\n");
	
	if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
		SolverElectricSemiCoupled(CgnsRoot);
	else
		SolverElectricSemiCoupled2d(CgnsRoot);
	
	//switch(CgnsRoot.CgnsBase->EcoupTag)
	//{
		//case 1:
			
		//break;
	//}
	
	//Calculate Electric Field along X, Y, Z Directions
	if(CgnsRoot.CgnsBase->BaseParameters->TagEfield==1)
		ElectricFieldXYZ(CgnsRoot.CgnsBase);
	
	//Write the solution
	WriteCgnsFile(Potential,CgnsRoot);
	
	printf("\n\n\t\t****************************************************************************\n\t\tElectric field solution is written into %s file\n\t\t****************************************************************************\n", Potential);

	return 0; 

}



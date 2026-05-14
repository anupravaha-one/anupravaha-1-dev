#include "APKuser_defined_funcs2.h"

/**
  Original APKuser_defined_funcs.c file was spilt into two
  As user are allowed to modify only the thirteen functions, they are put in APKuser_defined_funcs.c file
  and the remaining subroutines are kept in this file. 
 **/

/**
These are user defined functions (1-16) in which user can write functions in terms of equations for any particular 
variable. This function takes input as X,Y,Z coordinates, Zone pointer and physical-time. Thus which provides
the facility to write equations in terms of coordinates and time. As Zone pointer is provided functions can also 
be written in terms of flow and scalar variables. It returns the value calulated from the the function. 
It is basically used for giving boundary conditions in the domain. 
*/

//X,Y and Z are face centres on a boundary.

double userfunc_Energysource(double X,double Y,double Z,double time,Zone *zone)
{
	double variable=500.0;
	variable = 1;
	return(variable);
}

double CalcLagrangeCoeff(double X, double XC, double B, double C)
{
	double Coeff;
	
	Coeff = ((X-B)*(X-C))/((XC-B)*(XC-C));
	
	return(Coeff);

}


//NOTE: Don't tamper with this function ... This is used in k-Omega and SST turbulence models ...
double user_func_14(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	//Above X,Y and Z are face centres on a boundary.
	int NI,NJ,NK,flagViscoeff=0,flagRho=0,index=0;
	double omega=0.0,beta=3.0/40.0,nu,rho=1.0,y1=0.0,mu;

	printf("----------------- Applying boundary condition for Omega -----------------\n");
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	//Warning!!! Call this function for Turbulent problems only and also if the boundary is Solid Wall!
	if(zone->flagSolidFluid == 0 && zonebc->Boundary_Tag == 22)
	{
		if(zonebc->Direction == 1 || zonebc->Direction == 2)
		{
			index = (NJ-1)*k+j;
			y1 = *(zonebc->deltaP+index);
		}

		if(zonebc->Direction == 3 || zonebc->Direction == 4)
		{
			index = (NI-1)*k+i;
			y1 = *(zonebc->deltaP+index);
		}

		if(zonebc->Direction == 5 || zonebc->Direction == 6)
		{
			index = (NI-1)*j+i;
			y1 = *(zonebc->deltaP+index);
		}
	
		flagViscoeff = zone->Flowsol->flagViscoeff;
		flagRho = zone->Flowsol->flagRho;
		rho = *(zonebc->ZSlab->Rho+index*flagRho);
		mu  = *(zonebc->ZSlab->Viscoeff+index*flagViscoeff);
		nu = mu/rho;
		omega = 60.0*nu/(beta*y1*y1);
		printf("nu=%e\ty1=%e\tOmega at the solid boundary = %lf\n", nu,y1,omega);
	}
	else
	{
		if(zone->flagSolidFluid == 1)
			omega = 0.0;
		else
		{
			printf("***** Warning!!! This boundary is not tagged as BCWallViscous ***** \nThis function should not be applied to this boundary.\n");
			exit(300);
		}
	}

	return (omega);
}

//NOTE: Don't tamper with this function ... This is used in k-Omega turbulence model ...
double user_func_15(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	//Above X,Y and Z are face centres on a boundary.
	int NI,NJ,NK,flagViscoeff=0,flagRho=0,index=0;
	double omega=0.0,beta=3.0/40.0,nu,rho=1.0,y1=0.0,mu;

	printf("----------------- Applying boundary condition for Omega -----------------\n");
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	//Warning!!! Call this function for Turbulent problems only and also if the boundary is Solid Wall!
	if(zone->flagSolidFluid == 0 && zonebc->Boundary_Tag == 22)
	{
		if(zonebc->Direction == 1 || zonebc->Direction == 2)
		{
			index = (NJ-1)*k+j;
			y1 = *(zonebc->deltaP+index);
		}

		if(zonebc->Direction == 3 || zonebc->Direction == 4)
		{
			index = (NI-1)*k+i;
			y1 = *(zonebc->deltaP+index);
		}

		if(zonebc->Direction == 5 || zonebc->Direction == 6)
		{
			index = (NI-1)*j+i;
			y1 = *(zonebc->deltaP+index);
		}
	
		flagViscoeff = zone->Flowsol->flagViscoeff;
		flagRho = zone->Flowsol->flagRho;
		rho = *(zonebc->ZSlab->Rho+index*flagRho);
		mu  = *(zonebc->ZSlab->Viscoeff+index*flagViscoeff);
		nu = mu/rho;
		omega = 6.0*nu/(beta*y1*y1);
		printf("nu=%e\ty1=%e\tOmega at the solid boundary = %lf\n", nu,y1,omega);
	}
	else
	{
		if(zone->flagSolidFluid == 1)
			omega = 0.0;
		else
		{
			printf("***** Warning!!! This boundary is not tagged as BCWallViscous ***** \nThis function should not be applied to this boundary.\n");
			exit(300);
		}
	}

	return (omega);
}

double user_func_16(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)// for saturated groundwater test case 
{
	double variable;

	variable = (1-((cosh(1.2*(Y)))/cosh(1.20/2)))/(1-((tanh(1.20/2))/(1.20/2)));
	return (variable);
}





	//double Tx,Ty,Tz;
	/*Tx=(X+Y+Z)/sqrt(3.0);
	Ty=(Y-Z)/sqrt(2.0);
	Tz=(2.0*X-Y-Z)/sqrt(6.0);
	variable=2*Ty+3*Tz+4;*/
	//variable = 150.0*(Y-0.1)*(0.3-Y);
	//variable = 6.0*Y*(1.0-Y);

	//	variable=10.0+0.5*X;  
// 	printf("in uiser_fun %e",X);   

	//variable = 150.0*(Y-0.1)*(0.3-Y);// Do not delete this <--- This function is used for Conjuagte Test case 3
	//	variable = 16*(Y-0.25)*(0.75-Y);// Do not delete this <--- This function is used for Backward facing Step Test Case
	//variable =  5.931528 * (1.0 - Y*Y);
	//variable = -600.0*Y*(0.1-Y);// Do not delete this <--- This function is used for Conjuagte Test case 1 and 2
	//return(variable);
	//variable = 6.0*Y*(6.0-Y)/36.0;
	/*	if(time<1.1)
	variable = 1.0;
	else
	variable = 2.0;
	/*FILE *Charge;
	
	double variable, Rmax, Rmin, Rad, X0, X1, Xm, c1, c2, c3;
	//The following is for charge density
	Charge = fopen("Settings/ESPInput/InputCharge.dat", "r");
	fscanf(Charge,"%lf", &Rmax);
	fscanf(Charge,"%lf", &Rmin);
	fscanf(Charge,"%lf", &Rad);
	fscanf(Charge,"%lf", &Xm);
	printf("\n*********************************************************\n");
	printf("\nRmax = %e\tRmin = %e\tRadius = %12.8lf\n", Rmax, Rmin, Rad);
	printf("\n*********************************************************\n");
	fclose(Charge);*/
	
	//Do not delete the above code for test case
	
	//Xm is the distance from origin to mid of wire in x-dirn
	
	//Xm = 0.105; X0 = Xm - Rad; X1 = Xm + Rad;//Huang
	//Xm = 0.07635;
	//X0 = Xm - Rad; X1 = Xm + Rad;
	
	//Calculate Lagrangian coefficients
	
	//c1 = CalcLagrangeCoeff(X, X0, X1, Xm);
	//c2 = CalcLagrangeCoeff(X, Xm, X0, X1);
	//c3 = CalcLagrangeCoeff(X, X1, X0, Xm);
	
	//variable = c1*Rmin + c2*Rmax + c3*Rmin;
	//printf("\nVariable = %e\n", variable);
	
	
	//variable = Rmax*(Y*Y)/(Rad*Rad)+Rmin*(1-(Y*Y)/(Rad*Rad));
	




	//double Tx,Ty,Tz;
	//Tx=(X+Y+Z)/sqrt(3.0);
	//Ty=(Y-Z)/sqrt(2.0);
	//Tz=(2.0*X-Y-Z)/sqrt(6.0);
	//	variable=1.5*Y*(2.0-Y);// Do not delete this <--- This function is used for Conjuagte Test case 1 and 2
	/*if(Y<0.2142 && Y>0.7858)
		variable = pow(0.0294,1.5) * pow(0.09,0.75)/(0.42 * Y);
	else
		variable = pow(0.0294,1.5) * pow(0.09,0.75)/(0.09);*/
	//variable = 6.08662253 * (1.0 - (Y*Y/16.0) + 0.6762633 * logl(Y/4.0));
	//variable = 0.0001039/(Y+1);
	//	variable = (4.0/(1-(tanh(25))/25))*(1-(cosh(100*(Y)))/cosh(25));
	//	variable = (1-((cosh(100*(Y)))/cosh(25)))/((1-(1.0/cosh(25))));
	//	variable = (1-((cosh(2.5*(Y)))/cosh(2.50/2)))/(1-((tanh(2.5/2))/(2.5/2)));
	//	if(X<=5)
	//  variable = 1./(1+exp(-(X-2)/0.15));
	  //	if(X>5)




	//variable=150.0*(Y+0.1)*(0.1-Y);

	
	//variable=4.0*(Y)*(1-Y);//Do not delete for ESP test case
       	//variable = 1.5*(Y-0.25)*(2.25-Y);

	//	variable = (1-((cosh(10*(Y)))/cosh(0.25)))/((1-(1.0/cosh(0.25))));
	//	variable = (1-((cosh(3.5*(Y)))/cosh(3.5/2)))/(1-((tanh(3.5/2))/(3.5/2)));
	//variable=-600.0*(Y)*(0.1-Y);

	//double variable,half,MeanDistance;

	
	//variable = (1-((cosh(5*(Y)))/cosh(5.0/2)))/(1-((tanh(5.0/2))/(5.0/2)));


	
	//	variable=7*X;
	//variable=0.5*(X+Y+Z);
	//variable=6*(Y)*(1.0-Y);
	//	variable=10*Z;
	//	variable = 1./(1+exp(-(X-2)/0.15));  
	//variable = (1-((cosh(6.0*(Y)))/cosh(6.0/2)))/(1-((tanh(6.0/2))/(6.0/2)));
	//	variable=11*X;
	//	variable = 1./(1+exp(-(Z-4)/0.15));  
	//variable = (1-((cosh(7.0*(Y)))/cosh(7.0/2)))/(1-((tanh(7.0/2))/(7.0/2)));

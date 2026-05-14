#include "APKuser_defined_funcs.h"

/**
These are user defined functions (1-13) in which user can write functions in terms of equations for any particular 
variable. This function takes input as X,Y,Z coordinates, Zone pointer and physical-time. Thus which provides
the facility to write equations in terms of coordinates and time. As Zone pointer is provided functions can also 
be written in terms of flow and scalar variables. It returns the value calulated from the the function. 
It is basically used for giving boundary conditions in the domain. 
*/

/**

X,Y and Z are face centres on a boundary.
User is allowed to define and use as many local variables as needed. 
But the function should return a single numerical value for given X, Y, Z coordinate value

User are not allowed to define new subroutines or change the subroutine names

**/


double user_func_1(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable = 6 * Y * (1-Y);		// Used for Laminar Backward Facing Step Tutorial
	return(variable);
}

double user_func_2(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable =  0.0;
	return(variable);
}

double user_func_3(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable = 0.0;
	return (variable);
}

double user_func_4(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable = 0.0;
	return (variable);
}

double user_func_5(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
  	double variable;

  	variable = 0.0;
	return (variable);
}

double user_func_6(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;
	
  	variable = 0.0;
	return (variable);
}

double user_func_7(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable = 0.0;
	return (variable);
}

double user_func_8(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{

	double variable;

	variable = 0.0;
	return (variable);
}
	
double user_func_9(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{

	double variable;

	variable = 0.0;
	return (variable);
}

double user_func_10(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable = 0.0;
	return (variable);
}	

double user_func_11(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;

	variable = 0.0;
	return (variable);
}	

double user_func_12(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double U;

	U = 0.0;
	return(U);
}

double user_func_13(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc)
{
	double variable;
	
	variable = 0.0;
	return(variable);
}	


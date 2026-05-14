#ifndef __APKData_Structures_h
#define __APKData_Structures_h

/* ******************************************************************************************************* *
 *  !!!!! IMPORTANT !!!!!! : The original Data Structure of version 3.2 is a standard and is 
 *                           recognised(Downloaded from CGNS website).CHANGES have been made to this
 *                           Data Structure to meet our needs and named as Data_Structures4.1.h 
 * ---------------------------------------------------------------------------------------------------------
 * This header file contains all the structure definitions(Version NO. 3.2).
 * Each element of all the structures has a tag enclosed in square brackets which identifies the type
 * of subroutine or file which determines that element. The Values of tag are as follows:
 *	+ C - CGNS File
 *	+ A - API Subroutines/Functions
 *	+ S - Solver Subroutines/Functions
 * ******************************************************************************************************* */

/*----------------Data Structure for Boundary conditions----------*/
typedef struct
{
	///stores the indices of the points to be monitored
	int *i,*j,*k;
	///stores the zone number of the point to be monitored
	int *zone_num;
	int number;
}monitors;

typedef struct
{
	/// Brief: It stores the name of the variable like u, v, w etc. 
	char Name_of_variable[50];   	/* [A]  Name_of_variable                                     */
	/// Brief: It stores the information about the type of the boundary condition like Dirichlet, Neumann etc.
        int  Bc_Type;                	/* [A]  Type of boundary condition(Dirichlet,Neumann etc.)   */
	/// Brief: It stores the information about the boundary condition to be a function or constant. There are four elements for implementation of Robin boundary condition. 0->Constant 1->Function
	int  identity[4];   		/* [A]  Refer Bc_Patch.c for details                         */
	/// Brief: It store   for time dependent boundary condition. 0->No time dependent 1->time dependent
	int  flagTime ;	
	/// Brief: It stores the Index of function for implementing of time dependent boundary conditions
	int  FunctionIndex[3];		/*      Time dependent boundarycondition		     */
	/// Brief: This is a pointer type which stores the exact value of boundary conditions at the face center of each cell. 
	double *Data_Array;		/* [A]  Refer Bc_Patch.c for details                         */
	/// Brief: If the boundary condition is a function, then it is stored here in funct. It is a pointer. 
	double *funct;               	/* [A]  Refer Bc_Patch.c for details                         */
	/// Brief: If the boundary condition is a function and to implement Robin boundary condition, then it is stored here in funct_beta. It is a pointer. 
	double *funct_beta;          	/* [A]  Refer Bc_Patch.c for details                         */
	/// Brief: If the boundary condition is a function and to implement Robin boundary condition, than it is stored here in funct_gama. It is a pointer. 
	double *funct_gama;          	/* [A]  Refer Bc_Patch.c for details                         */
}Dataset;
/** \struct Dataset
	This is Dataset Structure. It contains all information about the boundary condition implementation. 
	It is pointer type which stores information about boundary condition as a function or constant etc. 
*/

//----------------------------------- BaldwinLomax Structure ------------------------------------------
typedef struct
{
	///Brief: Pointer to Mixing length.
	double *MixingLength;
	///Brief: Pointers to turbulence viscosity in inner and outer layers respectively.
	double *mu_ti,*mu_to;
} BaldwinLomax;
/** \struct BaldwinLomax
	This is BaldwinLomax structure. It contains variables used in the Baldwin-Lomax turbulence model.  This structure is called in the "TCoeff" and "ZoneSlab" nodes of the data structure.
*/

/* **-----------------------------------Zone Slab Information------------------------------------------** */
typedef struct
{
	///Stores the target indices of the neighbouring zone incase of Interface  This points to NULL incase of Boundary.
	int *Target_Indices;
	///Stores the number of cells at the interface/boundary.
	int Total_Cells;
	
	/// Brief: It stores the Volume of boundary and interface slabs which stores Zero for boundary slabs. 
	double	*Volume;
	/// Brief: It stores the Density of the boundary or interface slab. 
	double   *Rho;			/* [S] Array for Density                */
	/// Brief: It stores the dynamic viscocity of the boundary or interface slab. 
	double   *Viscoeff;		/* [S] Array for Viscosity              */
	/// Brief: It stores the Effective viscosity of the boundary or interface slab. 
	double   *EffViscoeff;		/* [S] Array for EffectiveViscosity(Dynamic) */
	//This stores the effective diffusivity of turb. combustion scalars
        double   *EffScalcoeff; //Shyam /* [S], This variable stores effective diffusivity of combustion variables*/
	/// Brief: It stores the Turbulent viscosity of the boundary or interface slab.
	double   *TViscoeff;		/* [S] Array for TurbulentViscosity(Dynamic) */
	/// Brief: It stores the turbulent thermal conductivity of the boundary or interface slab. 
	double *TurbConductivity;
	///Brief:it stores the effective turbulent thermal conductivity of boundary or slab.
	double *EffTurbConductivity;
	/// Brief: It stores the thermal conductivity of the boundary or interface slab. 
	double   *Diffcoeff;		/* [S] Array for Thermal Diffusivity    */
	/// Brief: It stores the specific heat of the boundary or interface slab. 
	double   *SpecificHeat;		
	/// Brief: It stores the current time step U velocity in the boundary or interface slab. 
	double   *u;			/* [S] Array for X-Velocity             */
	/// Brief: It stores the Previous time step U velocity in the boundary or interface slab.
	double   *uOld;			/* [S] Array for X-Velocity             */
	/// Brief: It stores the current time step V velocity in the boundary or interface slab. 
	double   *v;			/* [S] Array for Y-Velocity             */
	/// Brief: It stores the Previous time step V velocity in the boundary or interface slab. 
	double   *vOld;
	/// Brief: It stores the current time step W velocity in the boundary or interface slab. 
	double   *w;			/* [S] Array for Z-Velocity 	        */
	/// Brief: It stores the Previous time step W velocity in the boundary or interface slab. 
	double   *wOld;
	/// Brief: It stores the current time step U mass velocity in the boundary or interface slab. 
	double   *umass;		
	/// Brief: It stores the Previous time step U mass  velocity in the boundary or interface slab. 
	double   *umassOld;
	/// Brief: It stores the current time step V mass  velocity in the boundary or interface slab. 
	double   *vmass;		
	/// Brief: It stores the Previous time step V mass  velocity in the boundary or interface slab. 
	double   *vmassOld;
	/// Brief: It stores the current time step W mass  velocity in the boundary or interface slab. 
	double   *wmass;		
	/// Brief: It stores the Previous time step W mass  velocity in the boundary or interface slab. 
	double   *wmassOld;
	/// Brief: It stores the Current time step Pressure in the boundary or interface slab. 
	double   *P;			/* [S] Array for Pressure               */
	/// Brief: It stores the Current time step Temperature in the boundary or interface slab. 
	double   *T;			/* [S] Array for Temperature            */
	/// Brief: It stores the Previous time step Temperature in the boundary or interface slab. 
	double   *TOld;			/* [S] Array for Temperature            */
	
	/// Brief: It stores the Current time step enthalpy values in the boundary or interface slab.
	double   *enthalpy   ;//Created by bhans It stores current time step enthalpy values
	/// Brief: It stores the Current time step liquid fraction values in the boundary or interface slab.
	double   *f          ;//Created by bhans It stores the liquid fraction values in the slabs
	/// Brief: It stores the Previous time step enthalpy values in boundary or interface slab
	double   *enthalpyOld;//Created by bhans It stores the Old enthalpy values for implementing Orlanski.
	/// Brief: It stores the "Latent heat content" of boundary and interface slab.
	double   *Latent_Heat_Content;//Created by bhans It stores Latent heat content of each cell
	/// Brief: It stores the Current time step 10 Scalar in the boundary or interface slab. 
	double   *Scal[10];		/* [S] Array for scalars                */
	/// Brief: It stores the Previous time step 10 Scalar in the boundary or interface slab.
	double   *ScalOld[10];		/* [S] Array for scalars                */
 	/// Brief: It stores the Current time step 10 Scalar Diffusion coefficent in the boundary or interface slab.
	double   *Scalcoeff[10];	/* [S] Array for diff coeff for scalars */
	/// Brief: It stores the normal temperature gradient in the boundary or interface slab. 
	double   *GradT;
	/// Brief: It stores the normal Scalar gradient in the boundary or interface slab. 
	double   *GradScal[10];
	/// Brief: This pointer stores the minimum distance from a slab face center to any arbitrary solid wall ... used in SST turbulence model
	double *minWallDistance;
	
	//Combustion Parameters
	///Brief: It stores current time step values of mixture fraction in slabs(Interface and Boundary)
        double   *MixFrac; //created by shyam
	///Brief: It stores previous time step values of mixture fraction in slabs(Interface and Boundary)
	double   *MixFracOld;  //created by shyam
	///Brief: It stores current time step values of fuel mass fraction in slabs(Interface and Boundary)
	double   *FuelMassFrac; //created by shyam
	///Brief: It stores previous time step values of fuel mass fraction in slabs(Interface and Boundary)
	double   *FuelMassFracOld; //created by shyam
	///Brief: It stores current time step values of oxygen mass fraction in slabs(Interface and Boundary)
        double	 *OxiMassFrac; //Shyam
	///Brief: It stores current time step values of inert gas mass fraction in slabs(Interface and Boundary)	
        double   *InMassFrac[10];  //Shyam
	///Brief: It stores current time step values of product mass fraction in slabs(Interface and Boundary)	
        double   *ProdMassFrac; //Shyam
	///Brief: It stores previous time step values of oxygen mass fraction in slabs(Interface and Boundary)	
        double	 *OxiMassFracOld; //Shyam
	///Brief: It stores previous time step values of inert gas mass fraction in slabs(Interface and Boundary)	
        double   *InMassFracOld[10];  //Shyam
	///Brief: It stores previous time step values of product mass fraction in slabs(Interface and Boundary)	
        double   *ProdMassFracOld; //Shyam
	///Brief: It stores net molecular weight.
        double   *MolWteq;  //Shyam
	///Brief: It stores  CPs of individual species
        double   *Cps;      //Shyam	
	
	// The following are the ELECTRICAL parameters
	
	/// Brief: It stores the current time step electric Potential in the boundary or interface slab. 
	double   *Potential;		/* [S] Array for Electric Potential     */
	/// Brief: It stores the current time step X-direction electric field in the boundary or interface slab. 
	double   *EfieldX;		/* [S] Array for X-Dirn Electric Field  */
	/// Brief: It stores the current time step Y-direction electric field in the boundary or interface slab.
	double   *EfieldY;		/* [S] Array for Y-Dirn Electric Field  */
	/// Brief: It stores the current time step Z-direction electric field in the boundary or interface slab.
	double   *EfieldZ;		/* [S] Array for Z-Dirn Electric Field  */
	/// Brief: It stores the current time step Ion Charge Density in the boundary or interface slab.
	double   *RhoIon;		/* [S] Array for Space Charge Density   */
	/// Brief: It stores the Previous time step electric potential in the boundary or interface slab.
	double   *PotentialOld;		/* [S] Array for Electric Potential     */ 
	/// Brief: It stores the Previous time step ion charge density in the boundary or interface slab. 
	double   *RhoIonOld;		/* [S] Array for Space Charge Density   */
	/// Brief: It stores the dielectric constant of the boundary or interface slab.
	double   *DiElecConst;		/* [S] Array for Di-Electric Constant   */
	/// Brief: It stores the ion mobility rate of the boundary or interface slab.
	double   *MobIon;		/* [S] Array for Ion Mobility Rate      */
	/// Brief: It stores the charge diffusion coefficient of the boundary or interface slab.
	double 	 *ChargeDiffcoeff;	/* [S] Stores the value of Charge Diffusion Coefficient   */
	
	// The following are the LMMHD parameters (Implemented by naveen)
	
	/// Brief: It stores the current time step electric Potential in the boundary or interface slab. 
	double   *PotentialMHD;		/* [S] Array for Electric Potential     */
	/// Brief: It stores the previous time step electric Potential in the boundary or interface slab. 
	double   *PotentialMHDold;	/* [S] Array for Electric Potential     */
	/// Brief: It stores the current time step X-direction Magnetic field in the boundary or interface slab. 
	double   *MfieldX;		/* [S] Array for X-Dirn Magnetic Field  */
	/// Brief: It stores the current time step Y-direction Magneticc field in the boundary or interface slab.
	double   *MfieldY;		/* [S] Array for Y-Dirn Magnetic Field  */
	/// Brief: It stores the current time step Z-direction Magnetic field in the boundary or interface slab.
	double   *MfieldZ;		/* [S] Array for Z-Dirn Magnetic Field  */
	/// Brief: It stores the Electrical conductivity of the liquid metal for boundary or interface slab. 
	double   *Eleconductivity;	/* [S] Array for Electrical conductivity */   
	double   *Permeability;//rajendra
	double   *UcrossBX;
	double   *UcrossBY;
	double   *UcrossBZ;
	double   *TotalMfieldX;//rajendra
	double   *TotalMfieldXold;
	double   *TotalMfieldY;
	double   *TotalMfieldYold;
	double   *TotalMfieldZ;
	double   *TotalMfieldZold;
	double   *Psi;
	double   *Psiold;
	double   *Magdiff;//rajendra
	
// 	Groundwater flow
	///Brief: used in ground water..
	///Brief: It store the pressure head in the current time step
	double  *H;
	///Brief: It store the Saturated hydraulic conductivity in the current time step
	double  *SatHydConductivity;
	///Brief: It store the Water content(theta) in the current time step
	double  *WaterContent;
	///Brief: It store the UnSaturated permeabilty in the current time step
	double  *UnSatPermeability;
	///Brief: It store the Gradient of scaler in x , y, z direction in the current time step
	double  *GradScalX[10];
	double  *GradScalY[10];
	double  *GradScalZ[10];
	///Brief: It will store value of Deltainverse(geometry parameters),in domain 
	double  *Delta;
	
	double *Flow_Potential;
	
	
	//	Added by sushil ....These are the multiphase variables...............
	
	///Brief:It stores current time step value of overall  momentum in x direction for boundary or interface slab.
	double  *Momentumx;
	///Brief:It stores current time step value of overall  momentum star in x direction for boundary or interface slab.
	double *Momentumxstar;
	///Brief:It stores previous time step value of overall  momentum in x direction for boundary or interface slab
	double  *MomentumxOld;
	///Brief:It stores previous time step value of overall  momentum star in x direction for boundary or interface slab
	double *MomentumxstarOld;
	///Brief:It stores current time step value of overall  momentum in y direction for boundary or interface slab.
	double  *Momentumy;
	///Brief:It stores current time step value of overall  momentum star in y direction for boundary or interface slab.
	double *Momentumystar;
	///Brief:It stores previous time step value of overall  momentum in y direction for boundary or interface slab
	double *MomentumyOld;
	///Brief:It stores previous time step value of overall  momentum star in y direction for boundary or interface slab
	double *MomentumystarOld;
	///Brief:It stores current time step value of overall  momentum in z direction for boundary or interface slab.
	double  *Momentumz;
	///Brief:It stores current time step value of overall  momentum star in z direction for boundary or interface slab.
	double *Momentumzstar;
	///Brief:It stores previous time step value of overall  momentum in z direction for boundary or interface slab
	double *MomentumzOld;
	///Brief:It stores previous time step value of overall  momentum star in z direction for boundary or interface slab
	double *MomentumzstarOld;
	///Brief:It stores current time step value of individual phase  momentum in x direction for boundary or interface slab.
	double  *Momentumofphasex[5];
	///Brief:It stores previous time step value of individual phase  momentum in x direction for boundary or interface slab.
	double  *MomentumofphasexOld[5];
	///Brief:It stores current time step value of  individual phase momentum in y direction for boundary or interface slab.
	double  *Momentumofphasey[5];
	///Brief:It stores previous time step value of  individual phase momentum in y direction for boundary or interface slab.
	double  *MomentumofphaseyOld[5];
	///Brief:It stores current time step value of individual phase momentum in z direction for boundary or interface slab.
	double  *Momentumofphasez[5];
	///Brief:It stores previous time step value of  individual phase momentum in z direction for boundary or interface slab.
	double  *MomentumofphasezOld[5];
	///Brief:It stores density of individual phase for boundary or interface slab.
	double   *RhoofPhase[5];
	///Brief:It stores viscosity of individual phase for boundary or interface slab.
	double  *ViscosityofPhase[5];
	///Brief:It stores each phase viscosity*volume fraction
	double  *ViscoeffDummy[5];
	
	///Brief:It stores current time step value of  individual phase velocity in x direction for boundary or interface slab.
	double  *U_phase[5];
	///Brief:It stores current time step value of  individual phase velocity in y direction for boundary or interface slab.
	double  *V_phase[5];
	///Brief:It stores current time step value of  individual phase velocity in z direction for boundary or interface slab.
	double  *W_phase[5];
	///Brief:It stores current time step value of  individual phase volume fraction for boundary or interface slab.
	double  *VarthetaPhase[5];
	
	///Brief:It stores previous time step value of  individual phase velocity in x direction for boundary or interface slab.
	double  *U_phaseOld[5];
	///Brief:It stores previous time step value of  individual phase velocity in y direction for boundary or interface slab.
	double  *V_phaseOld[5];
	///Brief:It stores previous time step value of  individual phase velocity in z direction for boundary or interface slab.
	double  *W_phaseOld[5];
	///Brief:It stores previous time step value of  individual phase volume fraction for boundary or interface slab.
	double  *VarthetaPhaseOld[5];
	///Brief: It holds the flag for Density to be constant or variable for phases. 0->constant Density; 1->variable Density
	int      flagRhoPhase[5];		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: It holds the flag for Dynamic Viscocity to be constant or variable for phases. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeffPhase[5];
	
	///  Intensity for radiative Problem
	double *RadIntensity;
	/// Refractive Index
	double *RefIndex;
	///Emissivity of Surface
	double epsilon;
	
	//Added by Adhiraj for dispersed phase flow...
	///Brief:It stores current time step value of dispersed phase  momentum in x direction for boundary or interface slab.
	double  *MomentumDispersex;
	///Brief:It stores previous time step value of dispersed phase  momentum in x direction for boundary or interface slab.
	double  *MomentumDispersexOld;
	///Brief:It stores current time step value of  dispersed phase momentum in y direction for boundary or interface slab.
	double  *MomentumDispersey; 
	///Brief:It stores previous time step value of  dispersed phase momentum in y direction for boundary or interface slab.
	double  *MomentumDisperseyOld;
	///Brief:It stores current time step value of dispersed phase momentum in z direction for boundary or interface slab.
	double  *MomentumDispersez;
	///Brief:It stores previous time step value of  dispersed phase momentum in z direction for boundary or interface slab.
	double  *MomentumDispersezOld;

//Mrunalini	
	///Brief:It stores current time step value of  dispersed phase velocity in x direction for boundary or interface slab.
	double  *U_Disperse;
	///Brief:It stores current time step value of  dispersed phase velocity in y direction for boundary or interface slab.
	double  *V_Disperse;
	///Brief:It stores current time step value of  dispersed phase velocity in z direction for boundary or interface slab.
	double  *W_Disperse;

	///Brief:It stores current time step value of  dispersed phase velocity in x direction for boundary or interface slab.
	double  *U_DisperseOld;
	///Brief:It stores current time step value of  dispersed phase velocity in y direction for boundary or interface slab.
	double  *V_DisperseOld;
	///Brief:It stores current time step value of  dispersed phase velocity in z direction for boundary or interface slab.
	double  *W_DisperseOld;
	
	///Brief:It stores current time step value of  dispersed phase mass velocity in x direction for boundary or interface slab.
	double  *U_Dmass;
	///Brief:It stores current time step value of  dispersed phase mass velocity in y direction for boundary or interface slab.
	double  *V_Dmass;
	///Brief:It stores current time step value of  dispersed phase mass velocity in z direction for boundary or interface slab.
	double  *W_Dmass;
	
	///Brief:It stores previous time step value of  dispersed phase mass velocity in x direction for boundary or interface slab.
	double  *U_DmassOld;
	///Brief:It stores previous time step value of  dispersed phase mass velocity in y direction for boundary or interface slab.
	double  *V_DmassOld;
	///Brief:It stores previous time step value of  dispersed phase mass velocity in z direction for boundary or interface slab.
	double  *W_DmassOld;

	///Brief:It stores density of dispersed  phase for boundary or interface slab.
	double   *RhoDisperse;
	///Brief:It stores viscosity of dispersed phase for boundary or interface slab.
	double  *ViscosityDisperse;
	
	///Brief:It stores viscosity of dispersed phase for boundary or interface slab.
	double  *ViscoeffDummyDisperse;
	
	///Brief:It stores current time step value of  dispersed phase volume fraction for boundary or interface slab.
	double  *VarthetaDisperse;
	///Brief:It stores previous time step value of  dispersed phase volume fraction for boundary or interface slab.
	double  *VarthetaDisperseOld;
	///Brief:It stores current  time step value of  dispersed phase number density for droplet phase for boundary or interface slab.
	double  *NumDensity;//shashwat
	
	
	//For Heat-Up period of evoparation
	double *Tdrop;
	
	double *TdropOld;
	
	double   *CPdrop;		
	//for evaporation droplet phase flow
	
	double *DropFuelMassFrac;
	
	double *DropFuelMassFracOld;
	/// Brief: It stores the mass diffusivity of the boundary or interface slab. 
	double   *DropDiffcoeff;		/* [S] Array for mass diffusivity for evaporated fuel mass fraction    */	
	//End of variables addded for Dispersed phase
	///Pointer to the parameters used in Baldwin-Lomax turbulence Model
	BaldwinLomax *BL;
	
	//Added by chandrasekhar --- for Turbulent MPFV
	/// Brief: It stores the Turbulent viscosity of individual phase at the boundary or interface slab.
	double  *TViscoeffofphase[5];//mu_tq(0),mu_tq(1).....
	/// Brief: It stores the Effective viscosity of individual phase at the boundary or interface slab. 
	double  *EffViscoeffofphase[5];//mu_effq(0),mu_effq(1)..........
	/// Brief: It stores the Dummyeffective viscosity(product of effective viscosity of phase and averaged volume fraction of phase) of individual phase at the boundary or interface slab, i.e. diffusion coefficient of the phasic momentum equation in case of Turbulent MPFV.
	double  *DummyeffViscoeffofphase[5];//mu_eff(0)*chi(0),mu_eff(1)*chi(1)..diffusion coefficient in the phasic momentum equation.
	///Brief: It stores the Diffusion coefficient for Volume fraction in case of turbulent MPFV.
	double *DiffVolfrac[5];//nu_tq/sigma_q..diffusion coefficient for volume fraction equation in case of Turbulent MPFV.
	//End --- chandrasekhar
	
	//for Turbulent gas-particle flow--chandrasekhar
	/// Brief: It is a pointer to an array which stores kinematic turbulent viscosity of disperse phase in case of Turbulent gas-particle flow.
	double   *TViscoeff_disperse;
	/// Brief: It is a pointer to an array which stores  diffusion coefficient for particle phase momentum equations in case of Turbulent gas-particle flow.
	double   *DummyViscoeff_disperse;
	/// Brief: It is a pointer to an array which stores diffusion coefficient of volumefraction equation in case of Turbulent Gas-Particle flow.
	double   *DiffVolfrac_disperse;
	/// Brief: It is a pointer to an array which stores eddy characteristic time T_e in case of Turbulent Gas-Particle flow.
	double   *T_e;
	

} ZoneSlab;
/** \struct ZoneSlab
	This is ZoneSlab Structure. It contains all information about the Zone Slab flow variable. 
	It is a pointer which has all Zone boundary or interface slab data, like flow variable, various coefficients, 
	Variable gradient etc. 
*/


/* ** ---------------------------- Data Structure for TurbCoeff Node ------------------------------ ** */

typedef struct
{
	/// Brief: It is a pointer to an array which stores turbulent viscosity.
	double   *TViscoeff;
	/// Brief: It is a pointer to a constant in Turbulent kinetic energy equation.
	double   *TKprandtl;
	/// Brief: It is a pointer to a constant in Turbulent Dissipation of energy equation.
	double   *TEprandtl;
	/// Brief: It is a pointer to an array which stores effective viscosity, which works as diffusion coefficient for RANS.
	double   *EffViscoeff;
        ///Brief: This stores the effective diffusivity of turb. combustion scalars
        double   *EffScalcoeff; //Shyam
	/// Brief: It is a pointer to an array which stores the production term for K and Epsilon equation.
	double   *TProduction;
	/// Brief: It is a pointer to a constant in Turbulent energy dissipation equation for RNG k-epsilon model.
	double   *C0epsi;
	/// Brief: It is a pointer to a constant in Turbulent energy dissipation equation.
	double   *C1epsi;
	/// Brief: It is a pointer to a constant in Turbulent energy Dissipation equation.
	double   *C2epsi;
	/// Brief: It is a pointer to a constant in Turbulent energy dissipation equation.
	double   *f1;
	/// Brief: It is a constant used in wall function treatment, having value as 0.42.
	double   Chi;
	/// Brief: It is a constant used in Turbulent kinetic energy equation, having value as 1.0.
	double   sigmaK;
	/// Brief: It is a constant used in the equation for dissipation of turbulent kinetic energy, having value as 1.39.
	double   sigmaEpsilon;
	/// Brief: It is a constant used in the equations governing turbulent fluid flow, having value as 0.09.
	double   Cmu;
	/// Brief: It is a constant used in volume fraction equation in case of Turbulent MPFV.
	double   sigmaVolfrac;
	/// Brief: It is a constant used  in case of Turbulent Gas-Particle, having value of 5.0
	double   sigma_f;
	/// Brief: It is a constant used in volume fraction equation in case of Turbulent MPFv, having value of 10^10
	double   sigma_r;
	//Constants for k-Omega Model
	/// Brief: It is a constant used in the equation for specific dissipation rate, having value as 0.5.
	double	sigmaOmega;
	/// Brief: It is a constant used in the equation for specific dissipation rate, having value as 3/40.
	double	beta;	//<----------- Coefficient of Dissipation term in Omega equation
	/// Brief: It is a constant used in the equation for specific dissipation rate, having value as 5/9.
	double	gamma;	//<----------- Coefficient of Production term in Omega equation
	/// Brief: It is a constant used in the equation for specific dissipation rate, having value as 0.09.
	double	betastar;	//<-------- Cmu in k-Omega Model and SST Model also
 	/// Brief: It is a constant used in the equation for specific dissipation rate, having value as 1.0.
	double	gammastar;

	// diffusion terms for Launder and Sharma model
	/// Brief: It is dissipation term in K equation for Launder and Sharma model.
	double  *PhiK;
	/// Brief: It is source term in Epsilon equation for Launder and Sharma model.
	double  *PhiEps;

	//Constants for SST turbulence Model
	double SigmaK1,SigmaK2;
	double SigmaOmega1,SigmaOmega2;
	double Beta1,Beta2;
	double Gamma1,Gamma2;
	double a1;
	//Parameters used in SST turbulence Model
	double *F1;

        // Parameters for Turbulent Combustion
	/// Brief: It is a model constant used in mean fuel mass fraction equation, having a value of 0.7.
        double   sigmacombvar; //shyam
        /// Brief: To store net molecular weight of mixture
	double   *MolWteq; //shyam
        /// Brief: It Stores Cps of individual species
        double   *Cps; //Shyam,

	//for turbulent thermal energy equation
	double Prt;//Turbulent Prandtl number
	double *TurbConductivity;
	double *EffTurbConductivity;

	///Pointer to the parameters used in Baldwin-Lomax turbulence Model
	BaldwinLomax *BL;
	
	//Added by chandrasekhar --- For turbulent MPFV 
	/// Brief: It is a pointer to an array which stores turbulent viscosity of individual phase.
	double  *TViscoeffofphase[5];//mu_tq(0),mu_tq(1).....
	/// Brief: It is a pointer to an array which stores effective viscosity of individual phase.
	double  *EffViscoeffofphase[5];//mu_effq(0),mu_effq(1)..........
	/// Brief: It is a pointer to an array which stores Dummy effective viscosity(product of effective viscosity of phase and averaged volume fraction of phase) of individual phase, i.e. diffusion coefficient of the phasic momentum equation in case of Turbulent MPFV.
	double  *DummyeffViscoeffofphase[5];//mu_eff(0)chi(0),mu_eff(1)chi(1)..diffusion coefficient in the phasic momentum equation.
	///Brief: It stores the Diffusion coefficient for Volume fraction in case of turbulent MPFV.
	double *DiffVolfrac[5];//nu_tq/sigma_q..diffusion coefficient for volume fraction equation in case of Turbulent MPFV.
	//End --- chandrasekhar
	//for Turbulent gas-particle flow--chandrasekhar
	/// Brief: It is a pointer to an array which stores kinematic turbulent viscosity of disperse phase in case of Turbulent gas-particle flow.
	double   *TViscoeff_disperse;
	/// Brief: It is a pointer to an array which stores diffusion coefficient for particle phase momentum equations in case of Turbulent gas-particle flow.
	double   *DummyViscoeff_disperse;
	/// Brief: It is a pointer to an array which stores diffusion coefficient of volumefraction equation in case of Turbulent Gas-Particle flow.
	double   *DiffVolfrac_disperse;
	/// Brief: It is a pointer to an array which stores eddy characteristic time T_e in case of Turbulent Gas-Particle flow.
	double   *T_e;
	
	
}TurbCoeff;

/* ** ---------------------------- Data Structure for Flow Solution Node ------------------------------ ** */
typedef struct
{
	/// Brief: It is a pointer to Turbulent Coefficients structure. It stores the information relating to turbulent flow quatities.
	TurbCoeff  *Tcoeff;
	/// Brief: It store the Current flow solution name.
	char   FlowSolutionName[33];	/* [C] Name of the Flow Solution Node   */
	/// Brief: `GridLocation' specifies the location of the solution data with respect to grid. For all possible values of GridLocation, see 'cgnslib.h' (line 85-93). For detailed description of GridLocation read the Developer's Manual.
	int    GridLocation;		/* [C] Location of the Grid (Read NOTE) */
	/// Brief: It is pointer type which stores the Density of Current time step at cell center of each cell in a Zone.
	double   *Rho;			/* [S] Array for Density                */
	/// Brief: It is pointer type which stores the Dynamic Viscocity of Current time step at cell center of each cell in a Zone.
	double   *Viscoeff;		/* [S] Array for Viscosity (Dynamic)    */
	/// Brief: It is pointer type which stores the Thermal Conductivity of Current time step at cell center of each cell in a Zone.
	double   *Diffcoeff;		/* [S] Array for Thermal Diffusivity    */
	/// Brief: It is pointer type which stores the Specific Heat at constant pressure of Current time step at cell center of each cell in a Zone.
	double	 *SpecificHeat;
	/// Brief: It is a array of pointer which stores all the 10 Sclar Diffusion Coefficent of Current time step at cell center of each cell in a Zone.
	double   *Scalcoeff[10];	/* [S] Array for diff coeff for scalars */
	/// Brief: It is a pointer which always stores 1. It doesen't have any physical sense but is used only in programming point of view.
	double   *PressureCoeff;
	/// Brief: It is pointer type which stores the U velocity of Current time step at cell center of each cell in a Zone.
	double   *u;			/* [S] Array for X-Velocity             */
	/// Brief: It is pointer type which stores the V velocity of Current time step at cell center of each cell in a Zone.
	double   *v;			/* [S] Array for Y-Velocity             */
	/// Brief: It is pointer type which stores the W Velocity of Current time step at cell center of each cell in a Zone.
	double   *w;			/* [S] Array for Z-Velocity 	        */
	/// Brief: It is pointer type which stores the Mass Flux of Current time step at each face of each cell in a Zone.
	double   *Massflux;		/* [S] Array for Massflux 	        */
	/// Brief: It is pointer type which stores the U mass velocity of Current time step at cell center of each cell in a Zone.
	double   *umass;
	/// Brief: It is pointer type which stores the V mass velocity of Current time step at cell center of each cell in a Zone.
	double   *vmass;
	/// Brief: It is pointer type which stores the W mass velocity of Current time step at cell center of each cell in a Zone.
	double   *wmass;
	/// Brief: It is pointer type which stores the Pressure of Current time step at cell center of each cell in a Zone.
	double   *P;			/* [S] Array for Pressure               */
	/// Brief: It is pointer type which stores the Temperature of Current time step at cell center of each cell in a Zone.
	double   *T;			/* [S] Array for Temperature            */
	/// Brief: It is a array of pointer which stores all the 10 Scalar of Current time step at cell center of each cell in a Zone.
	double   *Scal[10];		/* [S] Array for scalars                */
	/// Brief: It is pointer which stores the absolute value of Vorticity of Current time step at cell center of each cell in a Zone in SST turbulence model.
	double   *Vorticity;
	/// Brief: It is pointer type which stores the Vorticity in X direction of Current time step at cell center of each cell in a Zone.
	double   *vorticity_X;
	/// Brief: It is pointer type which stores the Vorticity in Y direction Current time step at cell center of each cell in a Zone.
	double   *vorticity_Y;
	/// Brief: It is pointer type which stores the Vorticity in Z direction Current time step at cell center of each cell in a Zone.
	double   *vorticity_Z;
	/// Brief: It is pointer type which stores the U Flux (1/2(convection+diffusion)) of the previous time step at cell center of each cell in a Zone.
	double   *PrevTStepFluxes_u;
	/// Brief: It is pointer type which stores the V Flux (1/2(convection+diffusion)) of the previous time step at cell center of each cell in a Zone.
	double   *PrevTStepFluxes_v;
	/// Brief: It is pointer type which stores the W Flux (1/2(convection+diffusion)) of the previous time step at cell center of each cell in a Zone.
	double   *PrevTStepFluxes_w;
	/// Brief: It is pointer type which stores the T Flux (1/2(convection+diffusion)) of the previous time step at cell center of each cell in a Zone.
	double   *PrevTStepFluxes_T;
	/// Brief: It is a array of pointer which stores all the 10 Scalar Flux (1/2(convection+diffusion)) of the previous time step at cell center of each cell in a Zone.
	double   *PrevTStepFluxes_Scal[10];
	/// Brief: It is pointer type which stores the enthalpy of Current time step at cell center of each cell in a Zone.
	double   *enthalpy; 		//Created by bhans
	/// Brief: It is pointer type which stores the liquid fraction of Current time step at cell center of each cell in a Zone.
	double   *f       ; 		//Created by bhans
	/// Brief: It is pointer type which stores the PrevTStepFluxes_e of Current time step at cell center of each cell in a Zone.
	double   *PrevTStepFluxes_e;	//Created by bhans
	/// Brief: It is pointer type which stores the Latent heat content of Current time step at cell center of each cell in a Zone.
	double   *Latent_Heat_Content;  //Created by bhans
	/// Brief: It is a pointer to an array which stores the value of heat generation source terms.
	double   *Energysource;
	/// Brief: It holds the flag for Turbulent flow. 1->Turbulent Flow, 0->Laminar Flow.
	int     flagTurb;
	/// Brief: It holds the tag to select turbulent combustion problem. 1->Turbulent combustion problem, 0->Not a combustion problem.
	int     flagTurbComb; // Shyam,
	/// Brief: It holds the flag for Turbulent flow model. 0-> Standard K-Epsilon, 1->RNG K-Epsilon.
	int 	flagRNG_KEPS;
	/// Brief: It holds the flag for Density to be constant or variable in a Zone. 0->constant Density; 1->variable Density
	int      flagRho;		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: /// Brief: It holds the flag for Dynamic Viscocity to be constant or variable in Zone. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeff;		/* [S] Flag storing 0->constant visccoeff; 1->variable visccoeff */
	/// Brief: /// Brief: It holds the flag for Thermal Conductivity to be constant or variable in Zone. 0->constant Thermal Conductivity; 1->variable Thermal Conductivity
	int      flagDiffcoeff;		/* [S] Flag storing 0->constant diffcoeff; 1->variable diffcoeff */
	/// Brief: /// Brief: It holds the flag for Specific heat to be constant or variable in Zone. 0->constant Specific Heat; 1->variable Specific Heat
	int      flagSpecificHeat;
	/// Brief: /// Brief: It holds the flag for scalar diffusion coefficient to be constant or variable in Zone. 0->constant Rho; 1->variable Rho
	int      flagScalcoeff[10];
	/// Brief: /// Brief: This stores (Scal[1]/Scal[2]) in k-epsilon, k-Omega, SST turbulence models.
	double *Gamma;
	/// Brief: This stores Cross Diffusion term in SST turbulence model.
	double *CrossDiff;

        //Combustion Variables	
	///Brief: It stores current time step values of mixture fraction in zones		
	double   *MixFrac;  //shyam
	///Brief: It stores current time step values of fuel mass fraction in zones			
	double   *FuelMassFrac;  //shyam	
	///Brief: It stores current time step values of oxygen mass fraction in zones				
        double	 *OxiMassFrac; //Shyam
	///Brief: It stores current time step values of inert mass fraction in zones				
        double   *InMassFrac[10];  //Shyam
	///Brief: It stores current time step values of product mass fraction in zones				
        double   *ProdMassFrac; //Shyam
	
	
	// The following are the ELECTRICAL parameters
	
	/// Brief: It is pointer type which stores the electric potential of Current time step at cell center of each cell in a Zone. 
	double	 *Potential;		/* [S] Array for Electric Potential     */
	/// Brief: It is pointer type which stores the X-direction electric field of Current time step at cell center of each cell in a Zone. 
	double   *EfieldX;		/* [S] Array for X-Dirn Electric Field  */
	/// Brief: It is pointer type which stores the Y-direction electric field of Current time step at cell center of each cell in a Zone. 
	double   *EfieldY;		/* [S] Array for Y-Dirn Electric Field  */
	/// Brief: It is pointer type which stores the Z-direction electric field of Current time step at cell center of each cell in a Zone. 
	double   *EfieldZ;		/* [S] Array for Z-Dirn Electric Field  */
	/// Brief: It is pointer type which stores the electric flux of Current time step at cell center of each cell in a Zone. 
	double   *Electricflux;		/* [S] Array for Electric Field Dlux    */
	/// Brief: It is pointer type which stores the Ion Charge Density of Current time step at cell center of each cell in a Zone. 
	double   *RhoIon;		/* [S] Array for Space Charge Density   */
	/// Brief: It is pointer type which stores the dielectric constant of Current time step at cell center of each cell in a Zone. 
	double   *DiElecConst;		/* [S] Array for Di-Electric Constant   */
	/// Brief: It is pointer type which stores the Ion Mobility rate of Current time step at cell center of each cell in a Zone. 
	double   *MobIon;		/* [S] Array for Ion Mobility Rate      */
	/// Brief: It is pointer type which stores the charge diffusion coefficient of Current time step at cell center of each cell in a Zone. 
	double 	 *ChargeDiffcoeff;	/* [S] Stores the value of Charge Diffusion Coefficient          */
	/// Brief: It holds the flag for Permittivity/Dielectric Coefficient to be constant or variable in Zone. 0->constant ; 1->variable
	int      flagDiElec;		/* [S] Flag storing 0->constant ; 1->variable  */
	/// Brief: It holds the flag for Ion mobility rate to be constant or variable in Zone. 0->constant ; 1->variable
	int      flagMobIon;		/* [S] Flag storing 0->constant ; 1->variable  */
	/// Brief: It holds the flag for Charge diffusion Coefficient to be constant or variable in Zone. 0->constant ; 1->variable
	int      flagChargeDiffCoeff;		/* [S] Flag storing 0->constant ; 1->variable  */
	
	double *ESP_Source[3];///Brief: This term stores the velocity source terms in ESP problems.  This is extremely useful incase of one-way coupling.
	
	
	
	// The following are the LMMHD parameters (Implemented by naveen)
	
	/// Brief: It is pointer type which stores the electric potential of Current time step at cell center of each cell in a Zone. 
	double	 *PotentialMHD;		/* [S] Array for Electric Potential     */
	/// Brief: It is a pointer which always stores 1. It doesen't have any physical sense but is used only in programming point of view.
	double   *PotentialMHDCoeff;
	/// Brief: It is pointer type which stores the X-direction Magnetic field of Current time step at cell center of each cell in a Zone. 
	double   *MfieldX;		/* [S] Array for X-Dirn Magnetic Field  */
	/// Brief: It is pointer type which stores the Y-direction Magnetic field of Current time step at cell center of each cell in a Zone. 
	double   *MfieldY;		/* [S] Array for Y-Dirn Magnetic Field  */
	/// Brief: It is pointer type which stores the Z-direction Magnetic field of Current time step at cell center of each cell in a Zone. 
	double   *MfieldZ;		/* [S] Array for Z-Dirn Magnetic Field  */	
	/// Brief: It is pointer type which stores the term in the RHS of the PotentialMHD poisson equation for Current time step at each face of each cell in a Zone.
	double   *PotentialMHDpoissonRHS; 
	/// Brief: It stores the Electrical conductivity of the liquid metal of Current time step at cell center of each cell in a Zone.  
	double   *Eleconductivity;	/* [S] Array for Electrical conductivity*/


	/// Brief: It holds the flag for Electrical Conductivity to be constant or variable in Zone. 0->constant Thermal Conductivity; 1->variable Thermal Conductivity
	int      flageleccond;		/* [S] Flag storing 0->constant electrical conductivity; 1->variable electrical conductivity*/

	
	/// Brief: It stores the X component of Electrical current density at the at cell center of each cell in a Zone.  
	double   *ElecurrentdensityX;	/* [S] Array for Electric current density*/
	/// Brief: It stores the Y component of Electrical current density at the at cell center of each cell in a Zone.  
	double   *ElecurrentdensityY;	/* [S] Array for Electric current density*/
	/// Brief: It stores the Z component of Electrical current density at the at cell center of each cell in a Zone.  
	double   *ElecurrentdensityZ;	/* [S] Array for Electric current density*/
	double   *UcrossBX;
	double   *UcrossBY;
	double   *UcrossBZ;
	double   *TotalMfieldX;//rajendra
	double   *TotalMfieldY;
	double   *TotalMfieldZ;
	double   *Permeability;
	double   *Magneticflux;//rajendra
	double   *Psi;
	double   *PsiCoeff;
	double   *Magdiff;
	int   flagpermeability,flagMagdiff;//rajendra

	double *LMMHD_Source[3];///Brief: This term stores the velocity source terms in LMMHD problems.  This is extremely useful incase of one-way coupling.
	double *LMMHDmass_Source[3];///Brief: This term stores the velocity source terms in LMMHD problems.  This is extremely useful incase of one-way coupling.
	///Brief: it stores the values for potential flow solution

	double *Flow_Potential;
	
	//Groundwater flow problems
	
	///Brief: used in ground water..
	///Brief: It stores the pressure head in the current time step
	double  *H;
	///Brief: It store the Saturated hydraulic conductivity in the current time step
	double  *SatHydConductivity;
	///Brief: It store the Water content(theta) in the current time step
	double  *WaterContent;
	///Brief: It store the UnSaturated permeabilty in the current time step
	double  *UnSatPermeability;
	/// Brief: It store flag for saturated and unsaturated flow 
	int flagTestSaturation;
	///Brief: as per analogus to Pressure coefficent
	double  *HeadCoeff;    
	///Brief: It store the gradient of theta(watercontent)w.r.t H(pressurehead)
	double  *GradTheta;     //
	///Brief: It store the storativity of the medium 
	double  Storativity;
	///Brief: It store the Gradient of scaler  in x, y, and z 
	double  *GradScalX[10];
	double  *GradScalY[10];
	double  *GradScalZ[10];
	///Brief: It store the Gradient of scaler  in x, y, and z 
	double  *GradU;
	double  *GradV;
	double  *GradW;
	
double *utemp;
double *vtemp;
double *wtemp;	
	///Brief: These variables are created to store the slab values of the variables directly in to them 
	///rather than trying to acess the values each and every time during calculation
	double *SlabUmass[6];
	double *SlabVmass[6];
	double *SlabWmass[6];
	double *SlabPress[6];
	double *SlabU[6];
	double *SlabV[6];
	double *SlabW[6];
	double *SlabT[6];
	double *SlabVol[6];
	double *Slab_Latent_Heat_Content[6];
	
	double *SlabMixFrac[6];
	double *SlabFuelFrac[6];
	double *SlabScalar[10][6];
	double *SlabScalar_old[10][6];
	double *prevIterVariable;
	double *VarDensityCoeff;

	double *Slab_Flow_Potential[6];
	double *Slab_Potential[6];
	double *Slab_RhoIon[6];

	double *Slab_Potential_MHD[6];
	double *Slab_Total_MfieldX[6];//rajendra
	double *Slab_Total_MfieldY[6];
	double *Slab_Total_MfieldZ[6];
	double *Slab_Psi_MHD[6];
	double *Slab_Enthalpy[6];
	
	double *SlabDensity[6];
	double *Slab_Eff_Visc[6];///For Turbulence effective viscosity.
	double *Slab_Tvisc[6];///For Temperature equations in Turbulence.
	double *Slab_Visc[6];
	double *Slab_K[6];
	double *Slab_CP[6];
	double *Slab_EffScal[6];
	double *Slab_TurbCond[6];
	
	double *Slab_H[6];
	double *SlabU_old[6];
	double *SlabV_old[6];
	double *SlabW_old[6];
	
	// Added by sushil......These are the multiphase variables parameter.
	double *SlabVartheta[5][6];
	double *SlabVxstar[6];
	double *SlabVystar[6];
	double *SlabVzstar[6];
	double *SlabVx[6];
	double *SlabVy[6];
	double *SlabVz[6];
	double *SlabVxPhase[5][6];
	double *SlabVyPhase[5][6];
	double *SlabVzPhase[5][6];
	
	double *Source_term;
	double *Diag_Source;
	double *Pressure_Grad[3];
	
	/// Brief: It is pointer type which stores the overall momentum in x direction for Current time step at cell center of each cell in a Zone.
	double  *Momentumx;
	/// Brief: It is pointer type which stores the overall momentum star values in x direction for Current time step at cell center of each cell in a Zone.
	double *Momentumxstar;
		
	/// Brief: It is pointer type which stores the overall momentum in y direction for Current time step at cell center of each cell in a Zone.
	double  *Momentumy;
	/// Brief: It is pointer type which stores the overall momentum star values in y direction for Current time step at cell center of each cell in a Zone.
	double *Momentumystar;
		
	///Brief:It stores current time step value of overall  momentum in z direction for boundary or interface slab.
	double  *Momentumz;
	/// Brief: It is pointer type which stores the overall momentum star values in z direction for Current time step at cell center of each cell in a Zone.
	double *Momentumzstar;
	/// Brief: It is pointer type which stores the individual momentum of each phase in x direction for Current time step at cell center of each cell in a Zone.
	double  *Momentumofphasex[5];
	
	/// Brief: It is pointer type which stores the individual momentum of each phase in y direction for Current time step at cell center of each cell in a Zone.
	double  *Momentumofphasey[5];
	
	//// Brief: It is pointer type which stores the individual momentum of each phase in z direction for Current time step at cell center of each cell in a Zone.
	double  *Momentumofphasez[5];
	
	/// Brief: It is pointer type which stores the Density(phasic density) of Current time step at cell center of each cell in a Zone.
	double   *RhoofPhase[5];
	/// Brief: It is pointer type which stores the viscosity(phasic viscosity) of Current time step at cell center of each cell in a Zone.
	double  *ViscosityofPhase[5];
	///Brief: It holds the flag for Density to be constant or variable for phases. 0->constant Density; 1->variable Density
	int      flagRhoPhase[5];		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: It holds  Dynamic Viscocity multiplied by volume fraction of each phase
	double *ViscoeffDummy[5];
	/// Brief: It holds the flag for Dynamic Viscocity to be constant or variable for phases. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeffPhase[5];
	
	/// Brief: It is pointer type which stores the individual phase velocity in x direction of Current time step at cell center of each cell in a Zone.
	double  *U_phase[5];
	/// Brief: It is pointer type which stores the individual phase velocity in y direction of Current time step at cell center of each cell in a Zone.
	double  *V_phase[5];
	/// Brief: It is pointer type which stores the individual phase velocity in z direction of Current time step at cell center of each cell in a Zone.
	double  *W_phase[5];
	///Brief:It stores current time step value of  individual phase volume fraction at cell center of each cell in a zone.
	double  *VarthetaPhase[5];
	
	double *diffu;
	double *diffv;
	double *diffw;
// 	double *diffu[5];
// 	double *diffv[5];
// 	double *diffw[5];
	/// Current Intensity for Radiation Problem
	double *RadIntensity;
	///Current Refractive Index
	double *RefIndex;
	/// flag finding whether refractive index is constant or variable for zone
	int flagRefIndex;
	///Absorption coefficient
	double *kappa;
	///Scattering Coefficient
	double *sigma;
	///Volumetric Source Term
	double *RadEnergysource;
	/// Fluid Initial Temperature as Source
	double *RadGasTemp;
	///Radiation Solution Temperature
	double *RadSolTemp;
	
	double* VAR_VERTEX;
	
	// Added by adhiraj......These are the dispersed phase variables parameter.
	double *SlabVarthetaDisperse[6];
	double *SlabVxDisperse[6];
	double *SlabVyDisperse[6];
	double *SlabVzDisperse[6];
	double *SlabRhoDisperse[6];
	double *SlabNumDensity[6];//shashwat
	
	double *SlabDropFuelMassFrac[6];//shashwat for evaporation
	
	double *SlabTdrop[6];//shash 
	
	/// Brief: It is pointer type which stores the momentum of Disperse phase in x direction for Current time step at cell center of each cell in a Zone.
	double  *MomentumDispersex;
	double  *MomentumDispersexOld;
	/// Brief: It is pointer type which stores the momentum of Disperse phase in y direction for Current time step at cell center of each cell in a Zone.
	double  *MomentumDispersey;
	double  *MomentumDisperseyOld;
	//// Brief: It is pointer type which stores the momentum of Disperse phase in z direction for Current time step at cell center of each cell in a Zone.
	double  *MomentumDispersez;
	double  *MomentumDispersezOld;
	/// Brief: It is pointer type which stores the Density(Disperse phase) of Current time step at cell center of each cell in a Zone.
	double   *RhoDisperse;
	/// Brief: It is pointer type which stores the viscosity(Disperse phase) of Current time step at cell center of each cell in a Zone.
	double  *ViscosityDisperse;
 	double  *ViscoeffDummyDisperse;
	///Brief: It holds the flag for Density to be constant or variable for Disperse phase. 0->constant Density; 1->variable Density
	int      flagRhoDisperse;		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: It holds the flag for Dynamic Viscocity to be constant or variable for Disperse phase. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeffDisperse;
	
	///Brief:It stores current time step value of  disperse phase volume fraction at cell center of each cell in a zone.
	double  *VarthetaDisperse;

//Mrunalini
	
	///Brief:It stores current time step value of  dispersed phase velocity in x direction for boundary or interface slab.
	double  *U_Disperse;
	///Brief:It stores current time step value of  dispersed phase velocity in y direction for boundary or interface slab.
	double  *V_Disperse;
	///Brief:It stores current time step value of  dispersed phase velocity in z direction for boundary or interface slab.
	double  *W_Disperse;

	///Brief:It stores previous iteration value of  dispersed phase velocity in x direction for boundary or interface slab.
	double  *U_DisperseOld;
	///Brief:It stores previous iteration value of  dispersed phase velocity in y direction for boundary or interface slab.
	double  *V_DisperseOld;
	///Brief:It stores previous iteration value of  dispersed phase velocity in z direction for boundary or interface slab.
	double  *W_DisperseOld;
	
	///Brief:It stores current time step value of  dispersed phase mass velocity in x direction for boundary or interface slab.
	double  *U_Dmass;
	///Brief:It stores current time step value of  dispersed phase mass velocity in y direction for boundary or interface slab.
	double  *V_Dmass;
	///Brief:It stores current time step value of  dispersed phase mass velocity in z direction for boundary or interface slab.
	double  *W_Dmass;
	
	///Brief:It stores previous iteration value of  dispersed phase mass velocity in x direction for boundary or interface slab.
	double  *U_DmassOld;
	///Brief:It stores previous iteration value of  dispersed phase mass velocity in y direction for boundary or interface slab.
	double  *V_DmassOld;
	///Brief:It stores previous iteration value of  dispersed phase mass velocity in z direction for boundary or interface slab.
	double  *W_DmassOld;

	double *SlabUDmass[6];
	double *SlabVDmass[6];
	double *SlabWDmass[6];

	//diameter of the dispersed particle
	double *ParticleDia;
	//Mass flux for the disperse phase
	double *MassfluxDisperse;
	/// Brief: It is pointer type which stores the value of  dispersed phase number density for droplet phase  for Current time step at cell center of each cell in a Zone.
	double  *NumDensity;//shashwat
	//for evaporation model
	/// Brief: It holds the flag for diffusion coefficient for mass fraction of evaporated fuel to be constant or 	variable in Zone. 0->constant mass diffusivity 1 ->variable constant mass diffusivity
	int flagDropDiffcoeff; /* [S] Flag storing 0->constant diffcoeff; 1->variable diffcoeff */
	
	double *DropFuelMassFrac;//shashwat for evaporation
	
	/// Brief: It stores the mass diffusivity of the boundary or interface slab. 
	double   *DropDiffcoeff;		/* [S] Array for mass diffusivity for evaporated fuel mass fraction    */
	
	double *EvapoConst;
	double *Evaporate;
	double *FuelSurfMassFrac;
	//For Heat-Up period of evoparation
	double *Tdrop;
	
	double *CPdrop;
	
	int flagCPdrop;
	
	double *QDroptoGas;
	
	double *QLatent;
	
	//end of disperse phase
	
	//Added by chandrasekhar --- For Turbulent MPFV
	double *SlabVx_Old[6];
	double *SlabVy_Old[6];
	double *SlabVz_Old[6];
	double *SlabDummyeffViscoeffofphase[5][6];
	
	double *SlabVcx_Old[6];
	double *SlabVcy_Old[6];
	double *SlabVcz_Old[6];
	//End --- chandrasekhar
	//Added to make slabs of Diffusion coefficient for particle phase in case of Turbulent gas-particle flow
	double *SlabDummyViscosity_disperse[6];
	//Added to make slabs of Diffusion coefficient for volume fraction equation in case of Turbulent Gas-Particle flow
	double *SlabDiffVolfrac_disperse[6];
	
	double *Diff_Disperse;
	
} FlowSolution;
/** \struct FlowSolution
	This is FlowSolution Structure. It contains Solution information about the Zone in the curent time
	step. All the flow variables are defined and stored in array. All the information in this structure is
	of Current time step.
*/
/* -------------------------------------------------------------------------------------------
      NOTE:
   -------------------------------------------------------------------------------------------
	1) 'GridLocation' specifies the location of the solution data with respect to grid.
	   For example, If the flow velocities are defined at the cell center then the
	   GridLocation for 'flow velocities' would be 'CellCenter'. Similarly if the
	   pressure is defined at the grid vertices then the GridLocation 'for pressure'
	   would be "Vertex".
	   For all possible values of GridLocation, see 'cgnslib.h' (line 85-93).
	   For detailed description of GridLocation read the Developer's Manual.
	2) The size of the arrays defined above depends upon
		- The GridLocation of the solution data.
		- Cell Dimension
  --------------------------------------------------------------------------------------------
*/

typedef struct
{
	/// Brief: It store the Previous flow solution name.
	char   PrevFlowSolutionName[33];/* [C] Name of the Flow Solution Node   */
	/// Brief: `GridLocation' specifies the location of the solution data with respect to grid. For all possible values of GridLocation, see 'cgnslib.h' (line 85-93). For detailed description of GridLocation read the Developer's Manual.
	int    GridLocation;		/* [C] Location of the Grid (Read NOTE) */
	/// Brief: It is a pointer which stores the density of the previous time step at the cell center of each cell in a Zone.
	double   *Rho,*Rho_prev;			/* [S] Array for Density                */
	/// Brief: It is a pointer which stores the Dynamic Viscocity of previous time step at cell center of each cell in a Zone.
	double   *Viscoeff;		/* [S] Array for Viscosity              */
	/// Brief: It is a pointer which stores the Thermal Conductivity of previous time step at cell center of each cell in a Zone.
	double   *Diffcoeff;		/* [S] Array for Thermal Diffusivity    */
	/// Brief: It is a pointer which stores the Specific Heat at constant pressure of previous time step at cell center of each cell in a Zone.
	double	 *SpecificHeat;
	/// Brief: It is a pointer which stores the U velocity of previous time step at cell center of each cell in a Zone.
	double   *u;			/* [S] Array for X-Velocity             */
	/// Brief: It is a pointer which stores the V velocity of previous time step at cell center of each cell in a Zone.
	double   *v;			/* [S] Array for Y-Velocity             */
	/// Brief: It is a pointer which stores the W Velocity of previous time step at cell center of each cell in a Zone.
	double   *w;			/* [S] Array for Z-Velocity             */
	/// Brief: It is a pointer which stores the Mass Flux of previous time step at each face of each cell in a Zone.
	double   *Massflux;
	/// Brief: It is a pointer which stores the U mass velocity of previous time step at cell center of each cell in a Zone.
	double   *umass;
	/// Brief: It is a pointer which stores the V mass velocity of previous time step at cell center of each cell in a Zone.
	double   *vmass;
	/// Brief: It is  a pointer which stores the W mass velocity of previous time step at cell center of each cell in a Zone.
	double   *wmass;
	/// Brief: It is a pointer which stores the Pressure of previous time step at cell center of each cell in a Zone.
	double   *P;			/* [S] Array for Pressure               */
	/// Brief: It is a pointer which stores the Temperature of previous time step at cell center of each cell in a Zone.
	double   *T;			/* [S] Array for Temperature            */
	/// Brief: It is a pointer which stores the enthalpy of previous time step at cell center of each cell in a Zone.
	double   *enthalpy;	  	//Created by bhans
	/// Brief: It is a pointer which stores the LatentheatContent of previous time step at cell center of each cell in a Zone.
	double   *Latent_Heat_Content;  //Created by bhans
	/// Brief: It is an array of pointers which stores all the 10 Scalars of previous time step at cell center of each cell in a Zone.
	double   *Scal[10];		/* [S] Array for scalars                */
	/// Brief: It is an array of pointers which stores all the 10 Scalar Diffusion Coefficents of previous time step at cell center of each cell in a Zone.
	double   *Scalcoeff[10];	/* [S] Array for diff coeff for scalars */

        //Combustion Variables
	///Brief: It stores previous time step values of mixture fraction in zones				
        double	 *MixFrac;  //shyam
	///Brief: It stores previous time step values of fuel mass fraction in zones					
        double   *FuelMassFrac; //shyam
	///Brief: It stores previous time step values of oxygen mass fraction in zones					
	double	 *OxiMassFrac; //Shyam
	///Brief: It stores previous time step values of inert gas mass fraction in zones					
        double   *InMassFrac[10];  //Shyam
	///Brief: It stores previous time step values of product mass fraction in zones					
        double   *ProdMassFrac; //Shyam

	
	// The following are the ELECTRICAL parameters	
	/// Brief: It is pointer type which stores the electric potential of previous time step at cell center of each cell in a Zone. 
	double	 *Potential;		/* [S] Array for Electric Potential     */
	/// Brief: It is pointer type which stores the X-direction electric field of previous time step at cell center of each cell in a Zone. 
	double   *EfieldX;		/* [S] Array for X-Dirn Electric Field  */
	/// Brief: It is pointer type which stores the Y-direction electric field of previous time step at cell center of each cell in a Zone. 
	double   *EfieldY;		/* [S] Array for Y-Dirn Electric Field  */
	/// Brief: It is pointer type which stores the Z-direction electric field of previous time step at cell center of each cell in a Zone. 
	double   *EfieldZ;		/* [S] Array for Z-Dirn Electric Field  */
	/// Brief: It is pointer type which stores the electric flux of previous time step at cell center of each cell in a Zone. 
	double   *Electricflux;		/* [S] Array for Electric Field Dlux    */
	/// Brief: It is pointer type which stores the Ion Charge Density of previous time step at cell center of each cell in a Zone. 
	double   *RhoIon;		/* [S] Array for Space Charge Density   */
	/// Brief: It is pointer type which stores the dielectric constant of previous time step at cell center of each cell in a Zone. 
	double   *DiElecConst;		/* [S] Array for Di-Electric Constant   */
	/// Brief: It is pointer type which stores the Ion Mobility rate of previous time step at cell center of each cell in a Zone. 
	double   *MobIon;		/* [S] Array for Ion Mobility Rate      */
	/// Brief: It is pointer type which stores the charge diffusion coefficient of previous time step at cell center of each cell in a Zone. 
	double 	 *ChargeDiffcoeff;	/* [S] Stores the value of Charge Diffusion Coefficient          */
	
	// The following are the LMMHD parameters (Implemented by naveen)	
	/// Brief: It is pointer type which stores the electric potential of previous time step at cell center of each cell in a Zone. 
	double	 *PotentialMHD;		/* [S] Array for Electric Potential     */
	/// Brief: It is pointer type which stores the X-direction Magnetic field of previous time step at cell center of each cell in a Zone. 
	double   *MfieldX;		/* [S] Array for X-Dirn Magnetic Field  */
	/// Brief: It is pointer type which stores the Y-direction Magnetic of previous time step at cell center of each cell in a Zone. 
	double   *MfieldY;		/* [S] Array for Y-Dirn Magnetic Field  */
	/// Brief: It is pointer type which stores the Z-direction Magnetic field of previous time step at cell center of each cell in a Zone. 
	double   *MfieldZ;		/* [S] Array for Z-Dirn Magnetic Field  */
	/// Brief: It is pointer type which stores the term in the RHS of the PotentialMHD poisson equation for Current time step at each face of each cell in a Zone.
	double   *PotentialMHDpoissonRHS; 	
	/// Brief: It stores the Electrical conductivity of the liquid metal of previous time step at cell center of each cell in a Zone.  
	double   *Eleconductivity;	/* [S] Array for Electrical conductivity*/
	double   *TotalMfieldX;//rajendra
	double   *TotalMfieldY;
	double   *TotalMfieldZ;
	double   *Permeability;
	double   *Magneticflux;//rajendra
	double   *Psi;
	double   *Magdiff;//rajendra

	
	//  Groundwater parameters 
	/// Brief{GROUNDWATER}It stores the value of pressure head used in groundwater flow
	double  *H;
	/// Brief{GROUNDWATER}It stores the value of water content used in groundwater flow
	double  *WaterContent;

	double *Flow_Potential;
	
	 //Added by sushil.....These are the multiphase variables parameter.
         /// Brief: It is pointer type which stores the overall momentum in x direction for previous time step at cell center of each cell in a Zone.
	double  *Momentumx;
	/// Brief: It is pointer type which stores the overall momentum star values in x direction for previous time step at cell center of each cell in a Zone.
	double *Momentumxstar;
		
	/// Brief: It is pointer type which stores the overall momentum in y direction for previous time step at cell center of each cell in a Zone.
	double  *Momentumy;
	/// Brief: It is pointer type which stores the overall momentum star values in y direction for previous time step at cell center of each cell in a Zone.
	double *Momentumystar;
		
	///Brief:It stores previous time step value of overall  momentum in z direction for boundary or interface slab.
	double  *Momentumz;
	/// Brief: It is pointer type which stores the overall momentum star values in z direction for previous time step at cell center of each cell in a Zone.
	double *Momentumzstar;
	/// Brief: It is pointer type which stores the individual momentum of each phase in x direction for previous time step at cell center of each cell in a Zone.
	double  *Momentumofphasex[5];
	
	/// Brief: It is pointer type which stores the individual momentum of each phase in y direction for previous time step at cell center of each cell in a Zone.
	double  *Momentumofphasey[5];
	
	//// Brief: It is pointer type which stores the individual momentum of each phase in z direction for previous time step at cell center of each cell in a Zone.
	double  *Momentumofphasez[5];
	
	/// Brief: It is pointer type which stores the Density(phasic density) of previous time step at cell center of each cell in a Zone.
	double   *RhoofPhase[5];
	/// Brief: It is pointer type which stores the viscosity(phasic viscosity) of previous time step at cell center of each cell in a Zone.
	double  *ViscosityofPhase[5];
	
	double *ViscoeffDummy[5];
	
	
	/// Brief: It is pointer type which stores the individual phase velocity in x direction of previous time step at cell center of each cell in a Zone.
	double  *U_phase[5];
	/// Brief: It is pointer type which stores the individual phase velocity in y direction of previous time step at cell center of each cell in a Zone.
	double  *V_phase[5];
	/// Brief: It is pointer type which stores the individual phase velocity in z direction of previous time step at cell center of each cell in a Zone.
	double  *W_phase[5];
	
	///Brief:It stores previous  time step value of  individual phase volume fraction at cell center of each cell in a zone.
	double  *VarthetaPhase[5];
	
	/// Previous value of Intensity for Radiation Problem
	double *RadIntensity;
	///Previous Refractive Index
	double *RefIndex;
	
	
	// Added by adhiraj......These are the dispersed phase variables parameter.
	double *SlabVarthetaDisperse[6];
	double *SlabVxDisperse[6];
	double *SlabVyDisperse[6];
	double *SlabVzDisperse[6];
	double *SlabRhoDisperse[6];
	double *SlabNumDensity[6];//shashwat
	
	double *SlabDropFuelMassFrac[6];//shashwat for evaporation
	/// Brief: It is pointer type which stores the momentum of Disperse phase in x direction for Current time step at cell center of each cell in a Zone.
	
	double *SlabTdrop[6];//shash for evaporation
	/// Brief: It is pointer type which stores the momentum of Disperse phase in x direction for Current time step at cell center of each cell in a Zone.
	double  *MomentumDispersex;
	double  *MomentumDispersexOld;
	/// Brief: It is pointer type which stores the momentum of Disperse phase in y direction for Current time step at cell center of each cell in a Zone.
	double  *MomentumDispersey;
	double  *MomentumDisperseyOld;
	//// Brief: It is pointer type which stores the momentum of Disperse phase in z direction for Current time step at cell center of each cell in a Zone.
	double  *MomentumDispersez;
	double  *MomentumDispersezOld;
	/// Brief: It is pointer type which stores the Density(Disperse phase) of Current time step at cell center of each cell in a Zone.
	double   *RhoDisperse;
	/// Brief: It is pointer type which stores the viscosity(Disperse phase) of Current time step at cell center of each cell in a Zone.
	double  *ViscosityDisperse;
	double  *ViscoeffDummyDisperse;
	///Brief: It holds the flag for Density to be constant or variable for Disperse phase. 0->constant Density; 1->variable Density
	int      flagRhoDisperse;		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: It holds the flag for Dynamic Viscocity to be constant or variable for Disperse phase. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeffDisperse;
	
	///Brief:It stores current time step value of  disperse phase volume fraction at cell center of each cell in a zone.
	double  *VarthetaDisperse;

//Mrunalini	
	///Brief:It stores current time step value of  dispersed phase velocity in x direction for boundary or interface slab.
	double  *U_Disperse;
	///Brief:It stores current time step value of  dispersed phase velocity in y direction for boundary or interface slab.
	double  *V_Disperse;
	///Brief:It stores current time step value of  dispersed phase velocity in z direction for boundary or interface slab.
	double  *W_Disperse;

	///Brief:It stores previous iteration value of  dispersed phase velocity in x direction for boundary or interface slab.
	double  *U_DisperseOld;
	///Brief:It stores previous iteration value of  dispersed phase velocity in y direction for boundary or interface slab.
	double  *V_DisperseOld;
	///Brief:It stores previous iteration value of  dispersed phase velocity in z direction for boundary or interface slab.
	double  *W_DisperseOld;
	///Brief:It stores previous time step value of  dispersed phase mass velocity in x direction for boundary or interface slab.
	double  *U_Dmass;
	///Brief:It stores previous time step value of  dispersed phase mass velocity in y direction for boundary or interface slab.
	double  *V_Dmass;
	///Brief:It stores previous time step value of  dispersed phase mass velocity in z direction for boundary or interface slab.
	double  *W_Dmass;

	double *SlabUDmass[6];
	double *SlabVDmass[6];
	double *SlabWDmass[6];
	
	double *RhoG_face;	//It stores the face values of previous time step density of gas-phase
	double *VarthetaD_face; //It stores the face values of previous time step normalized volume fraction of dispersed phase

	//diameter of the dispersed particle
	double *ParticleDia;
	//Mass flux for the disperse phase
	double *MassfluxDisperse;
	/// Brief: It is pointer type which stores the value of  dispersed phase number density for droplet phase  for Current time step at cell center of each cell in a Zone.
	double  *NumDensity;//shashwat
	double *DropFuelMassFrac;//shashwat for evaporation
	
	/// Brief: It stores the mass diffusivity of the boundary or interface slab. 
	double   *DropDiffcoeff;		/* [S] Array for mass diffusivity for evaporated fuel mass fraction    */
	//for evaporation model
	double *EvapoConst;
	double *Evaporate;
	double *FuelSurfMassFrac;
	//For Heat-Up period of evoparation
	double *Tdrop;
	
	double *CPdrop;
	
	double *QDroptoGas;
	
	double *QLatent;
	
	//end of disperse phase
} PrevFlowSolution;
/** \struct PrevFlowSolution
	This is PrevFlowSolution Structure. It contains Solution information about the Zone in Previous time
	step. All the flow variable are define and stored in array. All the information in this structure is
	one time step lagging.
*/

/* **-------------------------------Data Structure for Zone Geometry--------------------------------** */
typedef struct
{
	/// Brief: It is pointer type, stores the volume of each cell of Zone.
	double *Volume;			/*[A] Array for storing the Cell Volumes(3-D)*/
	/// Brief: It is a pointer type which, stores the surface area of all three components of each face and of each cell. Therefore each cell of a Zone requires 18 element array.
	double *SurfaceArea;		/*[A]Array for storing the Surface Areas of the cells(7-D)*/
	/// Brief: It is a pointer type which, stores the diffusion length of all three components of each face and of each cell. Therefore each cell of Zone requires 18 element array.
	double *Beta;			/*[A]Array for storing the Diffusion Coefficients(7-D)*/
	/// Brief: It is a pointer which, stores the cell center distance in all three direction.
	double *deltaX;
	/// Brief: It is a pointer which, stores the Inverse Matrix for calculating gradients.(See report)
	double *DeltaInverse;
	/// Brief: This pointer stores the min. distance from a cell center to any arbitrary solid wall ... used in SST turbulence model
	double *minWallDistance;
	
	double *DmF;
	///Calculation of Control angle
	double *Omega;
	double *VectorX, *VectorY, *VectorZ;
	double *NormalX, *NormalY, *NormalZ;
	
}ZoneGeometry;
/** \struct ZoneGeometry
	This is Zone Geometry structure. It contains the information of Zone geometry. It contain pointers of
	cell volumes, surface areas of cell faces and diffusion length of each cell in the zone.
*/

/* ** ----------------------------- Data Structure for Grid Coordinate Node --------------------------- ** */
typedef struct
{
	/// Brief: It Stores the name of Grid Coordinate Node. The name is read from the CGNS Grid file.
	char   GridcoordinateName[33];	/* [C] Name of the Grid Coordinate Node */
	/// Brief: It Stores the Number of Grid Coordinate Node. The number would be 2 for 2-D or 3 for 3-D cases.
	int    Number_of_Coordinates;	/* [C] Number of Data Arrays containing Grid Coordinate */
	/// Brief: It store information about the presision of the variable, Single or Double presision.
	int    Datatype;		/*[A]*/
	/// Brief: It store the Name of the first array of grid coordinate, generally it is X coordinate.
	char   Array_Name1[33];		/* [C] Name of the first Array */
	/// Brief: It is pointer type in which X coordinate is stored in a single long array of a Zone.
	double   *X;			/* [C] First Array  (for Coordinate i) */
	/// Brief: It store the Name of the second array of grid coordinate, generally it is Y coordinate.
	char   Array_Name2[33];		/* [C] Name of the Second Array */
	/// Brief: It is pointer type in which Y coordinate is stored in a single long array of a Zone.
	double   *Y;			/* [C] Second Array (for Coordinate j) */
	/// Brief: It store the Name of the third array of grid coordinate, generally it is Z coordinate.
	char   Array_Name3[33];		/* [C] Name of the Third Array */
	/// Brief: It is pointer type in which Z coordinate is stored in a single long array of a Zone.
	double   *Z;			/* [C] Third Array  (for Coordinate k) */
	/// Brief: It is an array which stores the minimum X,Y,Z coordinates, of all the coordinates making a zone.
	double Min[3];
	/// Brief: It is an array which stores the maximum X,Y,Z coordinates, of all the coordinates making a zone.
	double Max[3];
}GridCoordinate;
/** \struct GridCoordinate
	This is GridCoordinate Structure. It contains Grid information about the Zone.
	All the coordinate are define and stored here.
*/

/* ** ------------------------- Data Structure for the Zone Boundary ---------------------------------- ** */
typedef struct
{
	/// Brief: It is a structure which stores information for inplementation of various boundary conduction.
	Dataset *Data_Set; /*Contains Data that defines boundary conditions of a boundary patch*/
	/// Brief: It stores the RMS value at a boundary while implementing the Numann boundary condition. When we iterate the equation in boundary node it gives rms residual which is averaged and stored here. 
	double  NueRms;
	/// Brief: It stores the name of boundary patches. 
	char	BCName[33];		/* [C] Name of the boundary Patch Node */
	/// Brief: It is pointer type which stores the index range of boundary. It is a 6 element array. 
	int	*Point_Range;		/* [C] Point Index Range of the boundary patch */
	/// Brief: It is pointer type which stores the cell range of boundary. It is a 3 element array. 
	int	*Cell_Range;		/* [A] Cell Index Range of the boundary patch */
	/// Brief: It stores the boundary tag which give information about the type of boundary. for more information refer CGNS Manual.
	int	Boundary_Tag;		/* [C] Boundary Tag specifying the type of boundary */
	/// Brief: It stores the direction of boundary patch. for reference to direction see file `Mish.h'.
	int	Direction;
	/// Brief: It stores the distance of first grid point from the wall. This is calculated only in Turbulent flow case and only on walls.
	double *deltaP;
	/// Brief: It stores the friction velocity at the wall. This is calculated only in Turbulent flow case and only at the walls.
	double *uTau;
	/// Brief: It stores all the three components of pressure force over a boundary calculated using PressureForce() function. It is three element array. 
	double PressForce[3];		/*[S]*/
	/// Brief: It stores the sum of normal gradient of temperature at boundary. 
	double SumNorGradT;		/*[S]*/
	/// Brief: It stores the sum of normal gradient of all 10 sclar at boundary. It is 10 element array. 
	double SumNorGradScal[10];	/*[S]*/
	/// Brief: It is Zone Slab structure. It stores the flow variable, and boundary slab information. 
	ZoneSlab *ZSlab;		/*[A] Information of Zone Slab*/
	/// Brief: It stores the Average Nusselt number in the boundary for temperature calculated from function CalcNusseltNo() function. 
	double AvNusseltNoT;
	/// Brief: It stores the Average Nusselt number in the boundary for all 10 Sclar calculated from function CalcNusseltNo() function. 
	double AvNusseltNoScal[10];
	
	///It is used to store face center or projected face normal point coordinates of a boundary/conjugate interface 
	double *FaceCenter;
}ZoneBoundary;
/** \struct ZoneBoundary
	This is ZoneBoundary Structure. It contains all Boundary information about the Zone. 
	It is pointer type which have all-boundary information like name, index range, direction etc. 
*/

/* ** ------------------------------- Data Structure for Interface ------------------------------------ ** */
typedef struct Interface
{
	/// Brief: It store the interface name.
	char 	Interface_Name[33];	/* [C] Name of the interface 					*/
	/// Brief: It store the Adjacent Zone name. 
	char	Adjoining_Zone_Name[33];/* [C] Name of the adjoining zone 				*/
	/// Brief: It is a pointer which stores Index range of the interface in the current Zone index convention. it is 6 element array.
	int	*Point_Range;	   	/* [C] Point Index Range of the interface in the current zone 	*/
	/// Brief: It is a pointer which, stores Index range of the interface in the Adjacent Zone index convention. it is 6 element array.
	int 	*Point_Donor_Range;	/* [C] Point Index Range of the interface in the adjoining zone */
	/// Brief: It is a pointer which, stores transformation matrix in short to go from current Zone to adjacent Zone. It is 3 element array. 
	int	*Transform;		/* [C] Short Hand "current->adjoining" Transform rule 		*/
	/// Brief: It is a pointer which, stores transformation matrix in expanded form to go from current Zone to adjacent Zone. It is 9 element array. 
	int	*Transformation_Matrix;	/* [A] Expanded matrix for current->adjoining transformation 	*/
	/// Brief: It is a pointer which, stores constant to be added while transforming the indices. It is 3 element array. 
	int	addConstant[3];
	/// Brief: It store the direction of interface. For direction refer to``Mish.h" file.
	int 	Direction;		/* [A] Direction of the Interface (See NOTE)			*/
	/// Brief: It is structure of Zone Slab which contain Solution information of Intarface Slab.
	ZoneSlab *ZSlab;		/* [A] Information of Zone Slab					*/
	/// Brief: It stores the RMS value of interface while implementing interface condition in Conjugate Heat Transfer Problems. 
	double ConjugateRms;
	/// Brief: It stores the Average Nusselt Number at interface for Temperature. 
	double AvNusseltNoT;
	/// Brief: It is a pointer which, stores the Average Nusselt Number at interface for 10 Scalars.
	double AvNusseltNoScal[10];
	/// Brief: It stores the distance of first grid point from the wall. This is calculated only in Turbulent flow case and only on walls (Conjugate Heat Transfer).
	double *deltaP;
	/// Brief: It stores the friction velocity at the wall. This is calculated only in Turbulent flow case and only on walls (Conjugate Heat Transfer).
	double *uTau;
/**	This tag stores '1' or '2', depending on the zone type, for all the conjugate interfaces and '0' for the rest and is assigned with these values in the function initialize_T().
	'0' indicates -- Current and adjacent zones both are of the same type: either SOLID or FLUID.
	'1' indicates -- Current zone is FLUID and adjacent zone is SOLID; 
	'2' indicates -- Current zone is SOLID and adjacent zone is FLUID.
 */
	int ConjugateIfaceTag;
	///This pointer points to the adjoining zone at a given interface.
	struct Zone *AdjoiningZonePtr;
}Interface;
/** \struct Interface
	This is Interface Structure. It contains all information about the Zone Interface. 
*/
/* -------------------------------------------------------------------------------------------
      NOTE:
   -------------------------------------------------------------------------------------------
	1) Read "Misc.h" for description of Direction
   -------------------------------------------------------------------------------------------*/
   
  /* ** --------------------------- Data Structure for Multigrid ------------------------------------------- ** */
typedef struct
{
	/// Breif: Pointer for grid size.
	int *GSize;
	/// Breif: These are pointers to store 19 coefficients arising from Diffusion formulation.
	double *DiffusionCoeff_P,*DiffusionCoeff_nb[18];
	double *DiffusionCoeffPhi_P,*DiffusionCoeffPhi_nb[18];
	///Diffusion coefficients in case of problems having variable diffusion coefficient.
	double *Diffusion_Diag[6],*Diffusion_nb[36],*Diffusion_Nonortho[24];
	
	
	///These are the pointers to convective coefficients depending upon the UpwindScheme.  Storing them will enable us to remove unnecessary computations.

	double *Conv_Coeff_P,*Conv_Coeff_nb[6];
	
	/// Breif: These are pointers to store 19 coefficients of variable for CD and NS parts.
	double *ConvDiffusionCoeffPhi_P,*ConvDiffusionCoeffPhi_nb[18];
	
	
// 	double *DiffusionCoeffPhi_P_LMMHD,*DiffusionCoeffPhi_nb_LMMHD[18];//rajendra
// 	double *DiffusionCoeffT_P,*DiffusionCoeffT_nb[18];
// 
// 	double *DiffusionCoeffEnthalpy_P,*DiffusionCoeffEnthalpy_nb[18];//for enthalpy
// 
// 	double *DiffusionCoeffS_P[10],*DiffusionCoeffS_nb[10][18];
// 
 	double *DiffusionCoeff_Poten,*DiffusionCoeff_Poten_nb[18];
// 	double *DiffusionCoeffRhoIon_P,*DiffusionCoeffRhoIon_nb[18];
// 
// 	double *DiffusionCoeffMixF_P,*DiffusionCoeffMixF_nb[18];
	
	
// 	double *ConvDiffusionCoeffPhi_P_LMMHD,*ConvDiffusionCoeffPhi_nb_LMMHD[18];//rajendra
// 	double *ConvDiffusionCoeffT_P,*ConvDiffusionCoeffT_nb[18];//For Temperature
// 
// 	double *ConvDiffusionCoeffEnthalpy_P,*ConvDiffusionCoeffEnthalpy_nb[18];//For Enthalpy equation
// 
// 	double *ConvDiffusionCoeffS_P[10],*ConvDiffusionCoeffS_nb[10][18];// For Scalars
// 
// 	double *ConvDiffusionCoeffRhoIon_P,*ConvDiffusionCoeffRhoIon_nb[18];
// 
// 	double *ConvDiffusionCoeffMixF_P,*ConvDiffusionCoeffMixF_nb[18];
	/// Brief: These are the pointers to store the error and the right hand side at each grid point for all the grid levels.
	double *E,*GRHS,*Gresidual;
	/// Breif: Pointer to the grid variable and slab variable.
	double *V,*SlabV[6],*DummyV,*DummySlabV[6];
}Grid;

  /* ** ----- Data Structure for Parameters to be scanned from FILE (bin/Settings/InputParameters.dat)------------- ** */
typedef struct
{
	int iwrite;
	int CountFileWrite;
	double VelocityRelaxation;
	double TurbRelaxation;
	double OmegaRelaxation;
	double PressureRelaxation;
	double PotentialRelaxation;
	double ChargeRelaxation;
	double DummyViscocity;
	double PressureConvergenceCriteria;
	double ConvergenceCriteria;
	double EfieldTimeStep;
	double PotentialConvergence_TDMA;
	double ChargeConvergence_GS;
	double WireChargeUpdateData;
	double OverallmomentumRelaxation;
	double MomentumofPhaseRelaxation;
}Input;

/* ** --------------------------------- Data Structure for Zone --------------------------------------- ** */
typedef struct Zone
{
	/// Brief: It stores the name of the Zone given by the CGNS file genetaror software. 
	char   zonename[33];		/* [C] Name of the 'Zone' Node */
	/// Brief: User can give his own choise of name to a particular name. For conjugate heat transfer case the user is advised to enter the name of the Zone as Solid/Fluid.
	char	SolidFluidName[33];	/* [C] Name of the 'Zone' Node in term of Solid Liquid*/
	/// Brief: It holds the flag for Zone type. Zone Type - 1->Unstructured, 2->Structured 
        int    Zonetype;		/* [C] Zone Type - Structured or Unstructured */
	/// Brief: It contain the size of the zone in term of maximum vertex indices. i.e NI, NJ, NK. It is a pointer, from which you can acces all three indices in 3-D case.
	int    *Size;   		/* [C] Vertex Size, Cell Size & Boundary Cell Size */
	/// Brief: It hold flag for Solid or Fluid. Fluid->0 ; Solid->1
	int	flagSolidFluid;		/*fluid->0;solid->1*/
	/// Brief: It hold flag for Heat Generation for particular Zone. HeatGeneration->0 ; NO HeatGeneration->1
	int	flagHeatGeneration;	/*Yes->0;No->1*/
	/// Brief: Tag to denote whether the heat generation value is constant or not
	int	flagHeat_con_var;
	/// Brief: It stores the number of interfaces in a particular zone. 
	int	Number_of_Interfaces;	/* [C] Number of Interfaces in the zone */
	/// Brief: It is a pointer to Interfaces structure of a Zone, which contain all information of an interface in a Zone.  
	Interface *Iface;		/* [C/A]Interface info of all the interfaces in the zone */
	/// Brief: It stores number of grid coordinate nodes, 2 for 2-D and 3 for 3-D. 
        int    Number_of_Gridcoor_nodes;/* [C] Number of grid coordinate nodes */ 
	/// Brief: It is pointer to CgnsGrid of a Zone. It contain all the information about grid in a Zone.
	GridCoordinate *CgnsGrid;	/* [C] Grid Coordinate Info for the Zone */
	/// Brief: It is pointer to Flow Solution structure of a Zone. It contain solution information of a particular Zone.
        FlowSolution *Flowsol;		/* [S] Flow Solution Info for the Zone */
	/// Brief: It is pointer to Previous Flow Solution of a Zone. It contain solution information of previous time step for a Zone.
	PrevFlowSolution *PrevFlowsol;
	/// Brief: It stores the number of Boundaries in a particular Zone. 
        int    Number_of_BC_Patches;	/* [C] Number of boundary patches in the zone */
	/// Brief: It is a pointer to Zone Boundary of a Zone. It stores information about the Zone Boundary. 
	ZoneBoundary *ZoneBc;		/* [C/A] Zone Boundary Info for the Zone */
	/// Brief: It is a pointer to Zone Geometry of a Zone. It stores information about the Zone Geometry. 
	ZoneGeometry *ZGeom;
	
	// The following are the ELECTRICAL parameters	
	/// Brief: It is used for updating ion charge density on wire in each zone exactly once after specifie number of timesteps. 
	int	ZoneWireTag;	/*For updating charge on wire only once in the timestep*/
// 	Groundwater flow start
	///Brief: This one is stored as tag whether to calculate dispersion coefficent for zone from velocity  or take as given by input directly dispersion coeff.. 0 for not calculating and 1 for calcualtiong from velocity
       int CalDispersion;
       
	int Difftag;
       ///Brief: longitudinal dispersivity of porous medim 
       double alphaL;
       ///Brief: Transere dispersivity of porous medium
       double alphaT;
       /// Brief: porosity of the mdium
       double  Porosity;
       /// Brief: used to make choice over two model Haverkamp and  VGM model
        int flagtype_model;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double Alpha;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double Theta_residual;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double Theta_saturation;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double n;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double A;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double B;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double C;
       /// Brief: Paramters to be used in VGM model and Haverkamp in case of unsaturated groundwater fow
       double D;
       /// Brief: This dummy arry of 6 double pointers store the temperature values of all the slabs (both Interface and Boundary) in case of Conjugate Heat Transfer with Turbulence.
       double *dummySlabTemp[6];
       ///Brief: This tag stores '1' if a given zone contains at least one conjugate interface; otherwise stores '0'.
       int ConjugateZoneTag;
// end		
	Grid *grid;
	/// Brief: zonalCoarseGridIndex should be same as Base->CoarseGridIndex
	int zonalCoarseGridIndex; //CGI
	
	double *DeferredCorrCoeffPhi_P,*DeferredCorrCoeffPhi_nb[12];
	///All these coefficients for different variables are hashed and only one set of coefficients will be used	
	
// 	double *DeferredCorrCoeffPhi_P_LMMHD,*DeferredCorrCoeffPhi_nb_LMMHD[12];//rajendra
// 	double *DeferredCorrCoeffT_P,*DeferredCorrCoeffT_nb[12];
// 
// 	double *DeferredCorrCoeffEnthalpy_P,*DeferredCorrCoeffEnthalpy_nb[12];
// 
// 	double *DeferredCorrCoeffS_P[10],*DeferredCorrCoeffS_nb[10][12];
// 	double *DeferredCorrCoeffMixF_P,*DeferredCorrCoeffMixF_nb[12];
// 	
// 	double *DeferredCorrCoeffRhoIon_P,*DeferredCorrCoeffRhoIon_nb[12];
	
	double *DeferredCorrCoeffDummy_P;//Just for storing zero values for the solver which doesn't require convection terms.
	
	///Whethear Source is there in this zone
	int flagRadHeatGeneration;
	///If Yes Constant or variable
	int flagRadHeatConVar;
	/// Zone fluid Temperature
	int flagZoneTemp; 
	
	int BC_IFACE_TAG[6];///This Tag indicates whether a particular face is at interface or boundary tag=1 means boundary and tag=2 means it's an interface Bhans
	
	void *ZONE_FACES[6];///It stores the void pointer of the face of the zone( for both BC and IFACE) Bhans
				
	int ZONE_NEIGHBOUR_NUM[6];///Stores the value of the neighbouring zone number(index starting from zero) Bhans.  Stores -1 for boundaries
	
	//added by adhiraj for dispersed phase flow
	double *DragX;
	double *DragY;
	double *DragZ;
	double *t_p;
	//added by chandrasekhar for liquid-gas bubbly flow
	double *LiftX;
	double *LiftY;
	double *LiftZ;
	//Added by sushil
	int FlagMPFV;
	int Nphase;
	double *Dragcoeff;
	double *Bx;
	double *By;
	double *Bz;
	///Improved G-S Terms
	double *Residue,*Correction,*AP,*Correction_Prev,*AP_Prev;//sandeep
	///Bi-Conjgate Gradient terms... sandeep.....*Residue can be used by both the Algo.s.
	double *Residue_init,*Residue_Old,*Residue_New,*Pi,*Vi,*Si,*Ti;//these are Intermediate vectors used
	///Brief: It stores 0 for Mixture Turbulence model and 1 for dispersed turbulence model in case of Turbulent MPFV
	int flagTurbMPFVmodel;
	
}Zone;
/** \struct Zone
	This is Zone structure, it is define in Base. It holds the pointer to all zone in the computational 
	domain. It contains the information of interface, Grid coordinate, Flow solution, boundary, and geometry 
	in terms of pointer. It also contains the name of Zone from CGNS file and in terms of Solid/Fluid.  
*/

// ------------------------------ Data Structure for Turbulence Base Parameters --------------------------------
typedef struct
{
	///Flag for auto-initialization of turbulence variables viz., k and epsilon
	int flagAutoInit;
	///Flag for Durbin's realizability constraint on time-scale
	int flagRealizability;
	///Reference velocity for limiting turbulent kinetic energy and limiting values for turbulence variables
	double Uref,kmin,kmax,epsmin;
	///Variables used in context of limiting turbulence terms
	int TProductionMaxLimitCount,TimeScaleMinLimitCount,TimeScaleMaxLimitCount;
	///``numTurbScalars'' stores number of turbulence scalars to be solved
	int numTurbScalars;
	///``numDualTimeIter'' stores number of dual time stepping loop iterations to be performed in case of two-equation turbulence models. It is equal to the `1' for one-equation turbulence models
	int numDualTimeIter;
	///``DualTimeStep'' stores the dual time step for all the turbulence models. It is equal to the actual timestep for one-equation turbulence models
	double DualTimeStep;
	///flags for various limiters used in turbulence models
	int flagLimTimeScale,flagLimKEps,flagLimTVisc,flagLimProduction;
}ParametersTurb;
/** \struct ParametersTurb
	This is ParametersTurb structure. It is defined in Parameters. It contain information 
	about the turbulence base parameters and various flags for solving the problem in the computational domain.
 */

typedef struct
{
	/// Brief: This zone pointer points to the zone in which the particle resides currently.
	Zone *zone;
	/// Brief: This is the i index of the cell in which the particle resides currently.
	int i;
	/// Brief: This is the j index of the cell in which the particle resides currently.
	int j;
	/// Brief: This is the k index of the cell in which the particle resides currently.
	int k;
	/// Brief: It is pointer which stores the x-particle position and updates after every time step 
	double *Xp;
	/// Brief: It is pointer which stores the y-particle position and updates after every time step 
	double *Yp;
	/// Brief: It is pointer which stores the z-particle position and updates after every time step 
	double *Zp;
	///Brief: It is an array which stores u,v,w particle velocity and updates after every time step
	double ParticleVelocity[3];
	///Brief: It is an array which stores u,v,w velocity at actual particle position and updates after every time step
	double VelA[3];
	///Brief: It is an array which stores electric field in x,y,z direction at actual particle position and updates after every time step
	double EfieldA[3];
	///Brief: It is an array which stores Drag force in x,y,z direction at particle position and updates after every time step
	double DragForce[3];
	///Brief: It is an array which stores Electric force in x,y,z direction at particle position and updates after every time step
	double ElectricForce[3];
	///Brief: It stores the Modulus value of actual electric field and updates it after every time step
	double Efielda;
	///Brief: It stores the value of particle charge and updates it after every time step
	double ParticleCharge;
	///Brief: It stores the value of saturation particle charge and updates it after every time step
	double SaturationCharge;
	///Brief: It stores the velocity of the cell at inlet in which the particle resides at the start. This velocity is used in the collection efficiency calculations.
	double Uactual;
}Position;

/** \struct Position
	This is Position structure, it is defined in BaseParameters. It is used for particle tracking,
	where in the position of the particle, ie, the zone in which the particle lies and the indices i,j,k
	of the particular cell in which the particle lies are stored. Also the particles X,Y,Z coordiantes
	are stored.
 */
 
typedef struct
{
	/// Brief: It holds the flag for time step calculation. user can give his own choice of time or can ask the code to calculate the time step using CFL criterion. 0->code calculates time step, 1->user has to enter time step
	int TstepChoice;
	/// Brief: It holds Time stepping flag. 0->BTCS (Backward in time center in space), 1->CrankNicolson (Second order in time)
	int TimeStepingSchemeTag;
	/// Brief: It holds flag of calculating the Vorticity in the computational domain. 0->Not to calculate Vorticity 1->Calculate vorticity
	int tagVorticity;
	/// Brief: It stores the time step for current iteration. 
	double TimeStep;
	/// Brief: It stores the physical time of the problem. It is basically used in the unstedy problems. 
	double Time;
	/// Brief: It stores the flag to indicate whether the solution started afresh or from partially converged solution
	int restart;
        /// Brief: It stores the flag to indicate whether the solver reads initialdata.cgns file
	int read_initialdata_cgns;
	/// Brief: user can fix the Pressure poission convergence criteria, PressLoopConv stores this convergence criteria. It can be later changed from the file /MAIN/PressLoopConv.dat file.
	double PressLoopConv;
	/// Brief: user can fix Mass flux convergence criteria while giving input, it is stored in MassFluxConv.
	double MassFluxConv;
	/// Brief: It stores the ratio of central to upwind convection scheme. 0->Cenrtal convection scheme 1->Upwind convection Scheme.
	double RatioCenUP;
	/// Brief: It stores the ratio of central to upwind convection scheme in equations for TURBULENCE quantities. 0->Cenrtal convection scheme 1->Upwind convection Scheme.
	double RatioCenUPTurb;
	/// Brief: It store value of gravity in X-direction. 
	double gravityX;
	/// Brief: It store value of gravity in Y-direction. 
	double gravityY;
	/// Brief: It store value of gravity in Z-direction. 
	double gravityZ;
        double u_source;
        double v_source;

	/// Brief: It stores the Azimuthal angle for the implementation of Axi-Symettric angle. 
	double *AxisymAzimuthangle;
	/// Brief: It stores the flag of convection scheme (Deffered Correction). 0->Central Convection, 1->Quick convection 
	int ConvectionScheme;
	/// Brief: It stores flag of 2-D case. 0->3-D, 1->2-D
	int flag2d;
	/// Brief: It stores flag for gravity term. It is used as source term in the N-S equation. 0->No gravity, 1->Gravity term
	int flagGravity;
	/// Brief: It stores flag for Buoyancy term. It is used as source term in the N-S equation. 0->No Buoyancy, 1->Buoyancy term
	int flagBuoyancy;
	///Brief: It stores the thermal expansion coefficeint.(volumetric)
	double VolumetricExpCoeff;// Coefficient of Thermal Expansion.
	///Brief: Reference temperature to calculate bouyancy term in Boussinesq approximation.
	double TReference;// Reference temeprature.
	/// Brief: It stores flag for Energy source term. It is used as source term in the Energy equation. 0->No Energy term, 1->Energy term
	int flagEnergySource;
	/// Brief: It stores flag for Dissipation source term. It is used as source term in the N-S equation. 0->No Dissipation term, 1->Dissipation termc
	int flagDissipation;
	/// Brief: It is a pointer, it stores convergence criteria for all flow variable (i.e u, v, w), temperature and 10 sclar.
        int flagJouledissipation;
	double *convergenceCriterion;/*element 0->2 refer to u,v,w ----3 refers to temperature,----and 4->13 refer to the 10 scalars*/
	/// Brief: Differnt variables for Solidification and Melting problem.
	/// Brief: '_s' indicated solid property just before melting and '_f' indicates the liquid property just before solidification
	double K_s,K_f,Rho_s,Rho_f,CP_s,CP_f,E_sm,E_fm,Mu_f,Mu_s,C,P_o,Uc_l,Uc_v,Pc_l,Pc_v;//Created by bhans
	/// Brief: Stores the Melting temperature for isothermal solidification/melting problems
	double Melt_Temp;
	/// Brief: Stores the flag for mushy region problems
	int flag_mushy;
	/// Brief: Stores Liquidus and Solidus temperatures for mushy region problems
	double Solidus_Temp,Liquidus_Temp;//bhans

        //Combustion Variables  Shyam
	/// Brief: It stores Mol.Wts. of species.
        double	 *MolWt; //Shyam,
	/// Brief: It stores Specific heat of species at reference temperature.
        double	 *Cp_ref; //nitesh
	
	/// Brief: It stores no. of moles of each species in the balanced chemical reaction
        double	 *Moles; //Shyam,
	/// Brief: It stores   No. of species participates in the reaction
        int      Nspecies; //Shyam,
	/// Brief: It stores  mass stoichometric ratio of oxidant to fuel.
        double   *Stoich; //Shyam,
	/// Brief:  Cp is a second order function of temp, CpCoeff stores values of coefficients of that function.
        double   *CpCoeff; //Shyam,
	/// Brief:   It stores Fuel mass fraction in the fuel stream, it is defined in diffusion flames
        double   *FuelMassFracInlet; //Shyam,
	/// Brief:   It stores Oxidant mass fraction in the oxidant stream, it is defined in diffusion flames.
        double   *OxiMassFracInlet; //Shyam,
	/// Brief:  To store initial mass fractions of species.
       //double   *InitialMassFracs; //Shyam,
	/// Brief:  For storing initial temperature in the domain.
        double   *InitialTemp; //
	/// Brief:   It stores Heat of the reaction
        double   *delH; //Shyam,
	/// Brief:  It stores Stoichometric mixture fraction
        double   *Zst; //Shyam,
	/// Brief: It Stores enthalpy at fuel inlet
        double   *enthalpyFuelInlet; //Shyam,
	/// Brief:  It Stores enthalpy at oxygen inlet
        double   *enthalpyOxiInlet; // Shyam,
	/// Brief:  It Stores Ambient temperature density
        double   *rhoinfinity; //Shyam,
	/// Brief: It stores inert gas mass fraction values at fuel inlet
	double   *InertgasFuelInlet[10]; //Shyam
	/// Brief: It stores inert gas mass fraction values at oxygen inlet
        double   *InertgasOxiInlet[10]; //Shyam
	/// Brief: It stores number of inert gases present in the system
	int      NIngases;
	/// Brief: It stores the flag for determining whether the fule is a single one or mixture of fuels
        int	 flagfueltype; //0->pure fuel;1->mixture of fuels
	
	// The following are the ELECTRICAL parameters
	/// Brief: It holds Time stepping flag for electric field equations. 0->BTCS (Backward in time center in space), 1->CrankNicolson (Second order in time)
	int EfieldTSSTag;
	/// Brief: It stores the flag of convection scheme for electric field equations, (Deffered Correction). 0->Central Convection, 1->Quick convection 
	int EfieldConvScheme;
	/// Brief: It stores the convergency criteria for electric field equations.
	double EfieldConvCri;
	/// Brief: It holds the flag for time step calculation for electric field equations.
	int EfieldTstepChoice;
	/// Brief: It stores the time step for current iteration for electric field equations.
	double EfieldTimeStep;
	/// Brief: It stores the ratio of central to upwind convection scheme for electric field equations.
	double EfieldRatioCenUP;
	/// Brief: Used for the source term of momentum equations, 1 -> Non-dimensional, 1000 -> Dimensional*/
	double	ScaleFactor;
	/// Brief: It holds flag of calculating the Electric field in 3-directions (Ex, Ey, Ez). 0->Not to calculate Electric field 1->Calculate Electric field
	int TagEfield;
	
	// The following are the LMMHD parameters (Implemented by naveen)
	/// Brief: user can fix the PotentialMHD poisson convergence criteria, PotentialMHDLoopConv stores this convergence criteria. It can be later changed from the file /MAIN/Settings/potentialMHDconvergCriteria.dat file.
	double PotentialMHDLoopConv;
	///Brief: If user wants to apply thin wall condition then this is used to store the value of wall conductance ratio.
	double InducedMHDLoopConv;//rajendra
	double PsiMHDLoopConv;
	double wallconductanceratio;
// 	Groundwater flow sanz
	///it store the flag for 2-D domain in datum plane or in other plane 
	int flagdat;
	/// Brief:It stores the value of convergence head used in groundwater flow
	double  HeadLoopConv;
	/// Brief:It stores the time after which computations should stop. This is used in Ground Water flows.
	double RequiredTime;

	double*OMEGA;
	int flag_Potentialflow;
	//Added by sushil for multiphase flow
	///Brief: It holds the flag for Density to be constant or variable for phases. 0->constant Density; 1->variable Density
	int      flagRhoPhase[5];		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: It holds the flag for Dynamic Viscocity to be constant or variable for phases. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeffPhase[5];
	///Brief:It holds the constant value of Drag coefficient
	double ConstantCD;
	int  Bubbleapprox;
	int Dragflag;
	/// Brief: It stores the ratio of Diamond to Step scheme. 1->Diamond scheme 0->Step Scheme.
	double RatioCenUpRadiation;
	/// Brief: Choice for Radiative Model
	int flagRadiation;
	///Brief: Choice for energy equation with or without Radiation
	int flagflowRadiation;
	/// Number of Points in Theta and Phi direction
	int ntheta, nphi;
	///Convergence Criterion for Radiation Model
	double RadConCriterion;
	/// Decide either radiative flux or temperature.
	int flagRadiativeTemperature;
	///Decide wheather RTE is solved for Absorption
	int flagAbsorption;
	///Decide weather RTE is solved for scattering
	int flagScattering;
	/// Source If applicable in the domain or not
	int flagRadSource;
	/// Inclusion of Fluid Temperature
	int flagTemp;
	///Scattering Function Choice
	int Scatteringchoice;
	
	//Added by Adhiraj for disperse phase flow
	///Brief: It holds the flag for Density to be constant or variable for Disperse phase. 0->constant Density; 1->variable Density
	int      flagRhoDisperse;		/* [S] Flag storing 0->constant Rho; 1->variable Rho */
	/// Brief: It holds the flag for Dynamic Viscocity to be constant or variable for Disperse phase. 0->constant Dynamic Viscosity; 1->variable Dynamic Viscosity
	int      flagViscoeffDisperse;
	//flag for diameter of dispersed particles to be constant(0) or variable(1)
	int      flagVariableDia;
	
	int flagDropFuelMassFrac;
	
	double CPvapour;
	
	double EPSILONDrop;
	//end of disperse phase
	
	// The following are the Particle Tracking parameters
	/// Brief: Stores the number of particles to be tracked.
	int No_ofParticles;
	/// Brief: It stores the number of particles which are leaving the outlet without being collected on ESP plates.
	int ParticlesOut;
	/// Brief: It stores the diameter of the particles.
	double Diameter_ofParticles;
	/// Brief: It stores the mass of the particles.
	double Mass_ofParticles;
	/// Brief: It stores the Relative Permitivity of the medium.
	double RelativePermitivity;
	/// Brief: It stores the density of the particles.
	double Density_ofParticles;
	/// Brief: It stores the timestep used in particle tracking algorithm.
	double deltat;
	/// Brief: It stores the total lenght of the ESP.
	double ESPLength;
	/// Brief: It stores the spacing between wire and the collecting plates.
	double WirePlateSpacing;
	/// Brief: It stores the Collection Efficiency.
	double CollectionEfficiency;
	/// Brief: Variable used in collection efficiency calculation.
	double Weightage;
	/// Brief: Variable used in collection efficiency calculation.
	double TotalWeightage;
	/// Brief: Structure which stores the particle position,(Zone in which particle resides, the i,j,k index and particle x,y,z,position).
	Position *ParticlePosition;
	
	///All parameters for turbulence are stored in the structure "ParametersTurb".  If turbulence is not choosen memory will not be allocated to this structure thus saving memory required for all variables in this structure.
	ParametersTurb *Turb;

	//Added by Mrunalini
	double InitialVolumeFraction;
	double InitialDropletDiameter;

}Parameters;
/** \struct Parameters
	This is BaseParameters structure. It is defined in Base. It contain information 
	about the parameters and various flags for solving the problem in the computational domain.
*/

/* ** ----------------------------- Data Structure for Base ------------------------------------------- ** */
typedef struct
{
	/// Brief: It store the base name. while reading the CGNS file code read the name of base and store here.
	char	basename[33];		/* [C] Name of the 'Base' Node */
	/// Brief: It store the dimension of the problem. i.e 3-D or 2-D.
        int	Cell_Dimension;		/* [C] Cell Dimension */
	/// Brief: It give pointer to the BaseParameters.
	Parameters *BaseParameters;
	/// Brief: It store the dimension of the Grid which we have generated. i.e 3-D or 2-D.
        int	Physical_Dimension;	/* [C] Physical Dimension */
	/// Brief: It gives the information about the total number of Zones in a Base.
        int	Number_of_zones;		/* [C] Number of 'Zone' nodes under the 'Base' node */
	/// Brief: It hold the information about the solver scheme i.e. Implicit, Semi-Implicit and PseudoCompressiblity.
	int	SolverScheme;
	/// Brief: It hold information about the coupling of equations i.e. Coupled solver, semi-Coupled solver and Segregated solver.
	int 	couplingTag;
	/// Brief: It gives pointer of Zone in the base of computational domain.
	Zone   *CgnsZone;		/* [C] 'Zone' Nodes under the Base Node */
	/// Brief: The code calculates the total number of cell in the computational domain and store in TotalCellDomain. It contains both solid and fluid cell in the domain.
	int	TotalCellDomain;
	/// Brief: The code calculaP_o,Uc_l,Uc_v;tes the total number of Solid cell in the computational domain and store in TotalCellSolid.
	int	TotalCellSolid;
	/// Brief: The code calculate the total number of fluid cell in the computational domain and store in TotalCellfluid. 
	int	TotalCellFluid;
	/// Brief: It stores the flag of Steady state for Steady or Unsteady state. i.e. 0->SteadyState 1->UnSteady Problem
	int	flagSteadyState;/*0->SteadyState 1->UnSteady Problem*/
	/// Brief: It stores the flag of Conjugate Problem. i.e. 0->ConjugateProblem 1->Normal Flow Problem
	int	flagConjugateProblem;/*0->ConjugateProblem 1->Normal Flow Problem*/
	/// Brief: It stores the flag for VariableDensityProblem.
	int	flagVariableDensityProblem;
	/// Brief: It stores the flag for type of VariableDensityProblem 0--Temperature variation, 1 Fluid mixing.
	int	flagtypeVariableDensity;
	/// Brief: It stores the flag for Turbulence,i.e; whether turbulence problem has to be solved or not.
	int	flagTurb;
 	/// Brief: It stores the flag for Turbulence Model eg. k-epsilon  model, k-omega model etc.
	int	flagTurbModel;
	/// Brief: It stores the flag for solving logarithmic K and logarithmic Epsilon equations in turbulence modelling
	int	flagLogKEps;
	/// Brief: It stores the flag for Solidification/Meltng problem.  1-> Solidification problem 0-> Not a Solidification problem
	int     flagSolidification;//created by bhans. 
	/// Brief: It stores the flag for Two Phase Flow Boiling problems. 1-> TwoPhase problem 0-> Not a twophase problem
	int     flagTwoPhase;// created by anji
        //Combustion Tags
	/// Brief: Flags for deciding combustion solving, method for turbulent combustion.
        int	flagTurbComb,flagDiffKinetCont; //Shyam, 
	///Brief: Flag for deciding premixed or non-premixed type of combustion.
        int	flagPreNonPre; //Shyam,
	/// Brief: It stores the total number of variables being solved.
	int     total_var;
	/// Brief: It stores "Variable Tags" of each solvable variable in order of initialization.
	int     VARTAG[30];	
	// The following are the ELECTRICAL parameters	
	/// Brief: It stores the flag for electric field,i.e; whether electric field problem has to be solved or not.
	int 	flagEfield;             /*1->Electrically coupled, 0-> Normal flow*/
	/// Brief: It hold information about the type of coupling between flow and electric field
	int 	FtoEcouplingType; /*1->One-Way Coupling 2->Two-Way Coupling*/
	/// Brief: It hold information about the coupling of electric field equations i.e. Coupled solver, semi-Coupled solver and Segregated solver.
	int 	EcoupTag;
	/// Brief: It hold the information about the electric field solver scheme i.e. Implicit, Semi-Implicit and PseudoCompressiblity.
	int	EsolScheme;
	/// Brief: It hold the information about the charge density, i.e whether it is converged or not while implementing the boundary condition for wire
	int ChargeConv;		 /*1->Charge Not converged, 0-> Charge Converged*/
	
	//The following are the LMMHD parameters created by naveen.
	/// Brief: It stores the flag for LMMHD problem.  1-> LMMHD problem 0-> Not a LMMHD problem
	int     flagLMMHD;
	int     flaginduced;//rajendra
	int     flagWf;
// 	Groundwater flow parameters sanz
	/// Brief: It stores the flag for groundwater flow problem.  1-> Groundwater problem 0-> Not a Groundwater problem
	int    flagGroundwater; //used
	/// Brief: It stores the flag for groundwater flow problem.  1-> unsaturated 0-> saturated region problem
	int    flagSaturation; //used
	/// Brief: It stores the flag for groundwater flow problem.  1-> deformable porous media 0-> Not a deformable	
	int    flagdeformable; //used
	///Brief: It store flag whether to solve contaminant transport from given velcoity field or calculate velocity absed on pressure head distribution.
	int flagVelocity;
	
	//Added by sushil..... For multiphase flow variables.
	
	/// Brief: It stores the flag for MPFV problem.  1-> MPFV problem 0-> Not a MPFV problem
	int     flagMPFV;
	/// Brief: It stores the flag for VariableDensityProblem for each phase.
	int	flagVariableDensityPhaseProblem;
	/// Brief: It stores no of phase to be solved.
	int   nphase;

//	-----------------------------------------------------------
//	--------Parameters required for Multigrid technique--------
//	-----------------------------------------------------------
	///Brief: MultiGridflag : for Gauss seidel or Multigrid method selection.
	int MultiGridflag;
	///Brief: VWcycleflag: variable governing V or W cycle.
	int VWcycleflag;
	///Brief: FMG : for Full Multigrid Method.
	int FMG;
	///Brief: CoarseGridIndex stores the no of levels in a multigrid cycle.
	int CoarseGridIndex; //CGInd
//	-----------------------------------------------------------
	///Brief: inp is a structure which stores all input parameters to be scanned from bin/Settings/InputParameters.dat.
	Input inp;
	monitors * Monitors;
	
	/// Brief: stores the number of solvable variables for radiation.
	int	Radvar;
	
	//Added by adhiraj..... For disperse flow variables.
	/// Brief: It stores the flag for Disperse phase problem.  1-> disperse phase problem 0-> Not a disperse phase problem
	int     flagParticle;
	int     flagDroplet;
	///// Brief: It stores the flag for VariableDensityProblem for disperse phase.
	//int	flagVariableDensityPhaseProblem;
	///// Brief: It stores no of phase to be solved.
	//int   nphase;
	//end of disperse phase

	/// Brief: stores the tag for New GS Algo. >sandeep
	int     NewGS;
	///Brief: Bi-CgStab flag.
	int BiCgflag;
	
	/// Brief: It stores the flag for particle tracking,i.e; whether particle tracking is done or not.
	int flagParticle_tracking;             /*1->Particle tracking, 0-> No particle tracking*/
	///Brief: It stores 0-> for constant drag value and 1->user defined drag value.
	int constant;
	//Added by chandrasekhar --- for Turbulent MPFV
	///Brief: It stores tag to update the boundaries in case of Turbulent MPFV
	int masstag;
	///Brief: It stores flag for Turbulent dispersion term in case of Turbulent MPFV
	int flagDispersion;
	///Brief: It stores flag to include Lift force in case of Turbulent MPFV
	int flagLift;
	///Brief: It stores value of turbulent prandtl number for dispersion term in case of Turbulent MPFV
	double number;
	//End --- chandrasekhar
	///Brief: It stores flag for Turbulent modulation term in case of Turbulent Gas-Particle flow.
	int flagModulation;
	
} Base;
/** \struct Base
	This is CgnsBase Structure. It is define in Root structure. The computational domain can be divided 
	into a number of Bases, but normaly we have only one base in a Root. 
*/


/* ** ----------------------------- Data Structure for Root ------------------------------------------- ** */
typedef struct
{
	/// Brief: It stores the input grid filename.
	char inputgridname[333];//[S]
	/// Brief: It stores the output solution filename.
	char outputsolutionname[333];//[S]
	/// Brief: It stores CGNS Version from the CGNS library. 
	float  version;			/* [C] CGNS MLL Version */ 
	/// Brief: It gives the total number of bases in the Root.
	int    Number_of_bases;		/* [C] Number of 'Base' nodes under 'Root' Node */
	/// Brief: It gives pionter to CgnsBase. It gives way to each of the bases in the computational domain.
	Base   *CgnsBase;		/* [C] 'Base' Nodes under the 'Root' Node */
} Root;
/*! \struct Root
	This is CgnsRoot Structure. The computational domain has only one Root. It is the top most stucture 
	in the tree. 
*/
#endif


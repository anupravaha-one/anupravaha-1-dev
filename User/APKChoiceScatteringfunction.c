#include "APKChoiceScatteringfunction.h"

/**
These are user defined functions in which user can write functions in terms of equations for scattering phase function.
*/

double ChoiceScatteringfunction(int choice, int cellindex, double *IntOld, int l, int m, int nt, int np)
	{
	 int li,mi;
	 double value=0.0, PI, theta, phi, theta_i, phi_i, dtheta, dphi;
	 
	 PI=4.0*atan(1.0);
	 dtheta=PI/nt;
	 dphi=2.0*PI/np;
	 
	 switch(choice)
	 	{
		 case 1://Isotropic Scattering
		 
		 	for(li=0;li<nt;li++)
				{
				 theta_i=dtheta*(li+0.5);
				 theta=dtheta*(l+0.5);
				for(mi=0;mi<np;mi++)
					{
					 phi_i=dphi*(mi+0.5);
					 phi=dphi*(m+0.5);
					value=value+*(IntOld+cellindex*nt*np+np*li+mi)*(2.0*dphi*sin(theta_i)*sin(dtheta/2.0))*(2.0*dphi*sin(theta)*sin(dtheta/2.0));
					}
				}
		 break;
		 case 2://Linear forward Anisotropic Scattering
		 	for(li=0;li<nt;li++)
				{
				 theta_i=dtheta*(li+0.5);
				 theta=dtheta*(l+0.5);
				for(mi=0;mi<np;mi++)
					{
					 phi_i=dphi*(mi+0.5);
					 phi=dphi*(m+0.5);
					value=value+*(IntOld+cellindex*nt*np+np*li+mi)*(2.0*dphi*sin(theta_i)*sin(dtheta/2.0))*(2.0*dphi*sin(theta)*sin(dtheta/2.0))+*(IntOld+cellindex*nt*np+np*li+mi)*(cos(phi_i)*sin(dphi/2.0)*(dtheta-cos(2.0*theta_i)*sin(dtheta))*cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta))+sin(phi_i)*sin(dphi/2.0)*(dtheta-cos(2.0*theta_i)*sin(dtheta))*sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta))+dphi*sin(theta_i)*cos(theta_i)*sin(dtheta)*dphi*sin(theta)*cos(theta)*sin(dtheta));
					}
				}
		break;
		case 3://Linear Backwork Anisotroic Scattering
			for(li=0;li<nt;li++)
				{
				 theta_i=dtheta*(li+0.5);
				 theta=dtheta*(l+0.5);
				for(mi=0;mi<np;mi++)
					{
					 phi_i=dphi*(mi+0.5);
					 phi=dphi*(m+0.5);
					value=value+*(IntOld+cellindex*nt*np+np*li+mi)*(2.0*dphi*sin(theta_i)*sin(dtheta/2.0))*(2.0*dphi*sin(theta)*sin(dtheta/2.0))-*(IntOld+cellindex*nt*np+np*li+mi)*(cos(phi_i)*sin(dphi/2.0)*(dtheta-cos(2.0*theta_i)*sin(dtheta))*cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta))+sin(phi_i)*sin(dphi/2.0)*(dtheta-cos(2.0*theta_i)*sin(dtheta))*sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta))+dphi*sin(theta_i)*cos(theta_i)*sin(dtheta)*dphi*sin(theta)*cos(theta)*sin(dtheta));
					}
				}
		break;
		}
	 return(value);
	 }

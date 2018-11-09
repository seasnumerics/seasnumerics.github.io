#include "lbm.hh"

/* Initialize a fixed density inlet and outlet which
 * will evolve into a steady Poiseuille flow. */
void lbm::poiseuille() {
	double Ly=static_cast<double>(ny-1),Lx=static_cast<double>(nx-1);
	double dp=36*Re*nu*nu*Lx/D/Ly/Ly;
	int i,j;
	// bulk
	for(i=1;i<ny-1;i++) {
		for(j=1;j<nx-1;j++) {
			if(!b[i*nx+j]) {
				f[i*nx+j]->set_macro(1.,0,0);
				f[i*nx+j]->set_eq();
			}
		}
	}
	
	for(i=0;i<ny;i++) {
		// east edge and corners - outlet, constant pressure
		if(!b[i*nx+nx-1]) {
			f[i*nx+nx-1]->set_macro(1-dp/2,0,0,true,false);
			f[i*nx+nx-1]->set_eq();
		}
		// west edge and corners - inlet, constant pressure
		if(!b[i*nx]) {
			f[i*nx]->set_macro(1+dp/2,0,0,true,false);
			f[i*nx]->set_eq();
		}
	}

	for(j=1;j<nx-1;j++) {
		// top edge - constant velocity
		if(!b[j]) {
			f[j]->set_macro(1.,0,0,false,true);
			f[j]->set_eq();
		}
		// bottom edge - constant velocity
		if(!b[(ny-1)*nx+j]) {
			f[(ny-1)*nx+j]->set_macro(1.,0,0,false,true);
			f[(ny-1)*nx+j]->set_eq();
		}
	}
	
	// Print a summary of the simulation input.
	printf("Reynolds number: %.0f\n",Re);
	printf("Relaxation parameter: %.6f\n",tau);
	printf("Kinematic viscosity: %.6f\n",nu);
	printf("Characteristic length: %.0f\n",D);
	printf("Initial density at inlet: %.6f\n",f[ny/2*nx]->rho);
	printf("Initial density at outlet: %.6f\n",f[ny/2*nx+nx-1]->rho);
	printf("Centerline velocity: %.6f\n",3*Re*nu/2/D);
}

/* Initialize a steady state Poiseuille flow with fixed velocity profile
 * at the inlet, and constant density at the outlet. */
void lbm::pois_sst() {
	double p,u,y,Ly=static_cast<double>(ny-1),Lx=static_cast<double>(nx-1);
	double dp=36*Re*nu*nu*Lx/D/Ly/Ly;
	int i,j;
	// bulk
	for(i=1;i<ny-1;i++) {
		y=1-1./Ly*i;
		u=6*Re*nu/D*(y-y*y);
		for(j=1;j<nx-1;j++) {
			p=1+dp*(1./2-j/Lx);
			if(!b[i*nx+j]) {
				f[i*nx+j]->set_macro(p,u,0);
				f[i*nx+j]->set_eq();
			}
		}
	}
	
	for(i=0;i<ny;i++) {
		y=1-1./Ly*i;
		u=6*Re*nu/D*(y-y*y);
		// east edge and corners - outlet, constant pressure
		if(!b[i*nx+nx-1]) {
			f[i*nx+nx-1]->set_macro(1-dp/2,u,0,true,false);
			f[i*nx+nx-1]->set_eq();
		}
		// west edge and corners - inlet, constant velocity profile
		if(!b[i*nx]) {
			f[i*nx]->set_macro(1+dp/2,u,0,false,true);
			f[i*nx]->set_eq();
		}
	}

	for(j=1;j<nx-1;j++) {
		p=1+dp*(1./2-j/Lx);
		// top edge - constant velocity
		if(!b[j]) {
			f[j]->set_macro(p,0,0,false,true);
			f[j]->set_eq();
		}
		// bottom edge - constant velocity
		if(!b[(ny-1)*nx+j]) {
			f[(ny-1)*nx+j]->set_macro(p,0,0,false,true);
			f[(ny-1)*nx+j]->set_eq();
		}
	}
	
	// Print a summary of the simulation input.
	printf("Reynolds number: %.0f\n",Re);
	printf("Relaxation parameter: %.6f\n",tau);
	printf("Kinematic viscosity: %.6f\n",nu);
	printf("Characteristic length: %.0f\n",D);
	printf("Initial density at inlet: %.6f\n",f[ny/2*nx]->rho);
	printf("Initial density at outlet: %.6f\n",f[ny/2*nx+nx-1]->rho);
	printf("Centerline velocity: %.6f\n",3*Re*nu/2/D);
}

/* Initialize a cavity with no slip walls and bottom, and a lid
 * with constant horizontal velocity. */
void lbm::cavity_flow() {
	double Ly=static_cast<double>(ny-1),u=Re/Ly*nu;
	int i,j;
	// bulk
	for(i=1;i<ny-1;i++) {
		for(j=1;j<nx-1;j++) {
			if(!b[i*nx+j]) {
				f[i*nx+j]->set_macro(1,0,0);
				f[i*nx+j]->set_eq();
			}
		}
	}
	
	for(i=0;i<ny;i++) {
		// east edge and corners - no slip wall
		if(!b[i*nx+nx-1]) {
			f[i*nx+nx-1]->set_macro(1,0,0,false,true);
			f[i*nx+nx-1]->set_eq();
		}
		// west edge and corners - no slip wall
		if(!b[i*nx]) {
			f[i*nx]->set_macro(1,0,0,false,true);
			f[i*nx]->set_eq();
		}
	}

	for(j=1;j<nx-1;j++) {
		// top edge - constant velocity
		if(!b[j]) {
			f[j]->set_macro(1,u,0,false,true);
			f[j]->set_eq();
		}
		// bottom edge - no slip wall
		if(!b[(ny-1)*nx+j]) {
			f[(ny-1)*nx+j]->set_macro(1,0,0,false,true);
			f[(ny-1)*nx+j]->set_eq();
		}
	}
	
	// Print a summary of the simulation input.
	printf("Reynolds number: %.0f\n",Re);
	printf("Relaxation parameter: %.6f\n",tau);
	printf("Kinematic viscosity: %.6f\n",nu);
	printf("Characteristic length: %.0f\n",D);
	printf("Lid velocity: %.6f\n",u);
}

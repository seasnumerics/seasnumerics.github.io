#include "lbm.hh"

int main() {
    char odir[]="cyl_Re80.out"; // output directory
    char buf[]="cyl24.bin";  // file encoding barrier flags

    // Set up simulation parameters.
    double Re=80;      // Reynolds number
    double tau=0.6;    // Relaxation constant
    double D=24;	   // Reference length

    lbm_params par(Re,tau,D);

    // Set up grid properties.
    int nx=601;		   // Channel length
    int ny=101;		   // Channel width	

    lbm fl(par,nx,ny,odir,28,buf);

    // Initialize the type of flow.
    fl.pois_sst();	// Steady state poiseuille velocity profile
    //fl.poiseuille();  // Pressure boundaries evolving to poiseuille
    //fl.cavity_flow();   // lid-driven cavity flow

    // Evolve in time with regular output.
    fl.solve(60000,1001);
}

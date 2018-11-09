#ifndef LBM_HH
#define LBM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "lbm_params.hh"
#include "lattice.hh"

/** A class for the Lattice-Boltzmann method in 2D. */
class lbm : public lbm_params {
	public:
		// The number of grid points in the x direction.
		const int nx;
		// The number of grid points in the y direction.
		const int ny;
		// Number of threads to use in solve.
		const int nt;
		// Pointers to lattice node objects.
		node_base **f;
		/* Boolean array of barrier indicators of the same
		 * size as the simulation domain. */
		bool *b;
		lbm(lbm_params &par,const int nx_,const int ny_,char *outdir_,const int nt_=1,const char *obsfile=NULL);
		~lbm();
		void setup_grid();
		void solve(long int nsteps,int nout);
		void update();
		void collision();
		void streaming();
		void setup_output(char *outdir_);
		void output(int fr);
		void error();
		// Flow types
		void poiseuille();
		void pois_sst();
		void cavity_flow();
	private:
		/** The output directory filename. */
        char *outdir;
		/** The buffer for assembling output filenames. */
        char *outbuf;
        // Pointer to output file.
		FILE *fp;
		// Timing variables.
		double wt0,wt1;
		// Steady state error.
		double sse;
};

#endif

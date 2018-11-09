#include "lbm.hh"

lbm::lbm(lbm_params &par,const int nx_,const int ny_,char *outdir_,const int nt_,const char *obsfile)
    : lbm_params(par), nx(nx_), ny(ny_), nt(nt_),
    f(new node_base*[nx*ny]), b(new bool[nx*ny]) {

    // Set up output directory.
    setup_output(outdir_);

    // Set up obstacle flag array.
    if(obsfile!=NULL) {
		FILE *fb=safe_fopen(obsfile,"rb");
		safe_fread(b,sizeof(bool),nx*ny,fb,"obstacles");
		fclose(fb);
	}
	else {
	    for(bool *bp=b;bp<b+nx*ny;bp++) *bp=false;
	}

	// Set up grid.
	setup_grid();

}

/* Set up a grid of pointers to lattice objects, selecting the appropriate
 * subclass of the node base for each wall. */
void lbm::setup_grid() {

    int i,j,k;
	for(i=1;i<ny-1;i++) {
		for(j=1;j<nx-1;j++) {
			k=i*nx+j;
			if(b[k]) f[k]=new obstacle();
			else f[k]=new bulk_node(tau); // bulk
		}
	}
	for(i=1;i<ny-1;i++) {
		k=i*nx+nx-1;
		if(b[k]) f[k]=new obstacle();
		else f[k]=new east_edge(tau); // east edge
	}
	for(j=1;j<nx-1;j++) {
		k=j;
		if(b[k]) f[k]=new obstacle();
		else f[k]=new north_edge(tau); // north edge
	}

	for(i=1;i<ny-1;i++) {
		k=i*nx;
		if(b[k]) f[k]=new obstacle();
		else f[k]=new west_edge(tau); // west edge
	}
	for(j=1;j<nx-1;j++) {
		k=(ny-1)*nx+j;
		if(b[k]) f[k]=new obstacle();
		else f[k]=new south_edge(tau); // south edge
	}

	if(b[nx-1]) f[nx-1]=new obstacle();
	else f[nx-1]=new ne_corner(tau); // northeast corner
	if(b[0]) f[0]=new obstacle();
	else f[0]=new nw_corner(tau); // northwest corner
	if(b[(ny-1)*nx]) f[(ny-1)*nx]=new obstacle();
	else f[(ny-1)*nx]=new sw_corner(tau); // southwest corner
	if(b[ny*nx-1]) f[ny*nx-1]=new obstacle();
	f[ny*nx-1]=new se_corner(tau); // southeast corner

	// Populate neighbors.
	// bulk
	for(i=1;i<ny-1;i++) {
		for(j=1;j<nx-1;j++) {
			k=i*nx+j;
			f[k]->set_neighs(f[k+1],f[k-nx],f[k-1],f[k+nx],
						f[k-nx+1],f[k-nx-1],f[k+nx-1],f[k+nx+1]);
		}
	}
	// east edge
	for(i=1;i<ny-1;i++) {
		k=i*nx+nx-1;
		f[k]->set_neighs(NULL,f[k-nx],f[k-1],f[k+nx],
						NULL,f[k-nx-1],f[k+nx-1],NULL);
	}
	// north edge
	for(j=1;j<nx-1;j++) {
		k=j;
		f[k]->set_neighs(f[k+1],NULL,f[k-1],f[k+nx],
						NULL,NULL,f[k+nx-1],f[k+nx+1]);
	}
	// west edge
	for(i=1;i<ny-1;i++) {
		k=i*nx;
		f[k]->set_neighs(f[k+1],f[k-nx],NULL,f[k+nx],
						f[k-nx+1],NULL,NULL,f[k+nx+1]);
	}
	// south edge
	for(j=1;j<nx-1;j++) {
		k=(ny-1)*nx+j;
		f[k]->set_neighs(f[k+1],f[k-nx],f[k-1],NULL,
						f[k-nx+1],f[k-nx-1],NULL,NULL);
	}
	// northeast corner
	k=nx-1;
	f[k]->set_neighs(NULL,NULL,f[k-1],f[k+nx],
						NULL,NULL,f[k+nx-1],NULL);
	// northwest corner
	k=0;
	f[k]->set_neighs(f[k+1],NULL,NULL,f[k+nx],
						NULL,NULL,NULL,f[k+nx+1]);
	// southwest corner
	k=(ny-1)*nx;
	f[k]->set_neighs(f[k+1],f[k-nx],NULL,NULL,
						f[k-nx+1],NULL,NULL,NULL);
	// southeast corner
	k=ny*nx-1;
	f[k]->set_neighs(NULL,f[k-nx],f[k-1],NULL,
						NULL,f[k-nx-1],NULL,NULL);

}

/** The class destructor frees the dynamically allocated memory. */
lbm::~lbm() {
	delete [] b;
	for(int j=nx*ny-1;j>=0;j--) {
            delete f[j];
        }
    delete [] f;
    delete [] outdir;
}

/* Step the Boltzmann transport equation forward in time.
 * \param[in] nsteps number of timesteps to take.
 * \param[in] nout number of equally-spaced output frames to make. */
void lbm::solve(long int nsteps,int nout) {
	long int skip=nsteps/(nout-1);
	int fr=0;
	// Output the initial frame.
	output(fr);
	fr++;
	wt0=wtime();
	for(long int i=1;i<=nsteps;i++) {
		// Collision step.
		collision();
		// Streaming step.
		streaming();
		// Apply BCs and pdate macroscopic quantities
		update();
		// Output the fields to a file.
		if(i%skip==0) {
			output(fr);
			fr++;
		}
	}
}

/* Run the macroscopic quantity and equilibrium distribution updates
 * for the grid. */
void lbm::update() {
	#pragma omp parallel num_threads(nt)
	{
		#pragma omp for
		for(int i=0;i<ny;i++) {
			for(int j=0;j<nx;j++) {
				if(i==0) {
					if(j>0&&j<(nx-1)) f[i*nx+j]->update();
				}
				else if(i==ny-1) {
					if(j>0&&j<(nx-1)) f[i*nx+j]->update();
				}
				else f[i*nx+j]->update();
			}
		}
	}
	// corners
	f[0]->update();
	f[nx-1]->update();
	f[(ny-1)*nx]->update();
	f[ny*nx-1]->update();
}

/* Run the collision step for the grid. */
void lbm::collision() {
	#pragma omp parallel num_threads(nt)
	{
		#pragma omp for
		for(int i=0;i<ny;i++) {
			for(int j=0;j<nx;j++) {
				f[i*nx+j]->collide();
			}
		}
	}
}

/* Run the streaming step for the grid. */
void lbm::streaming() {
	#pragma omp parallel num_threads(nt)
	{
		#pragma omp for
		for(int i=0;i<ny;i++) {
			for(int j=0;j<nx;j++) {
				f[i*nx+j]->stream();
			}
		}
	}
}

/* Set up the output directory. */
void lbm::setup_output(char *outdir_) {
	size_t l=strlen(outdir_)+1;
	outdir=new char[2*l+32];
    memcpy(outdir,outdir_,sizeof(char)*l);
    outbuf=outdir+l;
    // Create output directory, if it doesn't exist.
	mkdir(outdir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
}

/* Perform output of macroscopic quantities. */
void lbm::output(int fr) {
	// write the macroscopic density and velocities to a file.
	// header line: timestep and system dimensions.
	sprintf(outbuf,"%s/fr_%04d.txt",outdir,fr);
	fp=safe_fopen(outbuf,"w");
	fprintf(fp,"%d %d %d\n",fr,nx,ny);
	// fields: rho, u, v.
	double p,u,v;
	for(int i=0;i<ny;i++) {
		for(int j=0;j<nx;j++) {
			f[i*nx+j]->get_macro(p,u,v);
			fprintf(fp,"%g %g %g\n",p,u,v);
		}
	}
	if(fr>0) {
		error();
		printf("Frame %d [%8.5f s], sse=%g\n",fr,(wt1=wtime())-wt0,sse);
		wt0=wt1;
	}
	else {
		printf("Frame 0\n");
	}
	fclose(fp);
}

/* Measure convergence to steady state as a relative error between
 * successive steps. */
void lbm::error() {
	double rho=0,vel=0,dvel=0,p,u,v,du,dv;
	for(int i=0;i<ny;i++) {
		for(int j=0;j<nx;j++) {
			f[i*nx+j]->get_macro(p,u,v);
			f[i*nx+j]->vel_diff(du,dv);
			rho+=p;
			vel+=sqrt(u*u+v*v);
			dvel+=sqrt(du*du+dv*dv);
		}
	}
	sse=dvel/vel;
}

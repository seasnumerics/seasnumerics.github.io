#include "lattice.hh"

/** Set the macroscopic variables. */
void node_base::set_macro(double p,double ux,double uy,bool cpf,bool cvf) {
    rho=p; u=ux; v=uy;
    uold=u; vold=v;
    cp=cpf; cv=cvf;
}

/** Get the macroscopic variables. */
void node_base::get_macro(double &p,double &ux,double &uy) {
    p=rho; ux=u; uy=v;
}

void node_base::set_neighs(node_base *e,node_base *n,node_base *w,node_base *s,
					node_base *ne,node_base *nw,node_base *sw,node_base *se) {
	fe=e; fn=n; fw=w; fs=s; fne=ne; fnw=nw; fsw=sw; fse=se;
}

/** Set the current distribution functions to their equilibrium values. */
void node_base::set_eq() {
	equilibrium();
    f0=feq0; f1=feq1; f2=feq2; f3=feq3; f4=feq4;
    f5=feq5; f6=feq6; f7=feq7; f8=feq8;
}

/** Return the macroscopic quantities. */
void node_base::vel_diff(double &du,double &dv) {
    du=u-uold;
	dv=v-vold;
}

/** Calculate the equilibrium distribution functions. */
void node_base::equilibrium() {
    double ui=3.*u/c, vi=3.*v/c,
    ui2=9.*u*u/c/c/2., vi2=9.*v*v/c/c/2.,
    uvi=9.*u*v/c/c, uvi2=3.*(u*u+v*v)/c/c/2.;

    feq0=4./9*rho*(1-uvi2);
    feq1=1./9*rho*(1+ui+ui2-uvi2);
    feq2=1./9*rho*(1+vi+vi2-uvi2);
    feq3=1./9*rho*(1-ui+ui2-uvi2);
    feq4=1./9*rho*(1-vi+vi2-uvi2);
    feq5=1./36*rho*(1+ui+vi+ui2+vi2+uvi-uvi2);
    feq6=1./36*rho*(1-ui+vi+ui2+vi2-uvi-uvi2);
    feq7=1./36*rho*(1-ui-vi+ui2+vi2+uvi-uvi2);
    feq8=1./36*rho*(1+ui-vi+ui2+vi2-uvi-uvi2);
}

/* Performs the collision step for a grid point. Assumes that the
equilibrium distribution functions of the input have been computed. */
void node_base::collide() {
    double r=1./tau;
    fi0=f0-r*(f0-feq0);
    fi1=f1-r*(f1-feq1);
    fi2=f2-r*(f2-feq2);
    fi3=f3-r*(f3-feq3);
    fi4=f4-r*(f4-feq4);
    fi5=f5-r*(f5-feq5);
    fi6=f6-r*(f6-feq6);
    fi7=f7-r*(f7-feq7);
    fi8=f8-r*(f8-feq8);
}

/* Computes new macroscopic quantities and 
 * equilibrium particle distributions. */
void bulk_node::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) rho=f0+f1+f2+f3+f4+f5+f6+f7+f8;
    if(!cv) {
    	u=c/rho*(f1-f3+f5-f6-f7+f8);
    	v=c/rho*(f2-f4+f5+f6-f7-f8);
    }
	equilibrium();
}

/* Applies the bounce-back boundary condition. */
void bulk_node::bounce_back() {
	if(fe->c==0) f3=fi1;
	if(fn->c==0) f4=fi2;
	if(fw->c==0) f1=fi3;
	if(fs->c==0) f2=fi4;
	if(fne->c==0) f7=fi5;
	if(fnw->c==0) f8=fi6;
	if(fsw->c==0) f5=fi7;
	if(fse->c==0) f6=fi8;
}

void east_edge::update() { // outlet
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known
		rho=c/(c+u)*(f0+f2+f4+2*(f1+f5+f8));
	}
	if(!cv) {
		// rho is known
		u=c*(1./rho*(f0+f2+f4+2*(f1+f5+f8))-1.); // momentum conservation.
		v=c/rho*(fw->f2-fw->f4+fw->f5+fw->f6-fw->f7-fw->f8); // Neumann-type condition.
	}
	equilibrium();
	f3=feq3+f1-feq1;
	f6=1./2*(rho/c*(v-u)+f1-f2-f3+f4+2*f8);
	f7=1./2*(-rho/c*(v+u)+f1+f2-f3-f4+2*f5);
}
void east_edge::bounce_back() {
	if(fn->c==0) f4=fi2;
	if(fw->c==0) f1=fi3;
	if(fs->c==0) f2=fi4;
	if(fnw->c==0) f8=fi6;
	if(fsw->c==0) f5=fi7;
}

void north_edge::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known
		rho=c/(c+v)*(f0+f1+f3+2*(f2+f5+f6));
	}
	if(!cv) { // rho is known
		v=c*(1./rho*(f0+f1+f3+2*(f2+f5+f6))-1.); // momentum conservation.
		u=c/rho*(fs->f1-fs->f3+fs->f5-fs->f6-fs->f7+fs->f8); // Neumann-type condition.
	}
	equilibrium();
	f4=feq4+f2-feq2;
	f7=1./2*(-rho/c*(u+v)+f1+f2-f3-f4+2*f5);
	f8=1./2*(rho/c*(u-v)-f1+f2+f3-f4+2*f6);
}
void north_edge::bounce_back() {
	if(fe->c==0) f3=fi1;
	if(fw->c==0) f1=fi3;
	if(fs->c==0) f2=fi4;
	if(fsw->c==0) f5=fi7;
	if(fse->c==0) f6=fi8;
}

void west_edge::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known
		rho=c/(c-u)*(f0+f2+f4+2*(f3+f6+f7));
	}
	if(!cv) { // rho is known
		u=c*(1-1./rho*(f0+f2+f4+2*(f3+f6+f7))); // momentum conservation.
		v=c/rho*(fe->f2-fe->f4+fe->f5+fe->f6-fe->f7-fe->f8); // Neumann-type condition.
	}
	equilibrium();
	f1=feq1+f3-feq3;
	f5=1./2*(rho/c*(u+v)-f1-f2+f3+f4+2*f7);
	f8=1./2*(rho/c*(u-v)-f1+f2+f3-f4+2*f6);
}
void west_edge::bounce_back() {
	if(fe->c==0) f3=fi1;
	if(fn->c==0) f4=fi2;
	if(fs->c==0) f2=fi4;
	if(fne->c==0) f7=fi5;
	if(fse->c==0) f6=fi8;
}

void south_edge::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known
		rho=c/(c-v)*(f0+f1+f3+2*(f4+f7+f8));
	}
	if(!cv) { // rho is known
		v=c*(1-1./rho*(f0+f1+f3+2*(f4+f7+f8))); // momentum conservation.
		u=c/rho*(fn->f1-fn->f3+fn->f5-fn->f6-fn->f7+fn->f8); // Neumann-type condition.
	}
	equilibrium();
	f2=feq2+f4-feq4;
	f5=1./2*(rho/c*(u+v)-f1-f2+f3+f4+2*f7);
	f6=1./2*(rho/c*(v-u)+f1-f2-f3+f4+2*f8);
}
void south_edge::bounce_back() {
	if(fe->c==0) f3=fi1;
	if(fn->c==0) f4=fi2;
	if(fw->c==0) f1=fi3;
	if(fne->c==0) f7=fi5;
	if(fnw->c==0) f8=fi6;
}

void ne_corner::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known - currently the only specification.
		rho=fs->rho;
	}
	equilibrium();
	f3=feq3+f1-feq1; f4=feq4+f2-feq2;
	f6=1./2*(rho+rho/c*v-f0-f1-2*f2-f3-2*f5);
	f7=1./2*(-rho/c*(u+v)+f1+f2-f3-f4+2*f5);
	f8=1./2*(rho+rho/c*u-f0-2*f1-f2-f4-2*f5);
}

void ne_corner::bounce_back() {
	if(fw->c==0) f1=fi3;
	if(fs->c==0) f2=fi4;
	if(fsw->c==0) f5=fi7;
}

void nw_corner::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known - currently the only specification.
		rho=fs->rho;
	}
	equilibrium();
	f1=feq1+f3-feq3; f4=feq4+f2-feq2;
	f5=1./2*(rho+rho/c*v-f0-f1-2*f2-f3-2*f6);
	f7=1./2*(rho-rho/c*u-f0-f2-2*f3-f4-2*f6);
	f8=1./2*(rho/c*(u-v)-f1+f2+f3-f4+2*f6);
}

void nw_corner::bounce_back() {
	if(fe->c==0) f3=fi1;
	if(fs->c==0) f2=fi4;
	if(fse->c==0) f6=fi8;
}

void sw_corner::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known - currently the only specification.
		rho=fn->rho;
	}
	equilibrium();
	f1=feq1+f3-feq3; f2=feq2+f4-feq4;
	f5=1./2*(rho/c*(u+v)-f1-f2+f3+f4+2*f7);
	f6=1./2*(rho-rho/c*u-f0-f2-2*f3-f4-2*f7);
	f8=1./2*(rho-rho/c*v-f0-f1-f3-2*f4-2*f7);
}

void sw_corner::bounce_back() {
	if(fe->c==0) f3=fi1;
	if(fn->c==0) f4=fi2;
	if(fne->c==0) f7=fi5;
}

void se_corner::update() {
	bounce_back();
	uold=u; vold=v;
	if(!cp) { // u and v are known - currently the only specification.
		rho=fn->rho;
	}
	equilibrium();
	f3=feq3+f1-feq1; f2=feq2+f4-feq4;
	f5=1./2*(rho+rho/c*u-f0-2*f1-f2-f4-2*f8);
	f6=1./2*(rho/c*(v-u)+f1-f2-f3+f4+2*f8);
	f7=1./2*(rho-rho/c*v-f0-f1-f3-2*f4-2*f8);
}

void se_corner::bounce_back() {
	if(fn->c==0) f4=fi2;
	if(fw->c==0) f1=fi3;
	if(fnw->c==0) f8=fi6;
}

/* Performs the streaming step. */
void bulk_node::stream() {
	f0=fi0; f1=fw->fi1; f2=fs->fi2; f3=fe->fi3; f4=fn->fi4;
	f5=fsw->fi5; f6=fse->fi6; f7=fne->fi7; f8=fnw->fi8;
}

void east_edge::stream() {
	f0=fi0; f1=fw->fi1; f2=fs->fi2; f4=fn->fi4; f5=fsw->fi5; f8=fnw->fi8;
}

void north_edge::stream() {
	f0=fi0; f1=fw->fi1; f2=fs->fi2; f3=fe->fi3; f5=fsw->fi5; f6=fse->fi6;
}

void west_edge::stream() {
	f0=fi0; f2=fs->fi2; f3=fe->fi3; f4=fn->fi4; f6=fse->fi6; f7=fne->fi7;
}

void south_edge::stream() {
	f0=fi0; f1=fw->fi1; f3=fe->fi3; f4=fn->fi4; f7=fne->fi7; f8=fnw->fi8;
}

void ne_corner::stream() {
	f0=fi0; f1=fw->fi1; f2=fs->fi2; f5=fsw->fi5;
}

void nw_corner::stream() {
	f0=fi0; f2=fs->fi2; f3=fe->fi3; f6=fse->fi6;
}

void sw_corner::stream() {
	f0=fi0; f3=fe->fi3; f4=fn->fi4; f7=fne->fi7;
}

void se_corner::stream() {
	f0=fi0; f1=fw->fi1; f4=fn->fi4; f8=fnw->fi8;
}

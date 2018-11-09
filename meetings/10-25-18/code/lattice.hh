#ifndef LATTICE_HH
#define LATTICE_HH

#include <cstdio>
#include <cmath>

/** Base class for storing the particle distribution function at grid points. */
class node_base {
    public:
        /** The lattice speed. */
        const double c;
        /** The relaxation rate. */
        const double tau;
        /** The macroscopic density. */
        double rho;
        /** The macroscopic horizontal velocity. */
        double u;
        /** The macroscopic vertical velocity. */
        double v;
        /** Whether the cell is at constant pressure. */
        bool cp;
        /** Whether the cell is at constant velocity. */
        bool cv;
        /** The particle distribution function with in-place stream. */
        double f0;
        /** The particle distribution function streaming east. */
        double f1;
        /** The particle distribution function streaming north. */
        double f2;
        /** The particle distribution function streaming west. */
        double f3;
        /** The particle distribution function streaming south. */
        double f4;
        /** The particle distribution function streaming northeast. */
        double f5;
        /** The particle distribution function streaming northwest. */
        double f6;
        /** The particle distribution function streaming southwest. */
        double f7;
        /** The particle distribution function streaming southeast. */
        double f8;
        /** The equilibrium particle distribution functions. */
        double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8;
        /** Intermediate values of f, before streaming. */
        double fi0,fi1,fi2,fi3,fi4,fi5,fi6,fi7,fi8;
        /* Tracking previous velocity for convergence to steady state. */
        double uold,vold;
        /** Pointers to lattice neighbors in order. */
        node_base *fe,*fn,*fw,*fs,*fne,*fnw,*fsw,*fse;
        node_base(double c_,double tau_) : c(c_), tau(tau_), cp(false), cv(false) {}
        virtual ~node_base() {}
        void set_macro(double p,double ux,double uy,bool cpf=false,bool cvf=false);
        void get_macro(double &p,double &ux,double &uy);
        void set_eq();
        void set_neighs(node_base *e,node_base *n,node_base *w,node_base *s,
                    node_base *ne,node_base *nw,node_base *sw,node_base *se);
        void vel_diff(double &du,double &dv);
        void equilibrium();
        void collide();
        virtual void update() = 0;
        virtual void stream() = 0;
        virtual void bounce_back() = 0;
};

class bulk_node : public node_base {
    public:
        bulk_node(double tau_) : node_base(1,tau_) {}
        virtual ~bulk_node() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class east_edge : public node_base {
    public:
        east_edge(double tau_) : node_base(1,tau_) {}
        virtual ~east_edge() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class north_edge : public node_base {
    public:
        north_edge(double tau_) : node_base(1,tau_) {}
        virtual ~north_edge() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class west_edge : public node_base {
    public:
        west_edge(double tau_) : node_base(1,tau_) {}
        virtual ~west_edge() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class south_edge : public node_base {
    public:
        south_edge(double tau_) : node_base(1,tau_) {}
        virtual ~south_edge() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class ne_corner : public node_base {
    public:
        ne_corner(double tau_) : node_base(1,tau_) {}
        virtual ~ne_corner() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class nw_corner : public node_base {
    public:
        nw_corner(double tau_) : node_base(1,tau_) {}
        virtual ~nw_corner() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class sw_corner : public node_base {
    public:
        sw_corner(double tau_) : node_base(1,tau_) {}
        virtual ~sw_corner() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class se_corner : public node_base {
    public:
        se_corner(double tau_) : node_base(1,tau_) {}
        virtual ~se_corner() {}
        virtual void update();
        virtual void stream();
        virtual void bounce_back();
};
class obstacle : public node_base {
    public:
        obstacle() : node_base(0,1) {set_macro(0,0,0);}
        virtual ~obstacle() {}
        virtual void update() {}
        virtual void stream() {}
        virtual void bounce_back() {}
};

#endif
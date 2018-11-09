#ifndef LBM_PARAMS_HH
#define LBM_PARAMS_HH

struct lbm_params {
    /** Reynolds number. */
    const double Re;
    /** The relaxation time. */
    const double tau;
    /** kinematic viscosity. */
    const double nu;
    /** Characteristic length. */
    const double D;

    lbm_params(double Re_,double tau_,double D_) : Re(Re_), tau(tau_), nu((2*tau-1)/6.), D(D_) {}
};

#endif

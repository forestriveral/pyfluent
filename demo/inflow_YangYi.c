#include "udf.h"
/* UDF for specifying steady-state velocity profile boundary condition */
/* INFLOW version: Standard k-E inflow profile (Yangyi shear stress) */

/* default parameters references */
/* u_star = 0.511     0.73407         */
/*     z0 = 0.00023   0.000225        */
/*    Cmu = 0.03      0.028           */
/*     C1 = -1.36    -0.55     -0.17  */
/*     C2 = 12.3      5.21     1.62   */
/*  sig_e = 1.3       2.51            */
/* C2e_1e = 0.48      0.42            */

/* fluent inflow settings */
/* RSM */
/* Roughness height:    0.006   */
/* Roughness constant:  0.9     */
/* C_mu (standard):     0.09    */
/* k-e */
/* Roughness height:    0.0025  */
/* Roughness constant:  0.9     */
/* C_mu (standard):     0.12    */
/* C_mu (RNG):          0.88    */
/* k-w */
/* Roughness height:    0.006   */
/* Roughness constant:  0.9     */


#define u_star 0.3256       /* Friction velocity (m/s) */
#define K 0.4187             /* Karman constant */
#define rho 1.225          /* Air density (kg/m3) */
#define z0 0.000221         /* Aerodynamic roughness length (m) */

#define Cmu 0.4000          /* Parameter in velocity TKE and dissipation profile */
#define C1 -0.9769          /* First constant in velocity TKE and dissipation profile */
#define C2 7.2829            /* Second constant in velocity TKE and dissipation profile */

#define sig_e 1.3         /* First constant in w source equation */
#define C2e_1e 0.88        /* Second constant in w source equation */

/* Optional parameters for profile modification */
#define H_hub 0.125         /* Hub height (m) */
#define D_rotor 0.125       /* Rotor diameter (m) */
#define H_star 1.86         /* Modification height ration (m) */


/*PROFILE FOR STREAMWISE VELOCITY*/
DEFINE_PROFILE(inlet_velocity, thread, index)
{
    real x[ND_ND];
    real z;
    face_t f;
    begin_f_loop(f, thread)
    {
        F_CENTROID(x, f, thread);
        z = x[2];
        F_PROFILE(f, thread, index) = (u_star / K) * log((z + z0) / z0);
    }
    end_f_loop(f, thread)
}


/*PROFILE FOR KINETIC ENERGY*/
DEFINE_PROFILE(k_profile, thread, index)
{
    real x[ND_ND];
    real z;
    face_t f;
    begin_f_loop(f, thread)
    {
        F_CENTROID(x, f, thread);
        z = x[2];
        if ((C1 * log((z + z0) / z0) + C2) >= 0) {
            F_PROFILE(f, thread, index) = pow(u_star, 2) / sqrt(Cmu) * sqrt(C1 * log((z + z0) / z0) + C2); /* Yangyi SS */
        }
        else {
            F_PROFILE(f, thread, index) = 0.3 * exp(-1.5 * z / D_rotor) + 0.0;
        }
    }
    end_f_loop(f, thread)
}


/*PROFILE FOR DISSIPATION RATE*/
DEFINE_PROFILE(w_profile, thread, index)
{
    real x[ND_ND];
    real z;
    face_t f;
    begin_f_loop(f, thread)
    {
        F_CENTROID(x, f, thread);
        z = x[2];
        if ((C1 * log((z + z0) / z0) + C2) >= 0) {
            F_PROFILE(f, thread, index) = (pow(u_star, 3) / (K * (z + z0))) * sqrt(C1 * log((z + z0) / z0) + C2);  /* Yangyi SS */
        }
        else {
            F_PROFILE(f, thread, index) = pow(u_star, 3) / (K * (z + z0));
        }
    }
    end_f_loop(f, thread)
}


/*SOURCE TERM FOR TKE EQUATION*/
DEFINE_SOURCE(k_source, c, t, dS, eqn)
{
    real x[ND_ND];
    real source;
    real z;

    C_CENTROID(x, c, t);
    z = x[2];
    source = 0;  /* Yangyi SS */
    dS[eqn] = 0;

    return source;
}


/*SOURCE TERM FOR DISSIPATION RATE EQUATION*/
DEFINE_SOURCE(w_source, c, t, dS, eqn)
{
    real x[ND_ND];
    real source;
    real d, z;

    C_CENTROID(x, c, t);
    z = x[2];
    if ((C1 * log((z + z0) / z0) + C2) >= 0) {
        source = rho * pow(u_star, 4) / (pow(K, 2) * pow(z + z0, 2) * sig_e) * (pow(K, 2) * (1.5 * C1 - C1 * log((z + z0) / z0) - C2) + sqrt(Cmu) * sig_e * C2e_1e * sqrt(C1 * log((z + z0) / z0) + C2));
    }
    else {
        source = rho * pow(u_star, 4) / (pow(K, 2) * pow(z + z0, 2) * sig_e) * (pow(K, 2) * (1.5 * C1 - C1 * log((z + z0) / z0) - C2));
    }
    dS[eqn] = 0;

    return source;
}


DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{
    Message("\n Standard k-e turbulence model for inflow profile (Yi Yang et al.): ");
    Message("u_star= %f, K= %f, rho= %f, z_0= %f ", u_star, K, rho, z0);
    Message("C_mu= %f, C_1= %f, C_2= %f ", Cmu, C1, C2);
    Message("Sigma_e= %f, C2e_1e= %f\n", sig_e, C2e_1e);
}
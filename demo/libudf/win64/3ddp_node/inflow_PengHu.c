#include "udf.h"
/* UDF for specifying steady-state velocity profile boundary condition */
/* INFLOW version: SST k-W inflow profile (Peng Hu et al.) */

/* default parameters references */
/* u_star = 0.511     0.73407   0.73407   0.73407   */
/*     z0 = 0.000225  0.00612   0.00612   0.00612   */
/*    Cu1 = -0.172    -0.33847  -0.25096  -0.15417  */
/*    Cu2 = 3.088     2.61666   2.25328   1.85135   */
/*    C1k = 0.10      0.172     0.172     0.172     */
/*    C2k = 0.64      0.716     0.706     0.696     */
/*    C1w = 0.12      0.156     0.156     0.156     */
/*    C4w = 0.02      0.015     0.015     0.015     */

/* fluent inflow settings */
/* RSM */
/* Roughness height:    0.006   */
/* Roughness constant:  0.9     */
/* k-e */
/* Roughness height:    0.0025  */
/* Roughness constant:  0.9     */
/* C_mu (standard):     0.12    */
/* C_mu (RNG):          0.88    */
/* k-w */
/* Roughness height:    0.006   */
/* Roughness constant:  0.9     */


#define u_star 0.3256        /* Friction velocity (m/s) */
#define K 0.4187             /* Karman constant */
#define rho 1.225          /* Air density (kg/m3) */
#define z0 0.000221        /* Aerodynamic roughness length (m) */

#define Cu1 -0.5235         /* First constant in velocity TKE and dissipation profile */
#define Cu2 4.6026           /* Second constant in velocity TKE and dissipation profile */

#define C1k 0.2720          /* First constant in k source equation */
#define C2k 0.7160           /* Second constant in k source equation */

#define C1w 0.1560           /* First constant in w source equation */
#define C2w 0.0            /* Second constant in w source equation */
#define C3w 0.0            /* Third constant in w source equation */
#define C4w 0.0150           /* Second constant in w source equation */

/* Optional parameters for profile modification */
#define H_hub 0.125         /* Hub height (m) */
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
        F_PROFILE(f, thread, index) = pow(u_star, 2) * pow((Cu1 * log((z + z0) / z0) + Cu2), 2);
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
        F_PROFILE(f, thread, index) = (u_star / (K * (z + z0))) * pow((Cu1 * log((z + z0) / z0) + Cu2), 2);
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
    source = rho * pow(u_star, 3) / (z + z0) * (C1k * pow(Cu1 * log((z + z0) / z0) + Cu2, 4) - C2k);
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
    source = rho * pow(u_star, 2) / pow(z + z0, 2) * (C1w * pow(Cu1 * log((z + z0) / z0) + Cu2, 4) - C4w);
    dS[eqn] = 0;

    return source;
}


DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{
    Message("\n SST k-w turbulence model for inflow profile (Peng Hu et al.): ");
    Message("u_star= %f, K= %f, rho= %f, z_0= %f", u_star, K, rho, z0);
    Message("C_u1= %f, C_u2= %f, C_1k= %f, C_2k= %f", Cu1, Cu2, C1k, C2k);
    Message("C_1w= %f, C_2w= %f, C_3w= %f, C_4w= %f\n", C1w, C2w, C3w, C4w);
}
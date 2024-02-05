#include "udf.h"
/* UDF for specifying steady-state velocity profile boundary condition */
/* INFLOW version: RSM inflow profile (LI Weibo and WANG Guoyan) */

/* default parameters references */
/* u_star = 0.511       */
/*     z0 = 0.000225    */
/*      a = 0.8         */
/*    Cmu = 0.09        */
/*     C1 = 3.088       */
/*     C2 = 0.10        */


#define u_star 0.3256       /* Friction velocity (m/s) */
#define K 0.4187             /* Karman constant */
#define rho 1.225          /* Air density (kg/m3) */
#define z0 0.000221         /* Aerodynamic roughness length (m) */

#define a 0.8              /* Scale constant in velocity TKE and dissipation profile */
#define Cmu 1.09         /* First constant in velocity TKE and dissipation profile */
#define C1 0.288          /* First constant in velocity TKE and dissipation profile */
#define C2 0.10            /* Second constant in velocity TKE and dissipation profile */

// #define C1k 0.005          /* First constant in k source equation */
// #define C2k 0.10           /* Second constant in k source equation */

// #define C1w 1.44           /* First constant in w source equation */
// #define C2w 0.0            /* Second constant in w source equation */
// #define C3w 0.0            /* Third constant in w source equation */
// #define C4w 1.92           /* Second constant in w source equation */

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
        F_PROFILE(f, thread, index) = a * pow(u_star, 2) * pow((C1 * log((z + z0) / z0) + C2), 2);
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
        F_PROFILE(f, thread, index) = pow(a, 2) * Cmu * pow(u_star, 3) * pow((C1 * log((z + z0) / z0) + C2), 4) / (K * (z + z0));
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
    // source = rho * pow(u_star, 3) / (z + z0) * (C1k * pow(Cu1 * log((z + z0) / z0) + Cu2, 4) - C2k);
    source = 0;
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
    // source = rho * pow(u_star, 2) / pow(z + z0, 2) * (C1w * pow(Cu1 * log((z + z0) / z0) + Cu2, 4) - C4w);
    source = 0;
    dS[eqn] = 0;

    return source;
}


DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{
    Message("\n RSM inflow profile (LI Weibo and WANG Guoyan): \n");
    Message("u_star= %f, K= %f, rho= %f, z_0= %f", u_star, K, rho, z0);
    Message("a= %f, C_mu= %f, C_1= %f, C_2= %f", a, Cmu, C1, C2);
}
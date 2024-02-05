#include "udf.h"
/* UDF for specifying steady-state velocity profile boundary condition */
/* INFLOW version: SST k-W inflow profile (Yangyi shear stress) */

#define u0 0.115           /* Friction velocity */
#define k 0.42             /* Karman constant */
#define rho 1.225          /* Air density */
#define z0 0.00003         /* Aerodynamic roughness length */
#define h 0.125            /* Hub height */

#define Cmu 0.09
#define Cu1 -0.412          /* First constant in velocity TKE and dissipation profile */
#define Cu2 3.77           /* Second constant in velocity TKE and dissipation profile */

#define C1k 0.005           /* First constant in k source equation */
#define C2k 0.10           /* Second constant in k source equation */
#define sigma_k 1.0

#define C1w 1.44           /* First constant in w source equation */
#define C2w 1.92           /* Second constant in w source equation */
#define Cw 0.48
#define sigma_e 1.50

#define ht 1.86

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
        F_PROFILE(f, thread, index) = (u0 / k) * log((z + z0) / z0);
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
        // F_PROFILE(f, thread, index) = pow(u0, 2) * pow(Cu1 * log((z + z0) / z0) + Cu2, 2);  /* SST k-W */
        if (z > ht * h) {z = ht * h;}
        F_PROFILE(f, thread, index) = pow(u0, 2) / sqrt(Cmu) * sqrt(Cu1 * log((z + z0) / z0) + Cu2); /* Yangyi SS */
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
        // F_PROFILE(f, thread, index) = (u0 / (k * (z + z0))) * pow(Cu1 * log((z + z0) / z0) + Cu2, 2);  /* SST k-W */
        if (z > ht * h) {z = ht * h;}
        F_PROFILE(f, thread, index) = (pow(u0, 3) / (k * (z + z0))) * sqrt(Cu1 * log((z + z0) / z0) + Cu2);  /* Yangyi SS */
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
    // source = rho * pow(u0, 3) / (z + z0) * (C1k * pow(Cu1 * log((z + z0) / z0) + Cu2, 4) - C2k);  /* SST k-W */
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
    d = x[0];
    z = x[2];
    // source = rho * pow(u0, 2) / pow(z + z0, 2) * (C1w * pow(Cu1 * log((z + z0) / z0) + Cu2, 4) - C2w);  /* SST k-W */
    if (d <= 10.0)
    {
        if (z > ht * h) {z = ht * h;}  /* Yangyi SS */
        source = rho * pow(u0, 4) / (pow(k, 2) * pow(z + z0, 2) * sigma_e) * (pow(k, 2) * (1.5 * Cu1 - Cu1 * log((z + z0) / z0) - Cu2) + sqrt(Cmu) * sigma_e * (Cw) * sqrt(Cu1 * log((z + z0) / z0) + Cu2));
    }
    else
    {
        source = 0;
    }
    dS[eqn] = 0;

    return source;
}


// DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
// {
//     Message("\nCu1= %f, Cu2= %f, Cw= %f, sigma_e= %f", Cu1, Cu2, Cw, sigma_e);
// }
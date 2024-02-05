#include "udf.h"
#include <math.h>
/* UDF version:  Ut and aprime1 inside the loop cycle; No tip/hub loss; disc thickness 0.5 */

/* General parameters for multiple yawed turbine */
#define pi 3.141592653      /* Pi constant */
#define B 3                 /* Number of the blades */   
#define dx 0.015            /* Thickness of the rotor disk (m)*/
#define oft 0.5             /* Offset of thickness of the rotor disk */
#define N_BE 8              /* Number of blade elements divided */
#define NUM_UDM 14          /* Number of UDM for variable storage */
#define Cd_nac 100          /* Drag coefficient of nacell */
#define Cd_tow 100          /* Drag coefficient of tower */

// /* Case parameters for the default one yawed turbine */
// #define N_WT 1              /* Number of turbines (default number is 1) */
// #define theta 0.0           /* Yaw angle of the turbine (degree) */
// #define radius 0.075        /* Radius of the rotor (m) */
// #define pitch 0.0           /* Pitch angle of the blade (degree) */
// #define inflow 2.2          /* The freestream velocity (m/s) */
// #define tsr 4.0             /* Tip-speed ratio of the turbine */
// #define rpm 1120            /* Rotation speed of the turbine (round per minute) */
// #define N_LD 93             /* Number of lift/drag coefficient data of airfoils blade */
// #define scale 1.0           /* Scale factor of simulation */

// /* Coordinate of yawed turbine rotor center if N_WT = 1 */
// #define x0 0.0
// #define y0 0.0
// #define z0 0.125

/* Key simulated parameters for multiple yawed turbine */
#define N_WT 2                           /* Number of multiple yawed turbines (default number is 2) */
real theta_wt[N_WT] = {0.0, 0.0};        /* Yaw angle of multiple yawed turbines (degree) */
real radius_wt[N_WT] = {0.075, 0.075};   /* Rotor radius of multiple yawed turbines (m) */
real pitch_wt[N_WT] = {0.0, 0.0};        /* Blade pitch angle of multiple yawed turbines (degree) */
real inflow_wt[N_WT] = {2.2, 1.8};       /* The effective inflow velocity of multiple yawed turbines (m/s) */
real tsr_wt[N_WT] = {4.0, 4.0};          /* Tip-speed ratio of multiple yawed turbines */
real rpm_wt[N_WT] = {1120, 950};         /* Rotation speed of multiple yawed turbines (round per minute) */
#define N_LD 93                          /* Default number of lift/drag coefficient data of airfoils blade */
#define scale 1.0                        /* Scale factor of simulation */

/* Coordinate array of yawed turbine rotor center */
real x_wt[N_WT] = {0.0, 0.75};           /* X coordinate array of multiple yawed turbines */
real y_wt[N_WT] = {0.0, 0.0};            /* Y coordinate array of multiple yawed turbines */
real z_wt[N_WT] = {0.125, 0.125};        /* Z coordinate array of multiple yawed turbines */

/* Size of multiple yawed turbines nacelle (Same setting for all turbine) */
#define diameter_nac 0.015
#define length_nac 0.03

/* Size and circle center coordinate of multiple yawed turbines tower (Same setting for all turbine) */
#define diameter_tow 0.01
#define height_tow 0.118
#define x_offset_tow 0.03
#define y_offset_tow 0.0
#define z_offset_tow 0.0


static int udm_offset = 0;


DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{
    Message("\n%d UDMs have been reserved by the current library %s", NUM_UDM, libname);

    Set_User_Memory_Name(udm_offset, "f_n");
    Set_User_Memory_Name(udm_offset + 1, "f_t");
    Set_User_Memory_Name(udm_offset + 2, "f_x");
    Set_User_Memory_Name(udm_offset + 3, "f_y");
    Set_User_Memory_Name(udm_offset + 4, "f_z");
    Set_User_Memory_Name(udm_offset + 5, "f_nx");
    Set_User_Memory_Name(udm_offset + 6, "f_ny");
    Set_User_Memory_Name(udm_offset + 7, "f_nz");
    Set_User_Memory_Name(udm_offset + 8, "f_tx");
    Set_User_Memory_Name(udm_offset + 9, "f_ty");
    Set_User_Memory_Name(udm_offset + 10, "f_tz");
    Set_User_Memory_Name(udm_offset + 11, "total_f_n"); 
    Set_User_Memory_Name(udm_offset + 12, "total_f_t"); 

    Message("\nUDM Offset for Current Loaded Library = %d", udm_offset);
}


/* Enumeration of used User-Defined Memory Locations */
enum rotor_force
{
    f_n,
    f_t,
    f_x,
    f_y,
    f_z,
    f_nx,
    f_ny,
    f_nz,
    f_tx,
    f_ty,
    f_tz,
    total_f_n,
    total_f_t,
    N_REQUIRED_UDM
};


DEFINE_EXECUTE_AT_END(caclulate_inflow)
{
#if RP_NODE
    Domain *d;
    d=Get_Domain(1);
    Thread *t;
    cell_t c;
    
    real location[ND_ND];
    real x1,y1,z1,x0t1_mm, y0t1_mm, z0t1_mm,xt1, yt1, zt1, rt1_mm;
    real epsilon = 1e-10;
    real finally_v;
    real v;
    real vel;
    real average_v;
    int mm;
    real x01[N_WT] = {0,382};
    real y01[N_WT] = {0,0};
    real z01[N_WT] = {90,90};
    int ID=2;
    for(mm=0; mm<N_WT; mm++)
    {
        real sum_v=0.0;
        real sum_volume=0.0;
        Thread *t =Lookup_Thread(d,ID);
            begin_c_loop_int(c,t)
            {
                C_CENTROID(location, c, t);
                x1 = location[0];
                y1 = location[1];
                z1 = location[2];
                x0t1_mm = x01[mm] * cos(theta * pi / 180) + y01[mm] * sin(theta * pi / 180);
                y0t1_mm = y01[mm] * cos(theta * pi / 180) - x01[mm] * sin(theta * pi / 180);
                z0t1_mm = z01[mm];
                xt1 = x1 * cos(theta * pi / 180) + y1 * sin(theta * pi / 180);
                yt1 = y1 * cos(theta * pi / 180) - x1 * sin(theta * pi / 180);
                zt1 = z1;
                rt1_mm = sqrt(pow(yt1 - y0t1_mm, 2) + pow(zt1 - z0t1_mm, 2)) + epsilon;
                if (fabsf(xt1 - x0t1_mm) <= oft * dx && rt1_mm<=radius)
                    {
                        v=C_VOLUME(c,t);
                        sum_volume +=v;
                        vel=C_U(c,t);
                        average_v=v*vel;
                        sum_v +=average_v;
                    }
                    
            end_c_loop_int(c,t)     
            }
        sum_volume=PRF_GRSUM1(sum_volume);
        sum_v=PRF_GRSUM1(sum_v); 
        inflow=sum_v/sum_volume;
        Message0("inflow=%g,sum_volume=%g,sum_v=%g\n",inflow,sum_volume,sum_v);        
    }
    
#endif  
}


DEFINE_SOURCE(x_source, c, t, dS, eqn)
{
#if !RP_HOST

    real rho, loc[ND_ND], x, y, z, r, locAzimuth, epsilon = 1e-10;
    real Ux, Uy, Uz, Uret, omega;
    real xt, yt, zt, rt, x0t, y0t, z0t, locAzimuth_nt;
    real chord_j, beta_j, Cl_j, Cd_j, alpha_j, phi, Un, Ut, Uref, sin_phi, cos_phi;
    real a0 = 0.0, a1 = 1.0, aprime0 = 0.0, aprime1 = 1.0, Ftip = 1.0, Fhub = 1.0;
    real xc, yc, zc, d, eta = 1.0, eps = 0.03;
    real df_n, df_t, sin_f, cos_f;
    real ncon, tcon, source;
    real volume;
    int i, j, k, wt_i;    /* Cyclic variable */

    /* Chord and twist data from Wu & Porté-Agel, 2011 */
    real blade[N_BE] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.075};       /* Blade points (m) */
    real chord[N_BE] = {0.0139, 0.0147, 0.0147, 0.0141, 0.0131, 0.0118, 0.0102, 0.0060};  /* Chord length (m) */
    real twist[N_BE] = {20.50, 20.90, 19.83, 16.91, 13.19, 10.67, 9.12, 6.66};  /* Twist angle (degree) */
    
    /* Chord and twist data from XXXXXXXXXXXXX, 20XX */
    // real blade[N_BE] = {5.9057568, 9.043176, 10.5196575, 14.2107555, 17.3481465, 21.4083825, 24.545844, 28.4214405, 31.189764, 35.065431, 38.0183235, 42.447627, 45.5850885, 49.0916175, 52.5981465, 56.4738135, 59.611275, 63.8560095, 66.993471, 70.5};       /* Blade points (m) */
    // real chord[N_BE] = {2.57128, 3.23534, 3.38278, 3.59137, 3.54159, 3.35626, 3.12187, 2.83812, 2.60379, 2.30772, 2.08567, 1.80182, 1.66589, 1.49297, 1.38159, 1.23322, 1.10959, 0.887314, 0.579074, 0};  /* Chord length (m) */
    // real twist[N_BE] = {18.6945, 16.3388, 15.161, 12.99055, 11.1921, 9.26935, 7.7802, 6.1672, 5.11195, 3.932085, 3.124365, 2.12992, 1.508125, 1.071445, 0.696845, 0.383839, 0.069052, 0.0562465, 0.065   ,   0};  /* Twist angle (degree) */

    /* Lift and drag data of NACA0012 (defaul attack angle unit is degree) */
    real attack[N_LD] = {-180, -177.79, -172.18, -168.33, -165.1, -159.59, -154.91, -150.61, -145.49, -141.27, -134.54, -127.58, -124.17, -120.56, -115.24, -109.23, -104.33, -100.65, -96.72, -88.29, -82.73, -78.79, -74.68, -68.72, -62.9, -58.21, -53.67, -50.77, -44.08, -38.09, -32.99, -29.64, -24.44, -19, -12.78, -9.56, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 21.61, 26.07, 30.47, 36.24, 43.09, 47.56, 55.54, 60.67, 66.93, 71.07, 75.5, 80.75, 86.9, 91.87, 95.79, 99.84, 104.93, 110.44, 115.85, 124.05, 129.73, 136.34, 140.53, 144.46, 145.07, 150.54, 155.32, 160.62, 163.39, 168.46, 175.01, 180};
    real Cl[N_LD] = {-0.088, 0.083, 0.701, 1.081, 0.661, 0.477, 0.612, 0.697, 0.757, 0.716, 0.643, 0.572, 0.537, 0.466, 0.361, 0.287, 0.179, 0.094, 0.022, -0.099, -0.192, -0.261, -0.332, -0.437, -0.538, -0.597, -0.656, -0.673, -0.68, -0.698, -0.655, -0.535, -0.386, -0.354, -0.463, -0.6, -0.504, -0.382, -0.261, -0.139, -0.017, 0.104, 0.226, 0.347, 0.449, 0.574, 0.694, 0.825, 0.962, 1.073, 1.142, 1.211, 1.288, 1.355, 1.418, 1.396, 1.36, 1.319, 1.301, 1.25, 1.209, 0.806, 0.825, 0.96, 0.955, 0.926, 0.865, 0.697, 0.626, 0.495, 0.401, 0.303, 0.199, 0.099, -0.014, -0.112, -0.214, -0.317, -0.427, -0.539, -0.629, -0.695, -0.796, -0.77, -0.729, -0.721, -0.674, -0.597, -0.843, -1.324, -0.88, -0.475, -0.088};
    real Cd[N_LD] = {0.036, 0.039, 0.07, 0.18, 0.268, 0.344, 0.426, 0.536, 0.681, 0.772, 0.829, 0.886, 0.938, 1, 1.069, 1.108, 1.145, 1.167, 1.183, 1.18, 1.166, 1.148, 1.112, 1.044, 1.003, 0.965, 0.902, 0.825, 0.702, 0.602, 0.508, 0.429, 0.28, 0.191, 0.012, 0.011, 0.011, 0.01, 0.01, 0.009, 0.009, 0.008, 0.008, 0.008, 0.007, 0.007, 0.008, 0.008, 0.009, 0.011, 0.012, 0.013, 0.015, 0.019, 0.024, 0.03, 0.039, 0.049, 0.062, 0.077, 0.094, 0.288, 0.444, 0.578, 0.721, 0.811, 0.874, 1.022, 1.095, 1.144, 1.179, 1.235, 1.268, 1.272, 1.265, 1.26, 1.254, 1.205, 1.127, 1.094, 0.978, 0.872, 0.805, 0.732, 0.633, 0.616, 0.466, 0.308, 0.271, 0.195, 0.073, 0.028, 0.036};
    
    /* Lift and drag data of NACAXXXX (defaul attack angle unit is degree) */
    // real attack[N_LD] = {-4.98918,-3.92089,-2.54747,-1.3265,0.0464392,1.19119,2.48814,3.78557,5.0063,6.30185,7.59787,8.89318,10.0353,11.253,12.394,13.61,15.0526,16.4949,17.5566,18.7701,19.9835};
    // real Cl[N_LD] = {-0.215173,-0.0977905,0.0504762,0.186397,0.32228,0.452027,0.587929,0.736214,0.865943,0.964693,1.07583,1.16839,1.23002,1.27926,1.30993,1.31582,1.28451,1.247,1.19101,1.12879,1.06657};
    // real Cd[N_LD] = {0.0173875,0.01671,0.016646,0.0165895,0.0171429,0.0170826,0.0170186,0.0169546,0.0175118,0.0186975,0.0204931,0.0229136,0.0253491,0.0290081,0.0339206,0.0394504,0.0499341,0.095208,0.139872,0.164035,0.190066};

    for (i = 0; i < N_LD; i++)
    {
        attack[i] = attack[i] * pi / 180;
    }    /* Degree to rad for calculation */

    /* Obtain density of the fluid */
    rho = C_R(c, t);

    /* Obtain volume of the current cell */
    volume = C_VOLUME(c, t); 

    for (wt_i = 0; wt_i < N_WT; wt_i++)
    {
        /* Parameters selection by current simulated yawd turbine Index */
        theta_i = theta_wt[wt_i];
        radius_i = radius_wt[wt_i];
        pitch_i = pitch_wt[wt_i];
        inflow_i = inflow_wt[wt_i];
        tsr_i = tsr_wt[wt_i];
        rpm_i = rpm_wt[wt_i];
        x_i = x_wt[wt_i];
        y_i = y_wt[wt_i];
        z_i = z_wt[wt_i];

        real beta[N_BE];          /* Beta = twist angle + pitch angle */
        for (i = 0; i < N_BE; i++)
        {
            beta[i] = (twist[i] + pitch_i) * pi / 180;
        }    /* Degree to rad for calculation */

        /* Obtain the coordinates of each cell centroid  and compute the relative values */
        C_CENTROID(loc, c, t);
        x = loc[0];
        y = loc[1];
        z = loc[2];
        r = sqrt(pow(y - y_i, 2) + pow(z - z_i, 2)) + epsilon;
        // locAzimuth = atan2(z - z0, y - y0);

        /* Local velocity at the rotor disc */
        Ux = C_U(c, t);   
        Uy = C_V(c, t);
        Uz = C_W(c, t);
        Uret = sqrt(pow(Ux, 2) + pow(Uy, 2) + pow(Uz, 2));
        // omega = rpm_i * 2 * pi / 60 * cos(theta_i * pi / 180)
        omega = inflow_i * cos(theta_i * pi / 180) * tsr_i / radius_i;

        /* Yawed turbine rotor coordinate transformation */
        xt_i = x_i * cos(theta_i * pi / 180) + y_i * sin(theta_i * pi / 180);
        yt_i = y_i * cos(theta_i * pi / 180) - x_i * sin(theta_i * pi / 180);
        zt_i = z_i;
        xt = x * cos(theta_i * pi / 180) + y * sin(theta_i * pi / 180);
        yt = y * cos(theta_i * pi / 180) - x * sin(theta_i * pi / 180);
        zt = z;
        rt = sqrt(pow(yt - yt_i, 2) + pow(zt - zt_i, 2)) + epsilon;
        // locAzimuth_nt = atan2(zt - zt_i, yt - yt_i);

        /* Yawed turbine tower center coordinate transformation */
        x_tow_i = (x_i + x_offset_tow) * cos(theta_i * pi / 180) + (y_i + y_offset_tow) * sin(theta_i * pi / 180);
        y_tow_i = (y_i + y_offset_tow) * cos(theta_i * pi / 180) - (x_i + x_offset_tow) * sin(theta_i * pi / 180);
        z_tow_i = (z_i + z_offset_tow);

        /* Calculate rotor force and store it in UDM */
        if (fabs(xt - xt_i) <= oft * dx)
        {
            /* Determine whether it is within rotor area */
            if (rt <= blade[N_BE - 1])
            {
                /* Interpolation for finding local chord and beta at specific radical distance r */
                if (rt >= blade[0])
                {
                    for (i = 0; i < N_BE; i++)
                    {
                        if (rt >= blade[i] && rt <= blade[i + 1])
                        {
                            chord_j = (chord[i + 1] - chord[i]) * (rt - blade[i]) / (blade[i + 1] - blade[i]) + chord[i];
                            beta_j = (beta[i + 1] - beta[i]) * (rt - blade[i]) / (blade[i + 1] - blade[i]) + beta[i];
                        }
                    }
                }
                else {chord_j = chord[0]; beta_j = beta[0];}

                /* ----------------------  Solving BEM for each airfoil blade element  ----------------------- */
                /* local normal velocity at the rotor disc */
                Un = Ux * cos(theta_i * pi / 180) + Uy * sin(theta_i * pi / 180) + epsilon;
                for (j = 0; j <= 100; j++)
                {   
                    // Un = inflow_i * cos(theta_i * pi / 180) * (1 - a0);
                    Ut = omega * rt * (1 + aprime0) - inflow_i * sin(theta_i * pi / 180) * ((zt - zt_i) / rt);
                    // Ut = omega * rt * (1 + aprime0);
                    phi = atan(Un / Ut);                         /* Inflow angle */
                    Uref = sqrt(pow(Un, 2) + pow(Ut, 2));        /* Relative wind velocity */
                    sin_phi = Un / Uref;                         /* Sin of inflow angle */
                    cos_phi = Ut / Uref;                         /* Cos of inflow angle */

                    /* Interpolation for finding local lift and drag coefficient at specific radical distance r */
                    if (rt <= blade[0])
                    {
                        Cl_j = 0;      /*  Default Lift cofficient near rotor center */
                        Cd_j = 1;      /*  Default Drag cofficient near rotor center */
                    }
                    else
                    {
                        alpha_j = phi - beta_j;                                    /* Attack angle */
                        for (i = 0; i < N_LD; i++)
                        {
                            if (alpha_j >= attack[i] && alpha_j <= attack[i + 1])
                            {
                                Cl_j = (Cl[i + 1] - Cl[i]) * (alpha_j - attack[i]) / (attack[i + 1] - attack[i]) + Cl[i];
                                Cd_j = (Cd[i + 1] - Cd[i]) * (alpha_j - attack[i]) / (attack[i + 1] - attack[i]) + Cd[i];
                            }
                        }
                    }
                    // Ftip = (2 / pi) * acos(exp(- (B * (radius_i - rt)) / (2 * (rt + 1e-6) * sin_phi)));  /* Blade tip loss */
                    // Fhub = (2 / pi) * acos(exp(- (B * (rt - blade[0])) / (2 * rt * sin_phi)));         /* Hub loss */
                    
                    aprime1 = pow(((8 * pi* rt * Ftip * Fhub * sin_phi * cos_phi / (chord_j * B * (Cl_j * sin_phi - Cd_j * cos_phi))) - 1), -1);
                    if (fabs(aprime1 - aprime0) <= 1e-6) {break;} else {a0 = a1;}          /* Convergence condition */

                    // a1 = pow(((8 * pi* rt * Ftip * Fhub * sin_phi * sin_phi / (chord_j * B * (Cl_j * cos_phi + Cd_j * sin_phi))) + 1), -1);
                    // aprime1 = pow(((8 * pi* rt * Ftip * Fhub * sin_phi * cos_phi / (chord_j * B * (Cl_j * sin_phi - Cd_j * cos_phi))) - 1), -1);
                    // if (fabs(a1 - a0) <= 1e-6 && fabs(aprime1 - aprime0) <= 1e-6) {break;}
                    // else {a0 = a1; aprime0 = aprime1;}   /* Convergence condition */
                }
            }
            else
            {
                Uref = 0; chord_j = 0; sin_phi = 0; cos_phi = 0; Cl_j = 0; Cd_j = 0;
            }

            /* Gussian kernel convolution */
            // xc = - (y - y_i) * sin(theta_i * pi / 180);
            // yc = (y - y_i) * cos(theta_i * pi / 180);
            // d = sqrt(pow(((x - x_i) - xc), 2) + pow(((y - y_i) - yc), 2));
            // d = fabs(xt - xt_i);
            // eta = (1 / (eps * pow(pi, 1 / 2))) * exp(- pow(d / eps, 2));
            
            df_n = 0.5 * eta * rho * pow(Uref, 2) * chord_j * (Cl_j * cos_phi + Cd_j * sin_phi);  /* Thrust force on each blade element */
            df_t = 0.5 * eta * rho * pow(Uref, 2) * chord_j * (Cl_j * sin_phi - Cd_j * cos_phi);  /* Torque force on each blade element */
            sin_f = (zt - zt_i) / rt;            /* Sin value of cell Azimuth angle */
            cos_f = (yt - yt_i) / rt;            /* Cos value of cell Azimuth angle */

            C_UDMI(c, t, f_n) = - df_n * B / (2 * pi * rt * dx);      /* Thrust force in each cell */
            C_UDMI(c, t, f_t) = - df_t * B / (2 * pi * rt * dx);      /* Torque force in each cell */
            C_UDMI(c, t, total_f_n) = C_UDMI(c, t, f_n) * volume;
            C_UDMI(c, t, total_f_t) = C_UDMI(c, t, f_t) * volume * rt;
            C_UDMI(c, t, f_x) = C_UDMI(c, t, f_n) * cos(theta_i * pi / 180) + sin_f * C_UDMI(c, t, f_t) * sin(theta_i * pi / 180);   /* x force in each cell */
            C_UDMI(c, t, f_y) = C_UDMI(c, t, f_n) * sin(theta_i * pi / 180) - sin_f * C_UDMI(c, t, f_t) * cos(theta_i * pi / 180);   /* y force in each cell */
            C_UDMI(c, t, f_z) = C_UDMI(c, t, f_t) * cos_f;           /* z force in each cell */
        }

        /* Calculate nacelle force stored in UDM */
        if (((xt - xt_i) > oft * dx) && ((xt - xt_i) <= (oft * dx + length_nac)) && (rt <= (diameter_nac / 2)))
        {
            ncon = Cd_nac * 0.5 * rho;
            C_UDMI(c, t, f_nx) = -ncon * Uret * Ux;
            C_UDMI(c, t, f_ny) = -ncon * Uret * Uy;
            C_UDMI(c, t, f_nz) = -ncon * Uret * Uz;
        }

        /* Calculate tower force stored in UDM */
        if ((z >= 0.0 && z < z_tow_i) && (sqrt(pow(x - x_tow_i, 2) + pow(y - y_tow_i, 2)) <= (diameter_tow / 2)))
        {
            tcon = Cd_tow * 0.5 * rho;
            C_UDMI(c, t, f_tx) = -tcon * Uret * Ux;
            C_UDMI(c, t, f_ty) = -tcon * Uret * Uy;
            C_UDMI(c, t, f_tz) = -tcon * Uret * Uz;
        }

        source = C_UDMI(c, t, f_x) + C_UDMI(c, t, f_nx) + C_UDMI(c, t, f_tx);
        dS[eqn] = 0;
    }

    return source;

#endif  /*!RP_HOST*/
}


DEFINE_SOURCE(y_source, c, t, dS, eqn)
{
#if !RP_HOST

    real source;    
    source = C_UDMI(c, t, f_y) + C_UDMI(c, t, f_ny) + C_UDMI(c, t, f_ty);
    dS[eqn] = 0;
    
    return source;
    
#endif  /*!RP_HOST*/
}


DEFINE_SOURCE(z_source, c, t, dS, eqn)
{
#if !RP_HOST

    real source;
    source = C_UDMI(c, t, f_z) + C_UDMI(c, t, f_nz) + C_UDMI(c, t, f_tz);
    dS[eqn] = 0;

    return source;
    
#endif  /*!RP_HOST*/
}


/* calculate generated power for each turbine */
DEFINE_ON_DEMAND(on_demand_calc_total_power)
{
#if RP_NODE

    cell_t c;
    Thread *t;
    Domain *domain;
    domain = Get_Domain(1);
    real f_t_sum = 0;
    real f_n_sum = 0;
    real power;
    
    thread_loop_c(t, domain)
    {
        if (FLUID_THREAD_P(t))
        {
            begin_c_loop_int(c, t)
            {                                 
                f_t_sum += C_UDMI(c, t, total_f_t);
                f_n_sum += C_UDMI(c, t, total_f_n);              
            }
            end_c_loop_int(c, t)            
        }
    }
                
    power = f_t_sum * 2 * pi * rpm / 60 / 1000000;
    f_n_sum = PRF_GRSUM1(f_n_sum);
    f_t_sum = PRF_GRSUM1(f_t_sum);
    power = PRF_GRSUM1(power);

    Message0("f_n_sum=%g, f_t_sum=%g, p=%g\n",f_n_sum, f_t_sum, power);

#endif  
}
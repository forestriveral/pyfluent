/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_PROFILE(inlet_velocity, thread, index);
extern DEFINE_PROFILE(k_profile, thread, index);
extern DEFINE_PROFILE(w_profile, thread, index);
extern DEFINE_SOURCE(k_source, c, t, dS, eqn);
extern DEFINE_SOURCE(w_source, c, t, dS, eqn);
extern DEFINE_EXECUTE_ON_LOADING(on_loading, libname);
__declspec(dllexport) UDF_Data udf_data[] = {
{"inlet_velocity", (void(*)())inlet_velocity, UDF_TYPE_PROFILE},
{"k_profile", (void(*)())k_profile, UDF_TYPE_PROFILE},
{"w_profile", (void(*)())w_profile, UDF_TYPE_PROFILE},
{"k_source", (void(*)())k_source, UDF_TYPE_SOURCE},
{"w_source", (void(*)())w_source, UDF_TYPE_SOURCE},
{"on_loading", (void(*)())on_loading, UDF_TYPE_EXECUTE_ON_LOADING},
};
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}

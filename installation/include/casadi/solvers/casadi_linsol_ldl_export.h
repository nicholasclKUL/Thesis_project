
#ifndef CASADI_LINSOL_LDL_EXPORT_H
#define CASADI_LINSOL_LDL_EXPORT_H

#ifdef CASADI_LINSOL_LDL_STATIC_DEFINE
#  define CASADI_LINSOL_LDL_EXPORT
#  define CASADI_LINSOL_LDL_NO_EXPORT
#else
#  ifndef CASADI_LINSOL_LDL_EXPORT
#    ifdef casadi_linsol_ldl_EXPORTS
        /* We are building this library */
#      define CASADI_LINSOL_LDL_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define CASADI_LINSOL_LDL_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef CASADI_LINSOL_LDL_NO_EXPORT
#    define CASADI_LINSOL_LDL_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef CASADI_LINSOL_LDL_DEPRECATED
#  define CASADI_LINSOL_LDL_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CASADI_LINSOL_LDL_DEPRECATED_EXPORT
#  define CASADI_LINSOL_LDL_DEPRECATED_EXPORT CASADI_LINSOL_LDL_EXPORT CASADI_LINSOL_LDL_DEPRECATED
#endif

#ifndef CASADI_LINSOL_LDL_DEPRECATED_NO_EXPORT
#  define CASADI_LINSOL_LDL_DEPRECATED_NO_EXPORT CASADI_LINSOL_LDL_NO_EXPORT CASADI_LINSOL_LDL_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CASADI_LINSOL_LDL_NO_DEPRECATED
#    define CASADI_LINSOL_LDL_NO_DEPRECATED
#  endif
#endif

#endif /* CASADI_LINSOL_LDL_EXPORT_H */

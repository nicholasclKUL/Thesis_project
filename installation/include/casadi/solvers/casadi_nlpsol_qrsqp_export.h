
#ifndef CASADI_NLPSOL_QRSQP_EXPORT_H
#define CASADI_NLPSOL_QRSQP_EXPORT_H

#ifdef CASADI_NLPSOL_QRSQP_STATIC_DEFINE
#  define CASADI_NLPSOL_QRSQP_EXPORT
#  define CASADI_NLPSOL_QRSQP_NO_EXPORT
#else
#  ifndef CASADI_NLPSOL_QRSQP_EXPORT
#    ifdef casadi_nlpsol_qrsqp_EXPORTS
        /* We are building this library */
#      define CASADI_NLPSOL_QRSQP_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define CASADI_NLPSOL_QRSQP_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef CASADI_NLPSOL_QRSQP_NO_EXPORT
#    define CASADI_NLPSOL_QRSQP_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef CASADI_NLPSOL_QRSQP_DEPRECATED
#  define CASADI_NLPSOL_QRSQP_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CASADI_NLPSOL_QRSQP_DEPRECATED_EXPORT
#  define CASADI_NLPSOL_QRSQP_DEPRECATED_EXPORT CASADI_NLPSOL_QRSQP_EXPORT CASADI_NLPSOL_QRSQP_DEPRECATED
#endif

#ifndef CASADI_NLPSOL_QRSQP_DEPRECATED_NO_EXPORT
#  define CASADI_NLPSOL_QRSQP_DEPRECATED_NO_EXPORT CASADI_NLPSOL_QRSQP_NO_EXPORT CASADI_NLPSOL_QRSQP_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CASADI_NLPSOL_QRSQP_NO_DEPRECATED
#    define CASADI_NLPSOL_QRSQP_NO_DEPRECATED
#  endif
#endif

#endif /* CASADI_NLPSOL_QRSQP_EXPORT_H */

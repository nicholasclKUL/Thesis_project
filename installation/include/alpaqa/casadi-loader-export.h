
#ifndef CASADI_LOADER_EXPORT_H
#define CASADI_LOADER_EXPORT_H

#ifdef CASADI_LOADER_STATIC_DEFINE
#  define CASADI_LOADER_EXPORT
#  define CASADI_LOADER_NO_EXPORT
#else
#  ifndef CASADI_LOADER_EXPORT
#    ifdef casadi_loader_EXPORTS
        /* We are building this library */
#      define CASADI_LOADER_EXPORT 
#    else
        /* We are using this library */
#      define CASADI_LOADER_EXPORT 
#    endif
#  endif

#  ifndef CASADI_LOADER_NO_EXPORT
#    define CASADI_LOADER_NO_EXPORT 
#  endif
#endif

#ifndef CASADI_LOADER_DEPRECATED
#  define CASADI_LOADER_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CASADI_LOADER_DEPRECATED_EXPORT
#  define CASADI_LOADER_DEPRECATED_EXPORT CASADI_LOADER_EXPORT CASADI_LOADER_DEPRECATED
#endif

#ifndef CASADI_LOADER_DEPRECATED_NO_EXPORT
#  define CASADI_LOADER_DEPRECATED_NO_EXPORT CASADI_LOADER_NO_EXPORT CASADI_LOADER_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CASADI_LOADER_NO_DEPRECATED
#    define CASADI_LOADER_NO_DEPRECATED
#  endif
#endif

#endif /* CASADI_LOADER_EXPORT_H */

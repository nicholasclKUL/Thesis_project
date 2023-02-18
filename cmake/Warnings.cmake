function(add_warnings_target tgt_name warnings_as_errors with_analyzer)

    # GCC, Clang, AppleClang
    set(COMMON_WARNINGS
        -fdiagnostics-show-option
        -Wall
        -Wextra
        -pedantic
        -Wpedantic
        -pedantic-errors
        -Wdouble-promotion
        -Wswitch-default
        -Wswitch-enum
        -Wimplicit-fallthrough
        -Wuninitialized
        -Wno-missing-braces
        -fopenmp
    )
    # GCC
    set(GCC_WARNINGS
        -Wno-error=unused-but-set-variable
        $<$<COMPILE_LANGUAGE:CXX>:-Wsuggest-override>
        -Wno-error=attributes
        -Wno-psabi
    )
    # Clang, AppleClang
    set(CLANG_WARNINGS
        -Wno-error=unknown-warning-option
        -Wno-newline-eof
        -Wno-error=unused-but-set-variable
        $<$<COMPILE_LANGUAGE:CXX>:-Winconsistent-missing-override>
        -Wno-gnu-zero-variadic-macro-arguments
    )
    # MSVC (Microsoft)
    set(MSVC_WARNINGS
        /W4
        /wd4127 # conditional expression is constant
        /wd4458 # declaration of 'x' hides class member
        /permissive-
    )
    # Intel ICC
    set(INTEL_WARNINGS 
        -Wall
        -Wextra
    )

    # Enable warnings as errors
    if (warnings_as_errors)
        if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            list(APPEND MSVC_WARNINGS /WX)
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            list(APPEND INTEL_WARNINGS -Werror)
        else()
            list(APPEND COMMON_WARNINGS -Werror)
        endif()
    endif()

    # Static analyzer
    if (with_analyzer)
        list(APPEND COMMON_WARNINGS -fanalyzer)
    endif()

    # Add target that defines all the warning options in its interface.
    add_library(${tgt_name} INTERFACE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${tgt_name} INTERFACE
            ${COMMON_WARNINGS} ${GCC_WARNINGS})
    elseif (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
        target_compile_options(${tgt_name} INTERFACE
            ${COMMON_WARNINGS} ${CLANG_WARNINGS})
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        target_compile_options(${tgt_name} INTERFACE
            ${MSVC_WARNINGS})
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        target_compile_options(${tgt_name} INTERFACE
            ${INTEL_WARNINGS})
    else()
        message(WARNING "No known warnings for this compiler")
    endif()

endfunction()
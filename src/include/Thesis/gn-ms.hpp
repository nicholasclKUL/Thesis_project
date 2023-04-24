#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>

#include <functional>

namespace alpaqa {

void indexing (crvec x, crvec u, Box &X, Box &U) {
    
    index_t nx = x.size(),
            nu = u.size(), 
            nx_J = 0, 
            nx_K = 0, 
            nu_J = 0, 
            nu_K = 0;
    
    vec indexes_J, indexes_K;

    for (size_t i = 0; i < nx; ++i) {
        if ((x(i) == X.lowerbound) || x(i) == (X.upperbound)) {
            indexes_J.push_back(i);
            ++nx_J;
        } else {
            indexes_K.push_back(i);
            ++nx_K;
        }        
    }
    
    for (size_t i = 0; i < nu; ++i) {
        if ((u(i) == U.lowerbound) || u(i) == (U.upperbound)) {
            indexes_U.push_back(i);
            ++nu_J;
        } else {
            indexes_K.push_back(i);
            ++nu_K;
        }        
    }

}
    
}
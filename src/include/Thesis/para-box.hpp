#pragma once 

#include <alpaqa/problem/box.hpp>

namespace alpaqa{

template <Config Conf>
inline auto para_project(const auto &v,       ///< [in] The vector to project
                        const Box<Conf> &box, ///< [in] The box to project onto
                        unsigned int i
) {
    USING_ALPAQA_CONFIG(Conf);
    return std::min(std::max(v[i], box.upperbound[i]-v[i]), box.lowerbound[i]-v[i]);     
}

/// Get the difference between the given vector and its projection.
/// @f[ v - \Pi_C(v) @f]
/// @warning    Beware catastrophic cancellation!
template <Config Conf>
inline auto para_projecting_difference(const auto &v,        ///< [in] The vector to project
                                       const Box<Conf> &box, ///< [in] The box to project onto
                                       unsigned int i
) {
    return v[i] - para_project(v, box, i);
}

}
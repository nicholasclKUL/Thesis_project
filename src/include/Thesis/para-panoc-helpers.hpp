#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <Kokkos_Core.hpp>

namespace alpaqa{
    
using Conf      = DefaultConfig; 
using Problem   = alpaqa::TypeErasedControlProblem<config_t>;

USING_ALPAQA_CONFIG(Conf);

void team_estimation(int &nteams,    // number of workers
                     int &nthrds     // number of threads per worker
                     ){

};

} // alpaqa namespace


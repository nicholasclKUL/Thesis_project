#include <Thesis/Thesis.hpp>
#include <taskflow/taskflow.hpp> 
#include <iostream>

namespace Thesis {
    void hello(){
     std::cout <<"hello!" << std::endl;   
    }

    void test(){
        tf::Executor executor;
        tf::Taskflow taskflow;

        auto [A, B, C, D] = taskflow.emplace(  // create four tasks
        [] () { std::cout << "TaskA\n"; },
        [] () { std::cout << "TaskB\n"; },
        [] () { std::cout << "TaskC\n"; },
        [] () { std::cout << "TaskD\n"; } 
        );                                  
                                      
        A.precede(B, C);  // A runs before B and C
        D.succeed(B, C);  // D runs after  B and C
                                      
        executor.run(taskflow).wait(); 
    }

}
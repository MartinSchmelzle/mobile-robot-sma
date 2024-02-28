#include <iostream>
#include <limits>
#include "../custom_datatypes.h"
#include "../slimemould_func_prototypes.h"
#include <chrono>

void printVector(const std::vector<double>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << " ";
    }
}

void print2DVector(const std::vector<std::vector<double>>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec[i].size(); ++j) {
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


void printSMAsol(const SMAsol_struct& sol) {
    std::cout << "Violation: " << sol.Violation << std::endl;
    std::cout << "L: " << sol.L << std::endl;
    std::cout << "L2: " << sol.L2 << std::endl;
    std::cout << "XS: ";
    for (int i = 0; i < 10; ++i) {
        std::cout << sol.XS[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "YS: ";
    for (int i = 0; i < 10; ++i) {
        std::cout << sol.YS[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "internalpath:" << std::endl;
    std::cout << "  startTurn: " << sol.internalpath.startTurn << std::endl;
    std::cout << "  endTurn: " << sol.internalpath.endTurn << std::endl;
    std::cout << "  s: ";
    printVector(sol.internalpath.s);
    std::cout << std::endl;
    std::cout << "  a: ";
    printVector(sol.internalpath.a);
    std::cout << std::endl;
    std::cout << "  r: ";
    printVector(sol.internalpath.r);
    std::cout << std::endl;
    std::cout << "  m:" << std::endl;
    print2DVector(sol.internalpath.m);
    std::cout << "  dis: ";
    printVector(sol.internalpath.dis);
    std::cout << std::endl;
}

bool mr_sma_alg(model_struct *model, path_struct *path)
{

    //read out kinematics- and AGV- related data from config file
    kinematics_time* kin = new kinematics_time;
    alg_data* alg = new alg_data;
    kin-> a_acc = findinh("a_acceleration")[0];//todo: I'm going to just use the values for straight driving here. Insert values for circle
    kin->v_max = findinh("v_max")[0];
    kin->v_desired = kin->v_max*0.95;
    kin->a_break = -findinh("a_break")[0];
    kin->P_idle = findinh("P_idle")[0]; //in W
    kin->mass = findinh("mass_empty")[0]; //weight with no payload
    alg->max_iter=findinh("max_iter")[0];
    alg->distance_discrete=findinh("distance_discrete")[0];

    //  run SMA Algorithm
    //if the first computed solution is infeasible, we give it four more attempts
    double inf = std::numeric_limits<double>::infinity();
    output_SMA out;
    SMAsol_struct SMAsol;
    unsigned attempt = 1;
    double bestPositions[6];
    //printSMAsol(SMAsol);
    out.Vio=inf;
    while(attempt<=5 && out.Vio>=0.2) {
        out=SMA(model, &SMAsol, kin, alg);
        for (int i{0}; i < 6; i++) {bestPositions[i] = out.waypoints[i];}
        F00(bestPositions,model,&SMAsol,kin,1);//Call F00 once more to populate SMAsol struct

        if(attempt==5 &&out.Vio>=0.2)
        {
            std::cout<<"Path calculation unsuccessful after five attempts.\n";
            return -1;
        }
        attempt+=1;
    }

    // Convert path into different format for path processing
    //todo: change path processing to work with SMAsol instead of requiring reformatting and path struct
    path_struct p=format_path(&SMAsol,model);
    //std::cout<<"path formatted\n";

    //Assign value of p to path (pointer) with appropriate error handling
    if (path != nullptr) {
            *path = p; // Copy the contents of p into the object pointed to by path
        } else {
            std::cout << "Path Output: Invalid pointer to path_struct!\n";
        }

    std::cout<<"SMA done!";
    delete kin;
    delete alg;
    return 0;
}
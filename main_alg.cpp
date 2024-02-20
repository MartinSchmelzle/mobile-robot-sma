#include <iostream>
#include <limits>
#include "custom_datatypes.h"
//#include "custom_math_operations.cpp"
//#include "F00.cpp"
#include "SMA.cpp"
#include "format_path.cpp"
#include "getData.cpp"
#include "slimemould_func_prototypes.h"
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

bool main_alg(model_struct *model, path_struct *path, const unsigned max_iter)
{
    SMAsol_struct SMAsol;
    //printSMAsol(SMAsol);

    double bestPositions[6];

//read out kinematics-related data from h
kinematics_time* kin = new kinematics_time;
kin-> a_acc = findinh("s.Acc.a")[0];//todo: I'm going to just use the values for straight driving here. Insert values for circle
kin->v_max = findinh("s.Con.v_Max")[0];
kin->v_desired = kin->v_max*0.95;
kin->a_break = -findinh("s.Slw.a")[0];

//  run SMA Algorithm
//if the first computed solution is infeasible, we give it two more attempts
double inf = std::numeric_limits<double>::infinity();
output_SMA out;
out.Vio=inf;
unsigned attempt = 1;
while(attempt<=5 && out.Vio>=0.2)
{
  out=SMA(model, &SMAsol, kin, max_iter);
  for (int i{0}; i < 6; i++) {bestPositions[i] = out.waypoints[i];}
  F00(bestPositions,model,&SMAsol,kin,1);//Call F00 once more to populate SMAsol struct
  //printSMAsol(SMAsol);
  //std::cout<<"SMA output: ";
  //for(int w=0; w<6;w++){std::cout<<bestPositions[w]<<",";} 
  //std::cout<<"Violation: "<<out.Vio<<", Fitness: "<<out.res<<"\n";

  //std::cout<<"attempt "<<attempt<<". Best Violation: "<<out.Vio<<"\n";
  if(attempt==5 &&out.Vio>=0.2)
  {
    std::cout<<"Path calculation unsuccessful after five attempts.\n";
    return -1;
  }
  attempt+=1;
}
//std::cout<<"attempt successfull\n";
std::array<std::array<double,2>,4> centerPoints;
centerPoints[0]={SMAsol.internalpath.m[0][0],SMAsol.internalpath.m[0][1]};
centerPoints[1]={SMAsol.internalpath.m[1][0],SMAsol.internalpath.m[1][1]};
centerPoints[2]={SMAsol.internalpath.m[2][0],SMAsol.internalpath.m[2][1]};
centerPoints[3]={SMAsol.internalpath.m[3][0],SMAsol.internalpath.m[3][1]};

  // Convert path into different format for processing outside the code
  path_struct p=format_path(&SMAsol,model,centerPoints);
  //std::cout<<"path formatted\n";

  //Assign value of p to path (pointer)
  if (path != nullptr) {
        *path = p; // Copy the contents of p into the object pointed to by path
    } else {
        std::cout << "Path Output: Invalid pointer to path_struct!\n";
    }

  std::cout<<"SMA done!";
  delete kin;
  return 0;
}
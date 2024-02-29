//#include "get_Energy_v2.cpp"
//#include <iostream>
#include "../custom_datatypes.h"
#include "../slimemould_func_prototypes.h"
#include <cstdint>
#include <iostream>
#include <cmath>

//  Based on an algorithm developed by Yarpiz (Project Code: YPAP115)
//  Project Title: Path Planning using PSO in MATLAB
//  Publisher: Yarpiz (www.yarpiz.com)

//set accuracy to 0 for rough calculation, to 1 for fine calculation
struct_fitness F00(double x[6], const model_struct *model, SMAsol_struct *sol,kinematics_time* kin, const bool accuracy)
{
  double unchecked_XS[6],unchecked_YS[6];
  
  //This function gets called 1000s of times and we don't really care about what the historic solutions of F00 look like
  //delete everything in SMAsol_struct
  sol->xx.clear();
  sol->yy.clear();
  sol->internalpath.a.clear();
  sol->internalpath.dis.clear();
  sol->internalpath.m.clear();
  sol->internalpath.r.clear();
  sol->internalpath.s.clear();
  setArrayTo0(kin->a_array,10000);
  setArrayTo0(kin->v_array,10000);
  setArrayTo0(kin->x_array,10000);
  setArrayTo0(kin->y_array,10000);
  setArrayTo0(kin->s_array,10000);
  setArrayTo0(kin->w_array,10000);
  setArrayTo0(kin->dwdt_array,10000);
  setArrayTo0(kin->time_array,10000);

  //coordinate transfer for guide points
  //This function is the only one I adapted from MATLAB Coder, it is a bit difficult to read
  const std::array<double,8> xxyy = getPointstoXY(x, model);
  //std::cout<<"F00: getPointstoXY done\n";

  //Setup for checkPoints
  unchecked_XS[0] = model->xs;
  unchecked_YS[0] = model->ys;
  unchecked_XS[1] = xxyy[0];
  unchecked_YS[1] = xxyy[4];
  unchecked_XS[2] = xxyy[1];
  unchecked_YS[2] = xxyy[5];
  unchecked_XS[3] = xxyy[2];
  unchecked_YS[3] = xxyy[6];
  unchecked_XS[4] = xxyy[3];
  unchecked_YS[4] = xxyy[7];
  unchecked_XS[5] = model->xt;
  unchecked_YS[5] = model->yt;
  uint8_t k;
  
  //check guide points for unwanted behavior.
  checkPoints(unchecked_XS, unchecked_YS, sol, k);
  //std::cout<<"F00: checkPoints done\n";
  //printxsys(soldata_XS,soldata_YS); 

  //initialize path with radii (so we can optimize the radii later)
  initialization(sol, k);
  //std::cout<<"F00: initialization done\n";

  //optimize radii
  optimizationRad(sol, k);
  //std::cout<<"F00: optimizationRad done\n";

  //in earlier iterations, we don't need our calculations to be extremely accurate
  //in rough mode, code runs faster, in later iterations, it produces more accurate results
  double distance_discrete;
  if(accuracy){
    distance_discrete=0.2;  
    kin->delta_t=1;
    }
  else{
    distance_discrete=1;  
    kin->delta_t=3;
    }
  
  const double totalpathlength = discretize(sol,k, model->nrEl,distance_discrete); 
  //std::cout<<"F00: discretize done\n";

  //get Violation Score
  sol->Violation = getViolationRect(model, sol, k) 
                + getViolationCircBoundary(model,sol) 
                + getViolationChargingStation(sol, model,k) 
                + getSmallRadVio(sol);

  //calculate Energy demand
  sol->L = getEnergy_v2(sol, kin);
  //std::cout<<"L: "<<sol->L<<"\n";

  struct_fitness F00_res;
  F00_res.L=sol->L;
  F00_res.Vio=sol->Violation;
  F00_res.res=sol->L +(100000 * sol->Violation)+1;
  //std::cout<<"Violation: "<<sol->Violation<<", Fitness Value: "<<F00_res.res<<"\n";

  //some error handling
  if(std::isnan(F00_res.res))
  {
    std::cout<<"F00 returned nan\n";
    std::cout<<"res: "<<F00_res.res<<", sol->L: "<<sol->L<<", sol->Violation: "<<sol->Violation<<std::endl;
    std::cout<<"L: "<<sol->L<<", totalpathlength: "<<totalpathlength<<std::endl;
  }
  else if(std::isinf(F00_res.res))
  {F00_res.res=100000000 +rand01();}

  return F00_res;
}
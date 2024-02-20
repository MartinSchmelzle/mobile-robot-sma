#include "get_Energy_v2.cpp"
#include "F00_functions.cpp"
#include <cstdint>

//  Based on an algorithm developed by Yarpiz (Project Code: YPAP115)
//  Project Title: Path Planning using PSO in MATLAB
//  Publisher: Yarpiz (www.yarpiz.com)

//set accuracy to 0 for rough calculation, to 1 for fine calculation
struct_fitness F00(double x[6], const model_struct *model, SMAsol_struct *sol,kinematics_time* kin, const bool accuracy)
{
  double unchecked_XS[6],unchecked_YS[6],g_endTurn,g_startTurn;
  int nrEl;
  std::array<double,4> x_guide,y_guide;
  
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
  const std::array<double,8> xxyy = getPointstoXY(x, model);
  //std::cout<<"F00: getPointstoXY done\n";

  //assign x_guide and y_guide correctly
  for(int i=0;i<4;i++){x_guide[i]=xxyy[i];y_guide[i]=xxyy[i+4];}
  
  //amount of points the path is divided into for collision avoidance
  nrEl = static_cast<int>(model->nrEl);

  //Setup for checkPoints
  unchecked_XS[0] = model->xs;
  unchecked_YS[0] = model->ys;
  unchecked_XS[1] = x_guide[0];
  unchecked_YS[1] = y_guide[0];
  unchecked_XS[2] = x_guide[1];
  unchecked_YS[2] = y_guide[1];
  unchecked_XS[3] = x_guide[2];
  unchecked_YS[3] = y_guide[2];
  unchecked_XS[4] = x_guide[3];
  unchecked_YS[4] = y_guide[3];
  unchecked_XS[5] = model->xt;
  unchecked_YS[5] = model->yt;
  std::array<double,10> soldata_XS,soldata_YS;
  uint8_t k;
  
  //check guide points for unwanted behavior.
  checkPoints(unchecked_XS, unchecked_YS, soldata_XS, soldata_YS, g_startTurn, g_endTurn, k);
  for(int i=0;i<10;i++)
  {sol->XS[i]=soldata_XS[i];sol->YS[i]=soldata_YS[i];}
  //std::cout<<"F00: checkPoints done\n";
  //printxsys(soldata_XS,soldata_YS); 

  //write outputs into g
  sol->internalpath.startTurn=g_startTurn;
  sol->internalpath.endTurn=g_endTurn;

  //initialize path with radii (so we can optimize the radii later)
  initialization(sol, k);
  //std::cout<<"F00: initialization done\n";

  //optimize radii
  optimizationRad(sol, k);
  //std::cout<<"F00: optimizationRad done\n";

  //in earlier iterations, we don't need our calculations to be extremely accurate
  //in rough mode, code runs faster, in later iterations, it produces more accurate results
  double distance_discrete;
  if(accuracy){distance_discrete=0.2;  kin->delta_t=1;}
  else{distance_discrete=1;  kin->delta_t=3;}
  const double totalpathlength = discretize(sol,k, model->nrEl,distance_discrete); 
  //std::cout<<"F00: discretize done\n";
  //std::cout<<"totalpathlength: "<<totalpathlength<<"\n";

  //get Violation Score
  sol->Violation = getViolationRect(model, sol, k);
  //std::cout<<"Violation through rectangles: "<<sol->Violation<<std::endl;

  const double Violation_Circ = getViolationCircBoundary(model,sol);
  //std::cout<<"Violation through circles: "<<Violation_Circ<<std::endl;
  sol->Violation+=Violation_Circ;

  //get Violation Score for the Charging Station (Penalty for entering the charging station from the wrong side, I think)
  const double ViolationCharger=getViolationChargingStation(sol, model,k);
  //std::cout<<"Violation through Charging Stations: "<<ViolationCharger<<std::endl;
  sol->Violation+=ViolationCharger; //if(model->AisChargingStation){ //  Violation += getViolationChargingStation(&(sol->internalpath), model->xs_CS,model->ys_CS, 0, model->SA, sol, k, n); //} //if(model->BisChargingStation){ //  Violation += getViolationChargingStation(&(sol->internalpath), model->xs_CS,model->ys_CS, 0, model->SA, sol, k, n);

  //the algorithm sometimes converges against paths with tiny radii => This is dangerous and needs to be prevented.
  const double ViolationSmallRad = getSmallRadVio(sol);
  sol->Violation+=ViolationSmallRad;

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
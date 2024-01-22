#include "get_Energy_v2.cpp"
#include "F00_functions.cpp"
#include <cstdint>

//  Based on an algorithm developed by Yarpiz (Project Code: YPAP115)
//  Project Title: Path Planning using PSO in MATLAB
//  Publisher: Yarpiz (www.yarpiz.com)
using xycirc = std::vector<std::vector<double>>;
double F00(double x[6], const model_struct *model, SMAsol_struct *sol,xycirc* xCircle,xycirc* yCircle)
{
  std::vector<double> soldata_dx;
  std::vector<double> soldata_dy;    
  soldata_struct soldata;
  d_struct_T expl_temp;
  double unchecked_XS[6],unchecked_YS[6],g_degC_data[5],g_lengthC_data[5],g_endTurn,g_startTurn;
  int g_degC_size[2],g_lengthC_size[2],g_lengthL_size[2],sol_size[2],b_nrEl,nrEl;

  //dummy values for replacing SMA
  //two feasible paths, the first should have a igher fitness value than the second
  //
  //coordinate transfer for guide points
  std::array<double,8> xxyy = getPointstoXY(x, model);
  //std::cout<<"F00: getPointstoXY done\n";
  //assign x_guide and y_guide correctly
  std::array<double,4> x_guide,y_guide;
  for(int i=0;i<4;i++){x_guide[i]=xxyy[i];y_guide[i]=xxyy[i+4];}
  
  //amount of points the path is divided into for collision avoidance
  nrEl = static_cast<int>(model->nrEl);
  //preallocation
  soldata_dx.resize(nrEl);
  soldata_dy.resize(nrEl);
  for (int i{0}; i < nrEl; i++) {
    soldata_dx[i] = 0.0;
    soldata_dy[i] = 0.0;
  }
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
  //std::cout<<"F00: checkPoints done\n";
  //printxsys(soldata_XS,soldata_YS); 

  //write outputs into g
  sol->internalpath.startTurn=g_startTurn;
  sol->internalpath.endTurn=g_endTurn;

  //initialize path with radii (so we can optimize the radii later)
  initialization(&(sol->internalpath), k, soldata_XS,soldata_YS);
  //std::cout<<"F00: initialization done\n";

  //optimize radii
  optimizationRad(&(sol->internalpath), k, soldata_XS,soldata_YS, xCircle, yCircle);
  //std::cout<<"F00: optimizationRad done\n";
  //std::cout<<"After radius optimization: \n radii: "; 
  //for(int i=0; i<size(sol->internalpath.r);i++){std::cout<<sol->internalpath.r[i]<<",";}
  //std::cout<<" "<<std::endl;
  //std::cout<<"sizes: "; 
  //for(int i=0; i<size(sol->internalpath.s);i++){std::cout<<sol->internalpath.s[i]<<",";}
  //std::cout<<" "<<std::endl;
  //std::cout<<"angles: "; 
  //for(int i=0; i<size(sol->internalpath.a);i++){std::cout<<sol->internalpath.a[i]<<",";}
  //std::cout<<" "<<std::endl;
  //std::cout<<"center points: "<<std::endl; 
  //for(int i=0; i<size(sol->internalpath.m);i++){std::cout<<sol->internalpath.m[i][0]<<", "<<sol->internalpath.m[i][1]<<std::endl;}
  //std::cout<<" "<<std::endl;

  int n=model->nrEl;
  double totalpathlength = 0;
  totalpathlength = discretize(&(sol->internalpath),sol,k, n,soldata_XS,soldata_YS); 
  //std::cout<<"F00: discretize done\n";
  //std::cout<<"amount of points: "<<sol->xx.size()<<std::endl;

  //get Violation Score
  std::vector<double> xwidth=model->xwidth;
  std::vector<double> ywidth=model->ywidth;
  int nobs = xwidth.size();
  double Violation=0;
  Violation = getViolationRect(&(sol->internalpath),model->xobs_rect,model->yobs_rect,xwidth,ywidth,nobs,sol,k,n,model->tolerance);
  //std::cout<<"F00: getViolationRect done!\n";
  //std::cout<<"Violation through rectangles: "<<Violation<<std::endl;

  //std::cout<<"coordinates of circle: \n"<<model->xobs_circ[0]<<","<<model->yobs_circ[0]<<"\nradius of circle: "<<model->robs_circ[0]<<std::endl;
  nobs=model->robs_circ.size();
  double Violation_Circ = getViolationCirc(&(sol->internalpath), model->xobs_circ,model->yobs_circ,model->robs_circ,nobs,sol,model->tolerance);
  //std::cout<<"F00: getViolationCirc done!\n";
//  std::cout<<"Violation through circles: "<<Violation_Circ<<std::endl;
  Violation+=Violation_Circ;
//get Violation Score for the Charging Station (Penalty for entering the charging station from the wrong side, I think)
//todo: Make insides of this function
double ViolationCharger=getViolationChargingStation(&(sol->internalpath), sol, model,k,n);
//std::cout<<"F00: getViolationChargingStation done!\n";
//std::cout<<"Violation through Charging Stations: "<<ViolationCharger<<std::endl;
Violation+=ViolationCharger; //if(model->AisChargingStation){ //  Violation += getViolationChargingStation(&(sol->internalpath), model->xs_CS,model->ys_CS, 0, model->SA, sol, k, n); //} //if(model->BisChargingStation){ //  Violation += getViolationChargingStation(&(sol->internalpath), model->xs_CS,model->ys_CS, 0, model->SA, sol, k, n);
//}

//calculate Energy demand
//todo: Implement this function properly
double L = getEnergy_v2(sol, totalpathlength);
//std::cout<<"getEnergy done!\n";

// Real Output
sol->L = L;
sol->Violation=Violation;
for (int i = 0; i < 6; ++i) {
sol->XS[i] = soldata_XS[i];
sol->YS[i] = soldata_YS[i];
}
//assigned xx,yy in discretize function
double res=sol->L +(1000 * sol->Violation)+1;
//std::cout<<"Violation: "<<sol->Violation<<std::endl;
//std::cout<<"Fitness Value: "<<sol->L<<std::endl;
//std::cout<<"F00 output: "<<res<<std::endl;
if(std::isnan(res))
{
  std::cout<<"F00 returned nan\n";
  std::cout<<"res: "<<res<<", sol->L: "<<sol->L<<", sol->Violation: "<<sol->Violation<<std::endl;
  std::cout<<"L: "<<L<<", totalpathlength: "<<totalpathlength<<std::endl;
}
else if(std::isinf(res))
{
  //std::cout<<"F00 returned inf\n";
  //std::cout<<"res: "<<res<<", sol->L: "<<sol->L<<", sol->Violation: "<<sol->Violation<<std::endl;
  //std::cout<<"L: "<<L<<", totalpathlength: "<<totalpathlength<<std::endl;
  //this is happening because totalpathlength is inf, which, in return, is happening bcs the guide points are stronomically wrong
  res=100000000 +rand01(); //add rand so that the sorting algorithm could determine which is the worst
}

return res;
}
#include <iostream>
#include <limits>
#include "custom_datatypes.h"
#include "custom_math_operations.cpp"
//#include "F00.cpp"
#include "SMA.cpp"
#include "format_path.cpp"
#include "getData.cpp"
#include "slimemould_func_prototypes.h"

bool main_alg(model_struct *model, path_struct *path)
{
    SMAsol_struct SMAsol;

    double bestPositions[6];
    double bestPositions2[6];
    double AllFitness2;
    double d;

    //part of model creation that handles charging stations
  // This function requires h.updock, so it's separated from the rest of the model creation
  //Since h is literally too big to pass, we will pass updock directly here:
    struct_dock updock;
    updock.L=findinh("updock.L")[0];
    updock.E=findinh("updock.E")[0];
    updock.T=findinh("updock.T")[0];
    struct_dock dock;
    dock.L=findinh("dock.L")[0];;
    dock.E=findinh("dock.E")[0];;
    dock.T=findinh("dock.T")[0];;

  if (model->AisChargingStation != 0.0) {
    model->xs_CS = model->xs;
    model->ys_CS = model->ys;
    model->A_is_CS = 1.0;
    AllFitness2 = model->SA;
    d = AllFitness2;
    b_cosd(d);
    model->x_CS1 = model->xs - d * 0.75;
    d = AllFitness2;
    b_sind(d);
    model->y_CS1 = model->ys - d * 0.75;
    d = AllFitness2;
    b_cosd(d);
    model->xs += d * updock.L;
    b_sind(AllFitness2);
    model->ys += AllFitness2 * updock.L;
    model->SA -= 180.0;
    if (model->SA < 0.0) {
      model->SA += 360.0;
    }
  } else {
    model->A_is_CS = 0.0;
  }
  if (model->BisChargingStation != 0.0) {
    AllFitness2 = model->EA;
    d = AllFitness2;
    b_cosd(d);
    model->xt_CS = model->xt - d * dock.L;
    d = AllFitness2;
    b_sind(d);
    model->yt_CS = model->yt - d * dock.L;
    model->B_is_CS = 1.0;
    d = AllFitness2;
    b_cosd(d);
    model->x_CS2 = model->xt_CS - d * 0.75;
    b_sind(AllFitness2);
    model->y_CS2 = model->yt_CS - AllFitness2 * 0.75;
  } else {
    model->B_is_CS = 0.0;
  }

//Define xCircle, yCircle => represent points on arc, separated from SMAsol bcs accessing values is easier that way
xycirc xCircle,yCircle;

//  run SMA Algorithm
for (int i{0}; i < 6; i++) {bestPositions[i] = 0;bestPositions2[i]=0;}
double a__1[180];
double inf = std::numeric_limits<double>::infinity();
SMA(model, bestPositions2, a__1,xCircle,yCircle);
std::cout<<"main_alg: Fitness of last iteration: "<<AllFitness2;
std::cout<<"main_alg: overwriting fitness with that of SMA\n";
for (int i{0}; i < 6; i++) {bestPositions[i] = bestPositions2[i];}

std::cout<<"SMA output: ";
for(int w=0; w<6;w++){std::cout<<bestPositions[w]<<",";} 

  //*Test Cases for testing F00
  //SMAsol_struct testsol = SMAsol;
  //test case 0: cuts through circle, but not through rectangles
//  bestPositions[0]=17;bestPositions[1]=58;bestPositions[2]=20;bestPositions[3]=72;bestPositions[4]=40;bestPositions[5]=5;
//  F00(bestPositions, model, &testsol,&xCircle,&yCircle);
//  if(testsol.Violation<=1e-5){
//    std::cout<<"Test Case 0 failed: Path goes through circle, but Violation==0\n";
//  }
//  else{
//    std::cout<<"Test Case 0 passed - Violation: "<<testsol.Violation<<std::endl;
//  }
//
//  testsol = SMAsol;
//  //Test Case 1: no collision, somewhat efficient path, but collides with charging station
//  bestPositions[0]=22;bestPositions[1]=50;bestPositions[2] = 20; bestPositions[3] = 70; bestPositions[4] = 18; bestPositions[5] = 6; //first and last: distance from start/end point, following SA and EA => two guide points. Middle: x and y of other two guide points
//  F00(bestPositions, model, &testsol,&xCircle,&yCircle);
//  if(testsol.Violation>=1e-5){
//    std::cout<<"Test Case 1 passed; path collides with charging station\n";
//  }
//  else{
//    std::cout<<"Test Case 1 failed; Violation<=0 even though path collides with Charging Station"<<std::endl;
//  }
//  testsol = SMAsol;
  //test case 2: no collision, but worse fitness than in test case 1
//  bestPositions[0]=46;bestPositions[1]=58;bestPositions[2] = 48; bestPositions[3] = 58; bestPositions[4] = 19; bestPositions[5] = 7;
//  F00(bestPositions, model, &testsol,&xCircle,&yCircle);
//  if(testsol.Violation>=1e-5){
//    std::cout<<"Test Case 2 failed: Violation>0 even though path should be feasible\n";
//  }
//  else{
//    std::cout<<"Test Case 2 passed "<<std::endl;
//  }
//  testsol = SMAsol;
  //Test Case 3:  cuts through rectangles
//  bestPositions[0]=20.044536;bestPositions[1]=6.05022;bestPositions[2] = 40.9792; bestPositions[3] = 93.9738;
//  bestPositions[4] = 60.4157; bestPositions[5] = 70.3338;
//F00(bestPositions, model, &testsol,&xCircle,&yCircle);
//if(testsol.Violation<=1e-5){
//  std::cout<<"Test Case 3 failed: Path goes through rectangles, but Violation==0\n";
//}
//else{
//  std::cout<<"Test Case 3 passed - Violation: "<<testsol.Violation<<std::endl;
//}
//  testsol = SMAsol;
//
//calculate final solution
//bestPositions[0]=86;bestPositions[1]=58;bestPositions[2] = 88; bestPositions[3] = 58; bestPositions[4] = 19; bestPositions[5] = 7;
F00(bestPositions, model, &SMAsol,&xCircle,&yCircle);

std::array<std::array<double,2>,4> centerPoints;
centerPoints[0]={SMAsol.internalpath.m[0][0],SMAsol.internalpath.m[0][1]};
centerPoints[1]={SMAsol.internalpath.m[1][0],SMAsol.internalpath.m[1][1]};
centerPoints[2]={SMAsol.internalpath.m[2][0],SMAsol.internalpath.m[2][1]};
centerPoints[3]={SMAsol.internalpath.m[3][0],SMAsol.internalpath.m[3][1]};

  // Convert path into different format for processing outside the code
  path_struct p=format_path(SMAsol,model,xCircle,yCircle,centerPoints);
  //std::cout<<"format_path done!"<<std::endl;

  //Assign value of p to path (pointer)
  if (path != nullptr) {
        *path = p; // Copy the contents of p into the object pointed to by path
    } else {
        std::cout << "Path Output: Invalid pointer to path_struct!\n";
    }

  //std::cout<<"main_alg done!";
  return 0;
}
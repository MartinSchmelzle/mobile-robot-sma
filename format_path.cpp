#include "custom_datatypes.h"

path_struct format_path(SMAsol_struct* SMAsol,model_struct *model,std::array<std::array<double,2>,4> centerPoints)
{
  path_struct path; 
  std::array<double,2> startPoint = {model->xs, model->ys};
  std::array<double,2> endPoint = {model->xt, model->yt};
  std::array<double,4> arc_startpoints_x;  
  std::array<double,4> arc_startpoints_y;
  std::array<double,4> arc_endpoints_x;
  std::array<double,4> arc_endpoints_y;
  std::array<double,4> arc_centerpoints_x;
  std::array<double,4> arc_centerpoints_y;
  std::array<int,4> arc_counterclockwise;
  std::array<double, 4> arc_radii;

  //write guide points into path
  for(int i=0;i<10;i++){
  path.XS[i]=SMAsol->XS[i];
  path.YS[i]=SMAsol->YS[i];}
  //get xCircle and yCircle into path
  for(int i=0; i<4;i++)
  {
    arc_startpoints_x[i] = SMAsol->internalpath.arc_startpoints_x[i];
    arc_startpoints_y[i] = SMAsol->internalpath.arc_startpoints_y[i];
    arc_endpoints_x[i] = SMAsol->internalpath.arc_endpoints_x[i];
    arc_endpoints_y[i] = SMAsol->internalpath.arc_endpoints_y[i];
    arc_centerpoints_x[i] = centerPoints[i][0];
    arc_centerpoints_y[i] = centerPoints[i][1];
    arc_radii[i]=SMAsol->internalpath.r[i];
  }


  //calculate clockwise or counterclockwise
  using row=std::array<double,4>;
  row vecStart_x;
  row vecStart_y;
  row vecEnd_x;
  row vecEnd_y;
  double linebefore_x;
  double linebefore_y;
  double startend_x;
  double startend_y;
  double anglechecker;
  double crossproduct_z;
  int i;
  for(i=0;i<4;i++){
  // Calculate vectors from center to start and end points
  vecStart_x[i] = arc_startpoints_x[i] - arc_centerpoints_x[i];
  vecStart_y[i] = arc_startpoints_y[i] - arc_centerpoints_y[i];
  vecEnd_x[i] = arc_endpoints_x[i] - arc_centerpoints_x[i];
  vecEnd_y[i] = arc_endpoints_y[i] - arc_centerpoints_y[i];
  //use straight line before the arc to determine whether angle >180°
  if(i==0){
    linebefore_x=arc_startpoints_x[i] - model->xs;
    linebefore_y=arc_startpoints_y[i] - model->ys;
  }
  else{
    linebefore_x=arc_startpoints_x[i] - arc_endpoints_x[i-1];
    linebefore_y=arc_startpoints_y[i] - arc_endpoints_y[i-1];
  }
  //if dot product between linebefore and startend >0: angle<180°
  startend_x = arc_endpoints_x[i] - arc_startpoints_x[i];
  startend_y = arc_endpoints_y[i] - arc_startpoints_y[i];
  anglechecker = linebefore_x*startend_x+linebefore_y*startend_y;

  // Calculate the cross product (we only care about z value of cross product)
  crossproduct_z=vecStart_x[i]*vecEnd_y[i] - vecStart_y[i] * vecEnd_x[i];

  // Check the sign of the z-component of the cross product
  if(crossproduct_z>0){arc_counterclockwise[i]=-1;}
  else if(crossproduct_z<0){arc_counterclockwise[i]=1;}
  else{arc_counterclockwise[i]=0;}

  //Account for angles above 180°
  //* we assume angle<180° in other parts of the code
  if(anglechecker<0){arc_counterclockwise[i]*=-1;}
  }
  //debugging stuff
  //for(int i=0;i<4;i++){std::cout<<"arc_startpoints_x["<<i<<"]: "<<arc_startpoints_x[i]<<std::endl;
  //std::cout<<"arc_startpoints_y["<<i<<"]: "<<arc_startpoints_y[i]<<std::endl;}

  //write values into path
  for(int i=0;i<2;i++){
    path.startPoint[i] = startPoint[i];
    path.endPoint[i] = endPoint[i];
  }
  for(int i=0;i<4;i++)
  {
    path.arc_startpoints_x[i]=arc_startpoints_x[i];
    path.arc_startpoints_y[i]=arc_startpoints_y[i];
    path.arc_endpoints_x[i]=arc_endpoints_x[i];
    path.arc_endpoints_y[i]=arc_endpoints_y[i];
    path.arc_center_points_x[i]=arc_centerpoints_x[i];
    path.arc_center_points_y[i]=arc_centerpoints_y[i];
    path.arc_counterclockwise[i]=arc_counterclockwise[i];
    path.arc_radii[i]=arc_radii[i];
  }
  return path;
}
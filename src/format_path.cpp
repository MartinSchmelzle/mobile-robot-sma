#include "../custom_datatypes.h"

path_struct format_path(SMAsol_struct* SMAsol,model_struct *model)
{
  std::array<std::array<double,2>,4> centerPoints;
  centerPoints[0]={SMAsol->internalpath.m[0][0],SMAsol->internalpath.m[0][1]};
  centerPoints[1]={SMAsol->internalpath.m[1][0],SMAsol->internalpath.m[1][1]};
  centerPoints[2]={SMAsol->internalpath.m[2][0],SMAsol->internalpath.m[2][1]};
  centerPoints[3]={SMAsol->internalpath.m[3][0],SMAsol->internalpath.m[3][1]};

  path_struct path; 
  std::array<int,4> arc_counterclockwise;

  //write guide points into path
  for(int i=0;i<10;i++){
  path.XS[i]=SMAsol->XS[i];
  path.YS[i]=SMAsol->YS[i];}

  //compute clockwise or counterclockwise
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
  vecStart_x[i] = SMAsol->internalpath.arc_startpoints_x[i] - centerPoints[i][0];
  vecStart_y[i] = SMAsol->internalpath.arc_startpoints_y[i] - centerPoints[i][1];
  vecEnd_x[i] = SMAsol->internalpath.arc_endpoints_x[i] - centerPoints[i][0];
  vecEnd_y[i] = SMAsol->internalpath.arc_endpoints_y[i] - centerPoints[i][1];
  //use straight line before the arc to determine whether angle >180°
  if(i==0){
    linebefore_x=SMAsol->internalpath.arc_startpoints_x[i] - model->xs;
    linebefore_y=SMAsol->internalpath.arc_startpoints_y[i] - model->ys;
  }
  else{
    linebefore_x=SMAsol->internalpath.arc_startpoints_x[i] - SMAsol->internalpath.arc_endpoints_x[i-1];
    linebefore_y=SMAsol->internalpath.arc_startpoints_y[i] - SMAsol->internalpath.arc_endpoints_y[i-1];
  }
  //if dot product between linebefore and startend >0: angle<180°
  //* we can neglect angles above 180°, they never occur
  startend_x = SMAsol->internalpath.arc_endpoints_x[i] - SMAsol->internalpath.arc_startpoints_x[i];
  startend_y = SMAsol->internalpath.arc_endpoints_y[i] - SMAsol->internalpath.arc_startpoints_y[i];
  anglechecker = linebefore_x*startend_x+linebefore_y*startend_y;

  // Calculate the cross product (we only care about z value of cross product)
  crossproduct_z=vecStart_x[i]*vecEnd_y[i] - vecStart_y[i] * vecEnd_x[i];

  // Check the sign of the z-component of the cross product
  if(crossproduct_z>0){arc_counterclockwise[i]=-1;}
  else if(crossproduct_z<0){arc_counterclockwise[i]=1;}
  else{arc_counterclockwise[i]=0;}

  //for(int i=0;i<4;i++){std::cout<<"arc_startpoints_x["<<i<<"]: "<<arc_startpoints_x[i]<<std::endl;
  //std::cout<<"arc_startpoints_y["<<i<<"]: "<<arc_startpoints_y[i]<<std::endl;}

  //write values into path
  path.startPoint[0] = model->xs;
  path.startPoint[1] = model->ys;
  path.endPoint[0] = model->xt;
  path.endPoint[1] = model->yt;

  for(int i=0;i<4;i++)
  {
    path.arc_startpoints_x[i] = SMAsol->internalpath.arc_startpoints_x[i];
    path.arc_startpoints_y[i] = SMAsol->internalpath.arc_startpoints_y[i];
    path.arc_endpoints_x[i] = SMAsol->internalpath.arc_endpoints_x[i];
    path.arc_endpoints_y[i] = SMAsol->internalpath.arc_endpoints_y[i];
    path.arc_center_points_x[i] = centerPoints[i][0];
    path.arc_center_points_y[i] = centerPoints[i][1];
    path.arc_radii[i]=SMAsol->internalpath.r[i];
  }
  }
  return path;
}
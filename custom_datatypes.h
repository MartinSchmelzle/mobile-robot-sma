#ifndef MY_STRUCTS_H
#define MY_STRUCTS_H

#include <vector>
#include <array>

using xycirc = std::vector<std::vector<double>>;

struct turned_rectangle{
    double xpos;
    double ypos;
    double width;
    double height;
    double angle;

};

struct path_struct {
    std::array<double, 2> startPoint;
    std::array<double, 2> endPoint;
    std::array<double, 4> arc_startpoints_x;
    std::array<double, 4> arc_startpoints_y;
    std::array<double, 4> arc_endpoints_x;
    std::array<double, 4> arc_endpoints_y;
    std::array<double, 4> arc_center_points_x;
    std::array<double, 4> arc_center_points_y;
    std::array<int, 4> arc_counterclockwise;
    double XS[10];
    double YS[10];
};

struct struct_internalpath {
    std::vector<double> degC;
    std::vector<double> lengthC;
    std::vector<double> lengthL;
    double startTurn;
    double endTurn;
    std::vector<std::vector<double>> xCircle;
    std::vector<std::vector<double>> yCircle;
    std::vector<double> s;
    std::vector<double> a;
    std::vector<double> r;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<std::vector<double>> m;
    std::vector<double> dis;
};

struct struct_dock {
    double L;
    double E;
    double T;
};

struct struct_updock {
    double L;
    double E;
    double T;
};

struct model_struct {
    double AisChargingStation;
    double BisChargingStation;
    double SA;
    double EA;
    double tolerance;
    double xs_CS;
    double ys_CS;
    double A_is_CS;
    double x_CS1;
    double y_CS1;
    double xt_CS;
    double yt_CS;
    double B_is_CS;
    double x_CS2;
    double y_CS2;
    double xs;
    double ys;
    double xt;
    double yt;
    std::vector<double> xobs_rect;
    std::vector<double> yobs_rect;
    std::vector<double> xwidth;
    std::vector<double> ywidth;
    std::vector<double> xobs_circ;
    std::vector<double> yobs_circ;
    std::vector<double> robs_circ;
    double nrEl;
    double facV;
    double n;
    double lb;
    double ub;
    //newly implemented charging-related obstacles
    std::vector<turned_rectangle> blocked_area_down; //lower rectangles
    std::vector<turned_rectangle> blocked_area_left; //left rectangles
    std::vector<turned_rectangle> blocked_area_up; //upper rectangles
};

struct SMAsol_struct {
    double Violation;
    double L;
    double L2;
    double XS[10];
    double YS[10];
    struct_internalpath internalpath;
    std::vector<double> xx;
    std::vector<double> yy;
    // Add other members as needed
};

//define RGB value
struct RGBColor {
    unsigned red;
    unsigned green;
    unsigned blue;
};

struct soldata_struct {
  double XS[10];
  double YS[10];
  double L1;
  double T1;
  double L2;
  double L;
  double T;
  std::array<double, 2> dx;
  std::array<double, 2> dy;
  double Violation;
  bool IsFeasible;
  double length;
  std::array<double,2>xx;
  std::array<double,2>yy;
  std::array<std::array<double, 2>, 7> TS;
  std::array<double, 2> tt;
};


//I'm not sure what exactly this struct does, maybe I'll find out later
struct d_struct_T {
  double XS[10];
  double YS[10];
  double L1;
  double T1;
  double L2;
  double L;
  double T;
  std::array<double, 2> dx;
  std::array<double, 2> dy;
  double Violation;
  bool IsFeasible;
  double length;
  std::array<double, 2> xx;
  std::array<double, 2> yy;
};

#endif // MY_STRUCTS_H

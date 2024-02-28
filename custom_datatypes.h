#ifndef MY_STRUCTS_H
#define MY_STRUCTS_H

#include <vector>
#include <array>

//This struct is supposed to represent algorithm configuration
//todo: expand on it, e.g. with variable amount of arcs, variable amount of search agents etc.
struct alg_data{
    int max_iter;
    double distance_discrete;
};

struct output_SMA{
    std::array<double,6> waypoints;
    float Vio;
    double res;
};

struct struct_fitness{
  float Vio; //Violation
  double L; //Energy demand
  double res; //resulting combined fitness
};

struct kinematics_time{
    double a_acc; //represents acceleration of robot
    double v_max; //represents max speed
    double v_desired; //we don't want to be exactly on the maximum speed
    double a_break; //represents acceleration during deceleration. Must be <0!
    double delta_t; //time interval
    double P_idle; //power consumption with v=0 and a=0
    double mass; //weight of the mobile robot
    float time_array[10000];
    float a_array[10000] = {};
    float v_array[10000];
    float x_array[10000];
    float y_array[10000];
    float s_array[10000];
    float w_array[10000];
    float dwdt_array[10000];
};

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
    std::array<double, 4> arc_radii;
    std::array<int, 4> arc_counterclockwise;
    double XS[10];
    double YS[10];
};

struct struct_internalpath {
    double startTurn;
    double endTurn;
    std::array<double, 4> arc_startpoints_x;
    std::array<double, 4> arc_startpoints_y;
    std::array<double, 4> arc_endpoints_x;
    std::array<double, 4> arc_endpoints_y;
    std::vector<double> s;
    std::vector<double> a;
    std::vector<double> r;
    std::vector<std::vector<double>> m;
    std::vector<double> dis;
    std::array<std::array<double,2>,5> lineVector; //This is used by discretize and getEnergy. Important: This is a unit vector!
    std::array<double,5> lineLength;
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
    float Violation;
    double L;
    double L2;
    std::array<double,10> XS;
    std::array<double,10> YS;
    struct_internalpath internalpath;
    std::vector<float> xx;
    std::vector<float> yy;
    // Add other members as needed
};

using arr2 = std::array<int, 2>; //This makes the code less bulky. Unsigned datatype prevents negative values of position and size.

//define RGB value
struct RGBColor {
    unsigned red;
    unsigned green;
    unsigned blue;
};

//These lines are responsible for representing color => names, RGB values, mappings
enum class col {
    red,
    green,
    blue,
    orange,
    yellow,
    purple,
    cyan,
    magenta,
    pink,
    teal,
    lime,
    indigo,
    gold,
    silver,
    maroon,
    navy,
    brown,
    black,
    white,
    // Add more colors as needed...
};

#endif // MY_STRUCTS_H

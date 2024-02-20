#include<vector>
#include<array>

struct struct_fitness{
  double Vio; //Violation
  double L; //Energy demand
  double res; //resulting combined fitness
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
};

struct struct_internalpath {
    std::array<std::array<double, 2>, 5> degC;
    std::array<std::array<double, 2>, 5> lengthC;
    std::array<std::array<double, 2>, 6> lengthL;
    double startTurn;
    double endTurn;
    std::vector<std::vector<double>> xCircle; //*This used to be a cell wrap; I might need to rework this later
    std::vector<std::vector<double>> yCircle;
    std::array<std::array<double, 2>, 5> s;
    std::array<std::array<double, 2>, 5> a;
    std::array<std::array<double, 2>, 5> r;
    std::array<std::array<double, 2>, 5> alpha;
    std::array<std::array<double, 2>, 5> beta;
    std::vector<std::vector<double>> m;
    std::array<std::array<double, 2>, 6> dis;

    /*struct_internalpath() { //ChatGPT has automatically inserted this, I might need it later
        // Resize the vectors to the appropriate sizes
        xCircle.resize(5, std::vector<cell_wrap_0>(2));
        yCircle.resize(5, std::vector<cell_wrap_0>(2));
        m.resize(10, std::vector<double>(2));
    }*/
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
};


//*soldata
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

//time to take care of h
//the one for h itself is in the bottom 

// Type Definitions
struct struct3_T {
  double L;
  double T;
  double E;
  double D;
};

struct struct20_T {
  double L;
  double T;
  double E;
  double a;
};

struct struct_dock {
  double L;
  double E;
  double T;
};

struct struct4_T {
  double L[53595];
  double T[53595];
  double E[53595];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c100 {
  struct3_T Acc;
  struct3_T Slw;
  struct4_T Con;
};

struct struct6_T {
  double L[53496];
  double T[53496];
  double E[53496];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c10 {
  struct3_T Acc;
  struct3_T Slw;
  struct6_T Con;
};

struct struct8_T {
  double L[57050];
  double T[57050];
  double E[57050];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c1 {
  struct3_T Acc;
  struct3_T Slw;
  struct8_T Con;
};

struct struct10_T {
  double L[38208];
  double T[38208];
  double E[38208];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c20 {
  struct3_T Acc;
  struct3_T Slw;
  struct10_T Con;
};

struct struct12_T {
  double L[80536];
  double T[80536];
  double E[80536];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c2 {
  struct3_T Acc;
  struct3_T Slw;
  struct12_T Con;
};

struct struct14_T {
  double L[110674];
  double T[110674];
  double E[110674];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c3 {
  struct3_T Acc;
  struct3_T Slw;
  struct14_T Con;
};

struct struct16_T {
  double L[141760];
  double T[141760];
  double E[141760];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c4 {
  struct3_T Acc;
  struct3_T Slw;
  struct16_T Con;
};

struct struct18_T {
  double L[208581];
  double T[208581];
  double E[208581];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c6 {
  struct3_T Acc;
  struct3_T Slw;
  struct18_T Con;
};

struct struct21_T {
  double L[377802];
  double T[377802];
  double E[377802];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
};

struct struct_s {
  struct20_T Acc;
  struct20_T Slw;
  struct21_T Con;
};

struct struct25_T {
  double L[33171];
  double T[33171];
  double E[33171];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c40 {
  struct3_T Acc;
  struct3_T Slw;
  struct25_T Con;
};

struct struct27_T {
  double L[60795];
  double T[60795];
  double E[60795];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c5 {
  struct3_T Acc;
  struct3_T Slw;
  struct27_T Con;
};

struct struct29_T {
  double L[31331];
  double T[31331];
  double E[31331];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c7 {
  struct3_T Acc;
  struct3_T Slw;
  struct29_T Con;
};

struct struct31_T {
  double L[34904];
  double T[34904];
  double E[34904];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c8 {
  struct3_T Acc;
  struct3_T Slw;
  struct31_T Con;
};

struct struct33_T {
  double L[44849];
  double T[44849];
  double E[44849];
  double v_Max;
  double v_Mean;
  double b_EL;
  double b_TL;
  double b_ED;
  double b_TD;
};

struct struct_c9 {
  struct3_T Acc;
  struct3_T Slw;
  struct33_T Con;
};

struct struct_b {
  double TL[13];
  double TLs;
  double EL[13];
  double ELs;
  double TD[14];
  double ED[14];
};

struct struct_R {
  double l[13];
  double d[14];
};

struct struct_Acc {
  double D[14];
  double T[14];
  double E[14];
  double L[14];
  double a[14];
  double v_Max[14];
  double v_Mean[14];
};

struct struct_Slw {
  double D[14];
  double T[14];
  double E[14];
  double L[14];
  double a[14];
};

struct struct23_T {
  double L[22287];
  double T[22287];
  double E[22287];
  double v_Max;
  double v_Mean;
  double b_ED;
  double b_TD;
};

struct struct_c0 {
  struct3_T Acc;
  struct3_T Slw;
  struct23_T Con;
};

struct struct_h {
  struct_c100 c100;
  struct_c10 c10;
  struct_c1 c1;
  struct_c20 c20;
  struct_c2 c2;
  struct_c3 c3;
  struct_c4 c4;
  struct_c6 c6;
  struct_s s;
  struct_c0 c0;
  struct_c40 c40;
  struct_c5 c5;
  struct_c7 c7;
  struct_c8 c8;
  struct_c9 c9;
  struct_dock dock;
  struct_dock updock;
  struct_b b;
  struct_R R;
  struct_Acc Acc;
  struct_Slw Slw;
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
  std::vector<double> xobs;
  std::vector<double> yobs;
  std::vector<double> xwidth;
  std::vector<double> ywidth;
  double nrEl;
  double facV;
  double n;
  double lb;
  double ub;
};

//define RGB value
struct RGBColor {
    unsigned red;
    unsigned green;
    unsigned blue;
};
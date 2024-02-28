#pragma once
#include <cstdint>
#include "custom_datatypes.h"
#include <string>
//main
bool mr_sma_alg(model_struct *model, path_struct *path);

using point = std::array<double, 2>;
//custom math operations
std::vector<point> calculatePointsOnArc(const point start, const point center, const point end,
                                        int direction, double distance_betw_points,double addtonext, double &arcLength); 
std::vector<point> calculatePointsOnArc_legacy(const point start, const point center, const point end,
                                        int direction, double distance_betw_points,double addtonext, double &arcLength);
std::vector<double> linspace(double start, double end, int numPoints);

std::vector<double> diff(const std::vector<double>& input);

int sign(double n);

double rand01(); 
int randi(int lower_bound, int upper_bound);

void sort(double x[30], int idx[30]);

struct_fitness F00(double x[6], const model_struct *model, SMAsol_struct *sol, kinematics_time* kin, const bool accuracy);

double euclid_dist(std::array<double,2> start, std::array<double,2> end);

double arc_cross(std::array<double,2> start, std::array<double,2> center, std::array<double,2> end);

int arc_direction(std::array<double,2> start, std::array<double,2> center, std::array<double,2> end);

double getArcLength(std::array<double,2> start, std::array<double,2> center, std::array<double,2> end);

double b_cosd(double x);
double b_sind(double x);
//F00_functions
std::array<double,8> getPointstoXY(const double x[6],const model_struct* model);

void push_to_left(std::array<double, 10> &arr, const int index);

void push_to_left(std::array<double, 9> &arr,const int index);
void checkPoints(const double XS[6], const double YS[6], SMAsol_struct* sol, uint8_t &k);

void getAngRad(const std::array<double, 10> XS,const std::array<double, 10> YS,const int i,const double g_s, double& g_a, double& g_r);

void findMid(const std::array<double, 10> XS,const std::array<double, 10> YS, double dis[],const double g_r,const double g_s,const int i, double* m0, double* m1);

std::array<double,2> findMidAlphaBeta(const std::array<double, 10> XS,const std::array<double, 10> YS,const std::vector<double> dis, const double g_r,const double g_s,const int i, 
std::array<double,2> &startpoint, std::array<double,2> &endpoint);

void initialization(SMAsol_struct* sol,const int k);

void optimizationRad(SMAsol_struct* sol,const int k);

void printpoint(const point P);

double lastnonnan(const std::array<double,10> XS);

double discretize(SMAsol_struct* sol, const int k, const int n,const double distance_discrete); 

float getViolationRect(const model_struct* model,const SMAsol_struct* sol,const int k);

float getViolationCircBoundary(const model_struct* model,const SMAsol_struct* sol);

void printxsys(const std::array<double, 10> XS, const std::array<double, 10> YS);

bool pointInRotatedRectangle(const double rectX, const double rectY, const double rectWidth, const double rectHeight, double rotationAngle, const double px, const double py);

float getViolationChargingStation(const SMAsol_struct* sol, const model_struct* model, const int k);

float getSmallRadVio(const SMAsol_struct* sol);

void setArrayTo0(float* array, const size_t size);
//getEnergy
int find_i_start(const float array[10000]);

void printKin(const kinematics_time& kinematics);

double get_decel_dist(kinematics_time* kin, double length_last);

void get_kinematics_line(kinematics_time* kin, std::array<double,2> lineVector, const double lineLength, const bool breakatend);

void get_kinematics_arc(kinematics_time* kin, const double arc_break_length, const double centerPoint[2], const double endpoint[2], const double radius);

double getEnergy_v2(SMAsol_struct* sol, kinematics_time* kin);
//getData
std::vector<double> findinh(const std::string& searchString);
//SMA
output_SMA SMA(const model_struct *model, SMAsol_struct* sol1, kinematics_time* kin, alg_data* alg);
//format_path
path_struct format_path(SMAsol_struct* SMAsol,model_struct *model);
//processpath
void printPathStruct(const path_struct& path);
int drawpathintosvg(path_struct path,unsigned scale);
std::vector<point> calculatePointsOnArc_legacy(const point start, const point center, const point end,
                                        int direction, double distance_betw_points,double addtonext, double &arcLength);
std::vector<point> pathtopointcloud(path_struct path, double distance_betw_points, model_struct model);
int drawpointcloud_svg(std::vector<point> pointcloud, unsigned scale);
int pathtocsv(path_struct path);
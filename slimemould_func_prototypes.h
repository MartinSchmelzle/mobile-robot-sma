#pragma once
#include "custom_datatypes.h"
//main
bool main_alg(model_struct *model, path_struct *path, const unsigned max_iter);

using point = std::array<double, 2>;
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

inline double arc_cross(std::array<double,2> start, std::array<double,2> center, std::array<double,2> end);

inline int arc_direction(std::array<double,2> start, std::array<double,2> center, std::array<double,2> end);

double getArcLength(std::array<double,2> start, std::array<double,2> center, std::array<double,2> end);

double b_cosd(double x);
double b_sind(double x);
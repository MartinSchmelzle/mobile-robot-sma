#pragma once
#include "custom_datatypes.h"
//main
bool main_alg(model_struct *model, path_struct *path);

using point = std::array<double, 2>;
std::vector<point> calculatePointsOnArc(const point start, const point center, const point end,
                                        int direction, double distance_betw_points,double addtonext, double &arcLength); 

std::vector<double> linspace(double start, double end, int numPoints);

std::vector<double> diff(const std::vector<double>& input);

int sign(double n);

double rand01(); 
int randi(int lower_bound, int upper_bound);

void sort(double x[30], int idx[30]);

double F00(double x[6], const model_struct *model, SMAsol_struct *sol,xycirc* xCircle,xycirc* yCircle);
//#include "custom_math_operations.cpp"
#include <iostream>
#include <numeric>
#include "custom_datatypes.h"
#include <cassert>
#include <algorithm>
#include <cmath>
#include "slimemould_func_prototypes.h"
#include <cstdint>

//sets all entries of a given C-type array to 0
//void setArrayToNaN(double* array, size_t size) {
//    for (size_t i = 0; i < size; ++i) {
//        array[i] = std::nan("");
//    }
//}
void setArrayTo0(float* array, const size_t size) {
    for (size_t i = 0; i < size; ++i) {
        array[i] = 0;//std::nan("");
    }
}

//This function does some coordinate transformations.
//todo: Change how SA and EA are processed => The angles go clockwise, but should go counterclockwise => 90° and 270° are reversed
std::array<double,8> getPointstoXY(const double x[6],const model_struct* model)
{
  double xx[4],yy[4];
  double alpha, beta, RA_tmp,a,absxk,d,scale,t,v_rotatedA_idx_0,v_rotatedB_idx_0,y;

  alpha = 0.017453292519943295 * model->SA; 
  beta = 0.017453292519943295 * model->EA;
  //  Rotationsmatrix
  RA_tmp = std::sin(alpha);
  alpha = std::cos(alpha);
  const double RB_tmp = std::sin(beta);
  beta = std::cos(beta);
  //  Einheitsvektor vor der Drehung
  //  Berechnung des Einheitsvektors nach der Drehung
  v_rotatedA_idx_0 = alpha + -RA_tmp * 0.0;
  d = RA_tmp + alpha * 0.0;
  //  Normalisierung des resultierenden Vektors, um sicherzustellen, dass er ein
  //  Einheitsvektor bleibt
  a = std::abs(x[0]);
  scale = 3.3121686421112381E-170;
  v_rotatedB_idx_0 = beta + -RB_tmp * 0.0;
  absxk = std::abs(v_rotatedA_idx_0);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  alpha = RB_tmp + beta * 0.0;
  absxk = std::abs(d);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  y = scale * std::sqrt(y);
  beta = std::abs(x[5]);
  scale = 3.3121686421112381E-170;
  v_rotatedA_idx_0 = a * v_rotatedA_idx_0 / y;
  absxk = std::abs(v_rotatedB_idx_0);
  if (absxk > 3.3121686421112381E-170) {
    RA_tmp = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    RA_tmp = t * t;
  }
  absxk = std::abs(alpha);
  if (absxk > scale) {
    t = scale / absxk;
    RA_tmp = RA_tmp * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    RA_tmp += t * t;
  }
  RA_tmp = scale * std::sqrt(RA_tmp);
  v_rotatedB_idx_0 = beta * v_rotatedB_idx_0 / RA_tmp;
  xx[0] = model->xs + v_rotatedA_idx_0;
  xx[3] = model->xt + v_rotatedB_idx_0;
  yy[0] = model->ys + a * d / y;
  xx[1] = x[1];
  yy[1] = x[2];
  xx[2] = x[3];
  yy[2] = x[4];
  yy[3] = model->yt + beta * alpha / RA_tmp;
  std::array<double,8> returnvalue = {xx[0],xx[1],xx[2],xx[3],yy[0],yy[1],yy[2],yy[3]};
  return returnvalue;
}

//pushes numbers in a fixed-size array to the left, deleting a number
//dependency for checkPoints
void push_to_left(std::array<double, 10> &arr, const int index) {
    for (int i = index; i < arr.size() - 1; ++i) {
        arr[i] = arr[i + 1];
    }
}

void push_to_left(std::array<double, 9> &arr,const int index) {
    for (int i = index; i < arr.size() - 1; ++i) {
        arr[i] = arr[i + 1];
    }
}

//checks if guide points are in order (if they are too close to one another, if they are in a straight line, if they lead nowhere)
// Since MATLAB Coder forced me to use fixed size arrays, this function features a weird, C-like programming style.
//todo: Replace fixed-size arrays with vectors
void checkPoints(const double XS[6], const double YS[6], std::array<double, 10> &XS_out,
                 std::array<double, 10> &YS_out, double &startTurn, double &endTurn, uint8_t &k) {
    std::array<double, 6> temp = {XS[5], YS[5]};
    k = 6;
    XS_out.fill(NAN);
    YS_out.fill(NAN);
    for (int i = 0; i < k; ++i) {
        XS_out[i] = XS[i];
        YS_out[i] = YS[i];
    }
    double e = 1e-2;
    int p = 1;
    startTurn = 0;
    endTurn = 0;
    while (p) {
        p = 0;
        for (int i = 1; i < k; ++i) {
            if ((std::abs(XS_out[i] - XS_out[i - 1]) < e) && (std::abs(YS_out[i] - YS_out[i - 1]) < e)) {
                push_to_left(XS_out, i);
                push_to_left(YS_out, i);
                k--;
                p = 1;
                if (i == 1) {
                    startTurn = 1;
                } else if (i == k) {
                    endTurn = 1;
                }
            }
        }
        for (int i = 2; i < k - 1; ++i) {
            if ((std::abs(XS_out[i] - XS_out[i - 2]) < e) && (std::abs(YS_out[i] - YS_out[i - 2]) < e)) {
                push_to_left(XS_out, i);
                push_to_left(YS_out, i);
                k--;
                p = 1;
            }
        }
    }
    std::array<double, 9> XSV_out, YSV_out;
    for (int i = 0; i < k - 1; ++i) {
        XSV_out[i] = XS_out[i + 1] - XS_out[i];
        YSV_out[i] = YS_out[i + 1] - YS_out[i];
        double fac = (1.0 / std::sqrt(XSV_out[i] * XSV_out[i] + YSV_out[i] * YSV_out[i]));
        XSV_out[i] *= fac;
        YSV_out[i] *= fac;
    }
    p = 1;
    e = 1e-5;
    while (p) {
        p = 0;
        for (int i = 1; i < k - 1; ++i) {
            if (((std::abs(XSV_out[i] - XSV_out[i - 1]) < e) && (std::abs(YSV_out[i] - YSV_out[i - 1]) < e)) ||
                ((std::abs(XSV_out[i] + XSV_out[i - 1]) < e) && (std::abs(YSV_out[i] + YSV_out[i - 1]) < e))) {
                push_to_left(XS_out, i);
                push_to_left(YS_out, i);
                push_to_left(XSV_out, i);
                push_to_left(YSV_out, i);
                k--;
                p = 1;
            }
        }
    }
    if (XS_out[k - 1] != temp[0] || YS_out[k - 1] != temp[1]) {
        k++;
        XS_out[k - 1] = temp[0];
        YS_out[k - 1] = temp[1];
    }
}

//needed for initialization; this function takes in s (distance from straight line added to arc) and computes r (radius) and a (some kind of angle)
void getAngRad(const std::array<double, 10> XS,const std::array<double, 10> YS,const int i,const double g_s, double& g_a, double& g_r) {
    // Vektoren AB und BC
    double vectorBA[2] = { -(XS[i + 1] - XS[i]), -(YS[i + 1] - YS[i]) };
    double vectorBC[2] = { XS[i + 2] - XS[i + 1], YS[i + 2] - YS[i + 1] };

    // Berechne das Skalarprodukt von AB und BC
    double dotProduct = vectorBA[0] * vectorBC[0] + vectorBA[1] * vectorBC[1];

    // Berechne die Längen der Vektoren
    // Berechne den Winkel in Bogenmaß (Radian) zwischen AB und BC
    double lengthBA = std::sqrt(vectorBA[0] * vectorBA[0] + vectorBA[1] * vectorBA[1]);
    double lengthBC = std::sqrt(vectorBC[0] * vectorBC[0] + vectorBC[1] * vectorBC[1]);
    
    g_a = std::acos(dotProduct / (lengthBA * lengthBC));

    g_r = g_s * std::tan(g_a / 2.0);
    g_a *= 180 / 3.14159;
    //add a few asserts to prevent g_r==0
    assert(g_s!=0);
    assert(std::tan(g_a)!=0);
    assert(g_r>=0);
}

//version without alpha and beta => called by initialization
void findMid(const std::array<double, 10> XS,const std::array<double, 10> YS, double dis[],const double g_r,const double g_s,const int i, double* m0, double* m1) {
    // Mittelpunkt finden
    
    double vecAB[2] = {XS[i + 1] - XS[i], YS[i + 1] - YS[i]};
    double alpha, beta;
    
    if (g_r == 0) {
        alpha = 0;
        beta = 0;
        *m0 = XS[i] + (dis[i] - g_s) / dis[i] * vecAB[0];
        *m1 = YS[i] + (dis[i] - g_s) / dis[i] * vecAB[1];
    } else {
        // Kreuzprodukt finden (senkrechter Vektor auf AB)
        double cross_product[2] = {-vecAB[1], vecAB[0]};
        
        // Vektor auf Größe des Radius skalieren
        double a = (g_r / std::sqrt(cross_product[0] * cross_product[0] + cross_product[1] * cross_product[1])) * cross_product[0];
        double b = (g_r / std::sqrt(cross_product[0] * cross_product[0] + cross_product[1] * cross_product[1])) * cross_product[1];
    
        *m0 = XS[i] + (dis[i] - g_s) / dis[i] * vecAB[0] + a;
        *m1 = YS[i] + (dis[i] - g_s) / dis[i] * vecAB[1] + b;
    
        // Kreuzprodukt finden (senkrechter Vektor auf AB)
        double vecBC[2] = {XS[i + 2] - XS[i + 1], YS[i + 2] - YS[i + 1]};
        double cross_product_bc[2] = {-vecBC[1], vecBC[0]};
        
        // Vektor auf Größe des Radius skalieren
        double c = (g_r / std::sqrt(cross_product_bc[0] * cross_product_bc[0] + cross_product_bc[1] * cross_product_bc[1])) * cross_product_bc[0];
        double d = (g_r / std::sqrt(cross_product_bc[0] * cross_product_bc[0] + cross_product_bc[1] * cross_product_bc[1])) * cross_product_bc[1];
    
        double temp[2] = {XS[i + 1] + g_s / dis[i + 1] * vecBC[0] + c, YS[i + 1] + g_s / dis[i + 1] * vecBC[1] + d};
    
        if (std::sqrt(std::pow(temp[0] - *m0, 2) + std::pow(temp[1] - *m1, 2)) > 1e-5) {
            a = -a;
            b = -b;
            *m0 = *m0 + 2 * a;
            *m1 = *m1 + 2 * b;
        }
    }
}

//includes calculation of start and end point => called by optimizationRad
std::array<double,2> findMidAlphaBeta(const std::array<double, 10> XS,const std::array<double, 10> YS,const std::vector<double> dis, const double g_r,const double g_s,const int i, 
std::array<double,2> &startpoint, std::array<double,2> &endpoint) {
   // Mittelpunkt finden
   double middle0, middle1;
   
   double vecAB[2] = {XS[i + 1] - XS[i], YS[i + 1] - YS[i]};
    double vecBC[2] = {XS[i + 2] - XS[i + 1], YS[i + 2] - YS[i + 1]};
   if (g_r == 0) {
       middle0 = XS[i] + (dis[i] - g_s) / dis[i] * vecAB[0];
       middle1 = YS[i] + (dis[i] - g_s) / dis[i] * vecAB[1];
   } else {
       // Kreuzprodukt finden (senkrechter Vektor auf AB)
       double cross_product[2] = {-vecAB[1], vecAB[0]};
       
       // Vektor auf Größe des Radius skalieren
       double abs_cross=std::sqrt(cross_product[0] * cross_product[0] + cross_product[1] * cross_product[1]);
       double ax = (g_r / abs_cross) * cross_product[0];
       double ay = (g_r / abs_cross) * cross_product[1];
   
       middle0 = XS[i] + (dis[i] - g_s) / dis[i] * vecAB[0] + ax;
       middle1 = YS[i] + (dis[i] - g_s) / dis[i] * vecAB[1] + ay;
    
       // Kreuzprodukt finden (senkrechter Vektor auf AB)
       double cross_product_bc[2] = {-vecBC[1], vecBC[0]};
    
       // Vektor auf Größe des Radius skalieren
       double bx = (g_r / std::sqrt(cross_product_bc[0] * cross_product_bc[0] + cross_product_bc[1] * cross_product_bc[1])) * cross_product_bc[0];
       double by = (g_r / std::sqrt(cross_product_bc[0] * cross_product_bc[0] + cross_product_bc[1] * cross_product_bc[1])) * cross_product_bc[1];
   
       double temp[2] = {XS[i + 1] + g_s / dis[i + 1] * vecBC[0] + bx, YS[i + 1] + g_s / dis[i + 1] * vecBC[1] + by};
   

       if (std::sqrt(std::pow(temp[0] - middle0, 2) + std::pow(temp[1] - middle1, 2)) > 1e-5) {
           ax = -ax;
           ay = -ay;
           bx=-bx;
           by=-by;
           middle0 += 2 * ax;
           middle1 += 2 * ay;
       }
   }
    //extensively tested; correct until here
    //calculate start and end point of arc
    double len = sqrt(vecAB[0]*vecAB[0]+vecAB[1]*vecAB[1]);
    vecAB[0]/=len;
    vecAB[1]/=len;
    len=sqrt(vecBC[0]*vecBC[0]+vecBC[1]*vecBC[1]);
    startpoint[0]=XS[i+1] -g_s*vecAB[0];
    startpoint[1]=YS[i+1] -g_s*vecAB[1];
    vecBC[0]/=len;
    vecBC[1]/=len;
    endpoint[0]=XS[i+1] +g_s*vecBC[0];
    endpoint[1]=YS[i+1] +g_s*vecBC[1];
  
  std::array<double,2> result = {middle0,middle1};
  return result;
}

//initialize arcs with radii
void initialization(SMAsol_struct* sol,const int k)
{
  double g_s[k-2],g_a[k-2],g_r[k-2],alpha[k-2],beta[k-2],m[2][k-2],dis[k-1];
  //get distances between guide points
  double di;
  for (int i = 0; i < k - 1; ++i) {
    di = std::sqrt(std::pow(sol->XS[i + 1] - sol->XS[i], 2) + std::pow(sol->YS[i + 1] - sol->YS[i], 2));
    dis[i] = di; //distances between guide points
  }
  //for-loop for initialization
  for(int i=0;i<k-2;i++)
  {
    g_s[i] = std::min(dis[i],dis[i+1])*0.25; //in MATLAB Code: 0.25 is called 'start'; amount of space on lines used for arcs
    getAngRad(sol->XS,sol->YS,i,g_s[i],g_a[i],g_r[i]); //get angles and radii for a given g_s
    findMid(sol->XS,sol->YS,dis,g_r[i],g_s[i],i,&m[0][i],&m[1][i]); //find middle point
    sol->internalpath.a.push_back(g_a[i]); //write angles, sizes, center points, radii and distances into g
    sol->internalpath.s.push_back(g_s[i]);
    sol->internalpath.r.push_back(g_r[i]);
    sol->internalpath.m.push_back({m[0][i],m[1][i]});
    assert(g_r[i]>0);
  }
  for(int i=0;i<k-1;i++){sol->internalpath.dis.push_back(dis[i]);}
}

void optimizationRad(SMAsol_struct* sol,const int k) 
{
    struct_internalpath* g = &(sol->internalpath);
    double it, temp; int j; 
    const float lb=4, ub=9;

    for (int i = 0; i < k - 2; i++) {
        if (std::abs(g->r[i]) < lb) {
            j = 1;
            
            if ((i == 0) && (g->dis[0] == std::min(g->dis[0], g->dis[1]))) {
                it = 30;
                temp = g->dis[1] / 2;
            } else if ((i == k - 3) && (g->dis[k - 1] == std::min(g->dis[k - 2], g->dis[k - 1]))) {
                it = 30;
                temp = g->dis[k - 2]/ 2; 
            } else {
                it = 10;
                temp = std::numeric_limits<double>::infinity();
            }
            
            while ((std::abs(g->r[i]) < (lb - j * 0.1)) && (j <= it) && (g->s[i] < temp)) {
                j++;
                g->s[i] += g->s[i] / 10;
                g->s[i] = std::min(g->s[i], temp);
                getAngRad(sol->XS, sol->YS, i, g->s[i], g->a[i], g->r[i]);
            }
        } else if (std::abs(g->r[i]) > ub) {
            it=10;
            while (std::abs(g->r[i]) > ub && j<it) {
                j++;
                g->s[i] -= g->s[i] / 5;
                getAngRad(sol->XS, sol->YS, i, g->s[i], g->a[i], g->r[i]);
            }
        }

        std::array<double,2> startpoint,endpoint;
        std::array<double,2> mres = findMidAlphaBeta(sol->XS, sol->YS, g->dis, g->r[i], g->s[i], i,startpoint,endpoint); // Update center point
        assert(g->r[i]>0);
        g->m[i][0] = mres[0]; g->m[i][1] = mres[1]; //this should hopefully work since m is already initialized
        g->arc_startpoints_x[i] = startpoint[0];
        g->arc_startpoints_y[i] = startpoint[1];
        g->arc_endpoints_x[i] = endpoint[0];
        g->arc_endpoints_y[i] = endpoint[1];
     }
}
//prints out a point
void printpoint(const point P){std::cout<<P[0]<<","<<P[1]<<"\n";}

//since XS and YS have nan values in the back, I need a function that returns the last non-nan value
double lastnonnan(const std::array<double,10> XS)
{
    double lastNonNaN;
    for (size_t i = XS.size(); i > 0; --i) {
            if (!std::isnan(XS[i - 1])) {
                lastNonNaN = XS[i - 1];
                break;
            }
        }
    return lastNonNaN;
}

double discretize(SMAsol_struct* sol, const int k, const int n,const double distance_discrete) 
{
    //todo:Make this an optional algorithm parameter
    //const double distance_discrete = 0.2;
    struct_internalpath* g=&(sol->internalpath);
    std::vector<double> xx,yy;
    double totalpathlength=0;
    //calculate points on arc; reuse function from processpath here
    for(int i=0;i<4;i++)
    {
        //setup
        point startpoint = {g->arc_startpoints_x[i],g->arc_startpoints_y[i]};
        point endpoint = {g->arc_endpoints_x[i],g->arc_endpoints_y[i]};
        //Logic for appending the straight lines
        point lineVector;
        point runner;
        if(i==0)//first line: use XS(0) instead of endpoint; others: use endpoint from previous arc
        {
            lineVector={startpoint[0]-sol->XS[0],startpoint[1]-sol->YS[0]};
            runner={sol->XS[0],sol->YS[0]};}
        else
        {
            lineVector={startpoint[0]-g->arc_endpoints_x[i-1],startpoint[1]-g->arc_endpoints_y[i-1]};
            runner={g->arc_endpoints_x[i-1],g->arc_endpoints_y[i-1]};}
        //find vector of straight line
        double lineVectorLength = sqrt(lineVector[0]*lineVector[0] + lineVector[1]*lineVector[1]);
        totalpathlength+=lineVectorLength; //incrementally add up total path length
        lineVector[0]/=lineVectorLength;
        lineVector[1]/=lineVectorLength; //create unit vector
        //runner follows straight line, entries are added to xx,yy
        int n_linepoints = floor(lineVectorLength/distance_discrete);
        for(int j=0;j<n_linepoints;j++)
        {
            xx.push_back(runner[0]);
            yy.push_back(runner[1]);
            runner[0]+=lineVector[0]*distance_discrete;
            runner[1]+=lineVector[1]*distance_discrete;
        }

        //Logic for appending the arcs
        point centerpoint={g->m[i][0],g->m[i][1]};
        //!angles below 180° are assumed
        point vecCenterStart={centerpoint[0]-startpoint[0],centerpoint[1]-startpoint[1]};
        point vecCenterEnd={centerpoint[0]-endpoint[0],centerpoint[1]-endpoint[1]};
        double cross=vecCenterStart[0]*vecCenterEnd[1]-vecCenterStart[1]*vecCenterEnd[0];
        if(cross<0){cross=1;}else if(cross>0){cross=-1;}
        double arcLength;
        std::vector<point> arcpoints = calculatePointsOnArc_legacy(startpoint, centerpoint,endpoint, cross, distance_discrete,0,arcLength);
        totalpathlength+=arcLength;
        // write arc points into xx and yy
        for(int j=0;j<arcpoints.size();j++)
        {
            xx.push_back(arcpoints[j][0]);
            yy.push_back(arcpoints[j][1]);
        }
    }

    //last line from arc to endpoint
    point endpoint = {g->arc_endpoints_x[k-3],g->arc_endpoints_y[k-3]};
    point lineVector={lastnonnan(sol->XS)-endpoint[0],lastnonnan(sol->YS)-endpoint[1]};
    point runner = {endpoint[0],endpoint[1]};
    double lineVectorLength = sqrt(lineVector[0]*lineVector[0] + lineVector[1]*lineVector[1]);
    lineVector[0]/=lineVectorLength;
    lineVector[1]/=lineVectorLength; //create unit vector
    //runner follows straight line, entries are added to xx,yy
    int n_linepoints = floor(lineVectorLength/distance_discrete);
    for(int j=0;j<n_linepoints;j++)
    {
        xx.push_back(runner[0]);
        yy.push_back(runner[1]);
        runner[0]+=lineVector[0]*distance_discrete;
        runner[1]+=lineVector[1]*distance_discrete;
    }

    
    if (!(xx.back() == lastnonnan(sol->XS)) || !(yy.back() == lastnonnan(sol->YS)))
    {
        xx.push_back(lastnonnan(sol->XS));
        yy.push_back(lastnonnan(sol->YS));
    }

    sol->xx=xx;sol->yy=yy;
    return totalpathlength;
    }

//get Violation score for rectangles
//todo: While there are round edges in the map svg's corners, these are not implemented here yet
double getViolationRect(const model_struct* model,const SMAsol_struct* sol,const int k) {
    int amount_points=sol->xx.size();
    int i=0;
    double x,y,obs_left,obs_right,obs_up,obs_down,obs_left_tol,obs_right_tol,obs_up_tol,obs_down_tol;
    double Violation=0;
    const double add_collision=100.0/amount_points;//100.0 is important, otherwise we are dividing integers by integers, which returns wrong results
    bool collision_middle_side,collision_corner,collision_x,collision_y,coll_x_tol,coll_y_tol;
    while(i<amount_points) //Iterate over points generated in discretize
    {   
        x=sol->xx[i];
        y=sol->yy[i];
        for(int j=0;j<model->xwidth.size();j++) //Iterate over obstacles
        {
            if(abs(model->xobs_rect[j])<=10000000 && abs(model->yobs_rect[j])<=10000000)//there was a bug before where
            {
            //check if there is overlap in x and y direction without tolerance
            obs_left=model->xobs_rect[j]-0.5*model->xwidth[j];
            obs_right=model->xobs_rect[j]+0.5*model->xwidth[j];
            obs_up=model->yobs_rect[j]+0.5*model->ywidth[j];
            obs_down=model->yobs_rect[j]-0.5*model->ywidth[j];
            collision_x=(x-obs_left>=0) && (x-obs_right<=0);
            collision_y=(y-obs_down>=0) && (y-obs_up<=0);
            //check if point is on one of the sides (area covered by the tolerance, but not the round corners)
            obs_left_tol=obs_left-model->tolerance;
            obs_right_tol=obs_right+model->tolerance;
            obs_up_tol=obs_up+model->tolerance;
            obs_down_tol=obs_down-model->tolerance;
            coll_x_tol=(x-obs_left_tol>=0) && (x-obs_right_tol<=0);
            coll_y_tol=(y-obs_down_tol>=0) && (y-obs_up_tol<=0);
            collision_middle_side= ((coll_x_tol&&collision_y) || (coll_y_tol&&collision_x));
            //check if collision happened in one of the corners
            collision_corner=(coll_x_tol&&coll_y_tol) && !collision_middle_side;
            if(collision_middle_side)
            {
                Violation=Violation+add_collision;
            } 
            else if(collision_corner) //if collision happened in one of the corners
            {
                point rounded_centerp;//center point of the rounded corner the point collides with
                int upper_lower=sign(y-model->yobs_rect[j]);
                int right_left=sign(x-model->xobs_rect[j]);
                //determine center point of rectangle corner
                if(upper_lower==1 && right_left==1) //upper right corner
                {rounded_centerp={obs_right,obs_up};}
                else if(upper_lower==1&&right_left==-1)//upper left corner
                {rounded_centerp={obs_left,obs_up};}
                else if(upper_lower==-1&&right_left==1)//lower right corner
                {rounded_centerp={obs_right,obs_down};}
                else if(upper_lower==-1&&right_left==-1)//lower left corner
                {rounded_centerp={obs_left,obs_down};}
                double dist[2]={x-rounded_centerp[0], y-rounded_centerp[1]}; //distances of point to center point of corner
                if( sqrt(dist[0]*dist[0]+dist[1]*dist[1])<=model->tolerance ) //if distances small:point in circle => add Violation
                {Violation=Violation+add_collision;}
            }
            }
        }
        i++;
    }
    return Violation;
}

//getViolation for circles
double getViolationCircBoundary(const model_struct* model,const SMAsol_struct* sol) {
    int amount_points=sol->xx.size();
    int i=0;
    double Violation=0;
    double add_collision=100.0/amount_points;//100.0 is important, otherwise we are dividing integers by integers, which returns wrong results
    bool collision;
    while(i<amount_points)
    {
        for(int j=0;j<model->robs_circ.size();j++) //Iterate over obstacles
        {
            double dist[2] ={sol->xx[i]-model->xobs_circ[j],sol->yy[i]-model->yobs_circ[j]};
            collision= dist[0]*dist[0] + dist[1]*dist[1]<=(model->robs_circ[j]+model->tolerance)*(model->robs_circ[j]+model->tolerance);
            if(collision){Violation=Violation+add_collision;};
        }
        i++;

        //while this has nothing to do with circles, we can also check boundary violations here
        double ub = model->ub-model->tolerance;
        double lb = model->lb + model->tolerance;
        collision = sol->xx[i]>ub || sol->yy[i]<lb;
        if(collision){Violation=Violation+add_collision;};
    }
    return Violation;
}

//Print out current values for XS and YS
void printxsys(const std::array<double, 10> XS, const std::array<double, 10> YS){
  std::cout<<"XS, YS: "<<std::endl;
  for(int i=0;i<10;i++){
    std::cout<<XS[i]<<", "<<YS[i]<<std::endl;
  }
  std::cout<<" "<<std::endl;
}

// Function to check if a point (px, py) is inside a rotated rectangle
bool pointInRotatedRectangle(const double rectX, const double rectY, const double rectWidth, const double rectHeight, double rotationAngle, const double px, const double py) {
    // Convert the point coordinates to local coordinates of the rotated rectangle
    rotationAngle*=3.1415926/180;
    double localX = cos(-rotationAngle) * (px - rectX) - sin(-rotationAngle) * (py - rectY);
    double localY = sin(-rotationAngle) * (px - rectX) + cos(-rotationAngle) * (py - rectY);
    // Check if the point is inside the unrotated rectangle
    bool insideUnrotated = (localX >= 0 && localX <= rectWidth && localY >= 0 && localY <= rectHeight);

    return insideUnrotated;
}

//get Violation Score for Charging Station
double getViolationChargingStation(const SMAsol_struct* sol, const model_struct* model, const int k) {
    
    //setup
    double Violation=0;
    int amount_points=sol->xx.size();
    const double add_collision=100.0/amount_points;//100.0 is important, otherwise we are dividing integers by integers, which returns wrong results
    std::vector<turned_rectangle> blocked_area_down = model->blocked_area_down;
    std::vector<turned_rectangle> blocked_area_left = model->blocked_area_left;
    std::vector<turned_rectangle> blocked_area_up = model->blocked_area_up;
    std::vector<double> xx=sol->xx;
    std::vector<double> yy=sol->yy;
    bool collision;
    //for all the charging stations
    for(int j=0; j<blocked_area_down.size();j++)
    {
        turned_rectangle rect_down = blocked_area_down[j];
        turned_rectangle rect_left = blocked_area_left[j];
        turned_rectangle rect_up = blocked_area_up[j];
        //for all the points
        for(int i=0;i<amount_points;i++)
        {
            collision=pointInRotatedRectangle(rect_down.xpos,rect_down.ypos,rect_down.width,rect_down.height,rect_down.angle,xx[i],yy[i]);
            if(collision){Violation+=add_collision;}
            collision=pointInRotatedRectangle(rect_left.xpos,rect_left.ypos,rect_left.width,rect_left.height,rect_left.angle,xx[i],yy[i]);
            if(collision){Violation+=add_collision;}
            collision=pointInRotatedRectangle(rect_up.xpos,rect_up.ypos,rect_up.width,rect_up.height,rect_up.angle,xx[i],yy[i]);
            if(collision){Violation+=add_collision;}
        }
    }
  return Violation;
}

double getSmallRadVio(const SMAsol_struct* sol)
{
    const double addcollision = 2;
    double Violation=0;
    for(const auto &rad: sol->internalpath.r)
    {
        Violation += addcollision * (rad<1); //punish radii below 1m in radius
    }
    return Violation;
}
//#include "custom_math_operations.cpp"
#include <iostream>
#include <numeric>
#include "custom_datatypes.h"
#include <cassert>
#include <algorithm>
#include <cmath>
#include "slimemould_func_prototypes.h"
#include <cstdint>
//This function does some coordinate transformations.
//todo: Change how SA and EA are processed => The angles go clockwise, but should go counterclockwise => 90° and 270° are reversed
std::array<double,8> getPointstoXY(double x[6],const model_struct* model)
{
  double xx[4],yy[4];
 double model_SA=model->SA;
 double model_EA=model->EA;
 double model_xs=model->xs;
 double model_ys=model->ys;
 double model_xt=model->xt;
 double model_yt=model->yt;
  double RA_tmp,RB_tmp,a,absxk,alpha,beta,d,scale,t,v_rotatedA_idx_0,v_rotatedB_idx_0,y;

  alpha = 0.017453292519943295 * model_SA; 
  beta = 0.017453292519943295 * model_EA;
  //  Rotationsmatrix
  RA_tmp = std::sin(alpha);
  alpha = std::cos(alpha);
  RB_tmp = std::sin(beta);
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
  xx[0] = model_xs + v_rotatedA_idx_0;
  xx[3] = model_xt + v_rotatedB_idx_0;
  yy[0] = model_ys + a * d / y;
  xx[1] = x[1];
  yy[1] = x[2];
  xx[2] = x[3];
  yy[2] = x[4];
  yy[3] = model_yt + beta * alpha / RA_tmp;
  std::array<double,8> returnvalue = {xx[0],xx[1],xx[2],xx[3],yy[0],yy[1],yy[2],yy[3]};
  return returnvalue;
}

//pushes numbers in a fixed-size array to the left, deleting a number
//dependency for checkPoints
void push_to_left(std::array<double, 10> &arr, int index) {
    for (int i = index; i < arr.size() - 1; ++i) {
        arr[i] = arr[i + 1];
    }
}

void push_to_left(std::array<double, 9> &arr, int index) {
    for (int i = index; i < arr.size() - 1; ++i) {
        arr[i] = arr[i + 1];
    }
}

//checks if guide points are in order (if they are too close to one another, if they are in a straight line, if they lead nowhere)
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

//dependencies for initialization
void getAngRad(std::array<double, 10> XS, std::array<double, 10> YS, int i, double g_s, double& g_a, double& g_r) {
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

 //todo: deal with this edge case => g_r = 0 leads to problems further down in optimizationRad.
    // if (std::abs(g_a - 180) < 1e-4) {
    //     g_r = 0;
    // } else {
        // r bestimmen

        //std::cout<<"vectorBC:"<<vectorBC[0]<<", "<<vectorBC[1]<<std::endl;
        //std::cout<<"g_a:"<<g_a<<std::endl;
        g_r = g_s * std::tan(g_a / 2.0);
        g_a *= 180 / 3.14159;
        //add a few asserts to prevent g_r==0
        assert(g_s!=0);
        assert(std::tan(g_a)!=0);
        assert(g_r>=0);
   // }
}

//version without alpha and beta => called by initialization
void findMid(std::array<double, 10> XS, std::array<double, 10> YS, double dis[], double g_r, double g_s, int i, double* m0, double* m1) {
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
//dependency for findMidAlphaBeta; checks whether alpha or beta fulfills all required geometric conditions
bool checkalpha(double alphacandidate,double vecAB[2]){
    double len = sqrt(vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1]);
  bool alphacorrect = cos(alphacandidate) - vecAB[1]/len + sin(alphacandidate) - vecAB[0] / len < 1e-3; //if these conditions satisfied: alpha correct
  return alphacorrect;
}

//includes calculation of alpha and beta => called by optimizationRad
std::array<double,2> findMidAlphaBeta(std::array<double, 10> XS, std::array<double, 10> YS, std::vector<double> dis, double g_r, double g_s, int i, 
std::array<double,2> &startpoint, std::array<double,2> &endpoint) {
   // Mittelpunkt finden
   double middle0, middle1;
   
   double vecAB[2] = {XS[i + 1] - XS[i], YS[i + 1] - YS[i]};
    double vecBC[2] = {XS[i + 2] - XS[i + 1], YS[i + 2] - YS[i + 1]};
   if (g_r == 0) {
       //alpha = 0;
       //beta = 0;
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
  /* 
       double Vec[2];
       // get angle of alpha to x axis
       Vec[0] = -ax; Vec[1] = -ay;
       alpha = acos(Vec[0] / (std::sqrt(ax * ax + ay * ay) * 1))*sign(Vec[1]);// earlier version: * std::signbit(b);
       alpha = alpha /3.1415926 * 180;
   
       // get angle of beta to x axis
       Vec[0] = -bx; Vec[1] = -by;
       beta = acos(Vec[0] / (std::sqrt(bx * bx + by * by) * 1))*sign(Vec[1]);//earlier version: * std::signbit(d);
       beta = beta /3.1415926 *180;*/
   }
   //std::cout<<"arc number: "<<i<<", vecAB: "<<vecAB[0]<<","<<vecAB[1]<<", vecBC: "<<vecBC[0]<<", "<<vecBC[1]<<std::endl;
    //extensively tested; correct until here
    //calculate start and end point of arc
    double len = sqrt(vecAB[0]*vecAB[0]+vecAB[1]*vecAB[1]);
    vecAB[0]/=len;
    vecAB[1]/=len;
    //double startpoint[2];
    len=sqrt(vecBC[0]*vecBC[0]+vecBC[1]*vecBC[1]);
    startpoint[0]=XS[i+1] -g_s*vecAB[0];
    startpoint[1]=YS[i+1] -g_s*vecAB[1];
    vecBC[0]/=len;
    vecBC[1]/=len;
    //double endpoint[2];
    endpoint[0]=XS[i+1] +g_s*vecBC[0];
    endpoint[1]=YS[i+1] +g_s*vecBC[1];
    //std::cout<<"vecBC: "<<vecBC[0]<<", "<<vecBC[1]<<std::endl;
   
    //calculate alpha
    /*if(vecAB[0]>=0){ //positive in x axis
        if(vecAB[1]>0 || vecAB==0){alpha =acos(abs(vecAB[1])/dis[i])*180/3.14159; assert(alpha>=0 && alpha<=90);} //positive in x and y
        if(vecAB[1]<0){alpha = acos(abs(vecAB[1])/dis[i])*180/3.14159+90; assert(alpha>90 && alpha <=180);} //pos in x, neg in y
    }
    else{ //negative in x axis
        if(vecAB[1]>0|| vecAB==0){alpha = 270-acos(abs(vecAB[1]/dis[i]))*180/3.14159; assert(alpha>=270 && alpha<=360);} //positive in x and y
        if(vecAB[1]<0){alpha = 360-acos(abs(vecAB[1]/dis[i]))*180/3.14159; assert(alpha>180 && alpha <=270);} //pos in x, neg in y
    }
    //calculate beta
    if(vecBC[0]>=0){ //positive in x axis
        if(vecBC[1]>0 || vecBC==0){beta =acos(abs(vecBC[1])/dis[i+1])*180/3.14159; assert(beta>=0 && beta<=90);} //positive in x and y
        if(vecBC[1]<0){beta = acos(abs(vecBC[1])/dis[i+1])*180/3.14159+90; assert(beta>90 && beta <=180);} //pos in x, neg in y
    }
    else{ //negative in x axis
        if(vecBC[1]>0|| vecBC==0){beta = 270-acos(abs(vecBC[1]/dis[i+1]))*180/3.14159; assert(beta>=270 && beta<=360);} //positive in x and y
        if(vecBC[1]<0){beta = 360-acos(abs(vecBC[1]/dis[i+1]))*180/3.14159; assert(beta>180 && beta <=270);} //pos in x, neg in y
    }*/
  
  std::array<double,2> result = {middle0,middle1};
  return result;
}

//initialize arcs with radii
void initialization(struct_internalpath* g, int k, std::array<double, 10> XS, std::array<double, 10> YS)
{
  double g_s[k-2],g_a[k-2],g_r[k-2],alpha[k-2],beta[k-2],m[2][k-2],dis[k-1];
  //get distances between guide points
  double di;
  for (int i = 0; i < k - 1; ++i) {
    di = std::sqrt(std::pow(XS[i + 1] - XS[i], 2) + std::pow(YS[i + 1] - YS[i], 2));
    dis[i] = di; //distances between guide points
  }
  //for-loop for initialization
  for(int i=0;i<k-2;i++)
  {
    //if(i=k-2){g_s[i] = dis[i+1]*0.9;} //last arc should be as long as dis between guide point and end point
    g_s[i] = std::min(dis[i],dis[i+1])*0.25; //in MATLAB Code: 0.25 is called 'start'; amount of space on lines used for arcs
    getAngRad(XS,YS,i,g_s[i],g_a[i],g_r[i]); //get angles and radii for a given g_s
    findMid(XS,YS,dis,g_r[i],g_s[i],i,&m[0][i],&m[1][i]); //find middle point
    g->a.push_back(g_a[i]); //write angles, sizes, center points, radii and distances into g
    g->s.push_back(g_s[i]);
    g->r.push_back(g_r[i]);
    g->m.push_back({m[0][i],m[1][i]});
    assert(g_r[i]>0);
  }
  for(int i=0;i<k-1;i++){g->dis.push_back(dis[i]);}
}

//function for plotting arcs=> calculates points along arcs
//returns xunit, yunit, a
void circleRad(double x, double y, double r, double alpha, double beta,
               std::vector<double>& xunit, std::vector<double>& yunit, double& a) {
    //x,y: center point of current arc
    //r: radius of current arc
    //alpha, beta: angles at start, end
    //xunit, yunit: xCircle, yCircle => points on arc
    //a: angles of these points (I think)
    if (alpha < 0) alpha += 360;
    if (beta < 0) beta += 360;

    alpha = alpha * 3.14159 / 180.0; // degtorad
    beta = beta * 3.14159 / 180.0;

    std::vector<double> th1, th2, th; //th1 and th2 go around the arc in two different directions. Out of the two, th is the shorter one

    if (alpha < beta) {
        for (double th = alpha; th <= beta; th += 3.14159 / 100) //fill up th1
            th1.push_back(th);
    } else {
        for (double th = beta; th <= alpha; th += 3.14159 / 100)
            th1.push_back(th);
        std::reverse(th1.begin(), th1.end()); //flip array
    }

    if (th1.size() == 1) //arc smaller than pi/100
        th1 = { alpha, beta };

    if (alpha > beta) { //fill up th2
        for (double th = alpha; th <= 2 * 3.14159; th += 3.14159 / 100)
            th2.push_back(th);
        for (double th = 3.14159 / 100; th <= beta; th += 3.14159 / 100)
            th2.push_back(th);
    } else {
        for (double th = beta; th <= 2 * 3.14159; th += 3.14159 / 100)
            th2.push_back(th);
        for (double th = 3.14159 / 100; th <= alpha; th += 3.14159 / 100)
            th2.push_back(th);
        std::reverse(th2.begin(), th2.end());
    }

    if (th2.size() == 1) //if angle of arc <pi/100
        th2 = { alpha, beta };

    th = (th1.size() < th2.size()) ? th2 : th1;

    for (double th_val : th) {
        xunit.push_back(r * cos(th_val) + x);
        yunit.push_back(r * sin(th_val) + y);
    }
    assert(th.size()>1); //make sure that th is filled with values to avoid segfault
    a = std::abs(th.back() - th.front()) * 180.0 / 3.14159; // radtodeg conversion

    if (a > 180)
        a = 360 - a;
}

void optimizationRad(struct_internalpath* g, int k, std::array<double, 10> XS,
 std::array<double, 10> YS, std::vector<std::vector<double>>* xCirc, std::vector<std::vector<double>>* yCirc) 
{
    std::vector<std::vector<double>> xCircle(k - 2); // Using vector instead of cell arrays
    std::vector<std::vector<double>> yCircle(k - 2);
    std::vector<double> degC(k - 2, 0.0);
    std::vector<double> alpha(k-2),beta(k-2);
    double it, temp, lb, ub; int j; double os;
    lb=4;ub=9;

    for (int i = 0; i < k - 2; i++) {
        
        if (std::abs(g->r[i]) < lb) {
            j = 1;
            os = g->s[i];
            lb=4;ub=9; //lower, upper bounds
            
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
                g->s[i] += os / 10;
                g->s[i] = std::min(g->s[i], temp);
                getAngRad(XS, YS, i, g->s[i], g->a[i], g->r[i]);
            }
        } else if (std::abs(g->r[i]) > ub) {
            it=10;
            while (std::abs(g->r[i]) > ub && j<it) {
                j++;
                g->s[i] -= g->s[i] / 5;
                getAngRad(XS, YS, i, g->s[i], g->a[i], g->r[i]);
            }
        }


        double alpha,beta;
        std::array<double,2> startpoint,endpoint;
        std::array<double,2> mres = findMidAlphaBeta(XS, YS, g->dis, g->r[i], g->s[i], i,startpoint,endpoint); // Update center point
        //g->alpha.push_back(alpha);
        //g->beta.push_back(beta);
        assert(g->r[i]>0);
        //assert(!std::isnan(g->alpha[i]) && !std::isnan(g->beta[i]));
        g->m[i][0] = mres[0]; g->m[i][1] = mres[1]; //this should hopefully work since m is already initialized
        xCircle[i].push_back(startpoint[0]);
        yCircle[i].push_back(startpoint[1]);
        xCircle[i].push_back(endpoint[0]);
        yCircle[i].push_back(endpoint[1]);
     }

    //for (int i = 0; i < k - 2; i++) {
    //    circleRad(g->m[i][0], g->m[i][1], g->r[i], g->alpha[i], g->beta[i], xCircle[i], yCircle[i], degC[i]);
    //}

   // Update the fields of struct g
   g->xCircle = xCircle;
   g->yCircle = yCircle;
   g->degC = degC;
   *xCirc = xCircle;
   *yCirc = yCircle;
}
//prints out a point
void printpoint(point P){std::cout<<P[0]<<","<<P[1]<<"\n";}

//since XS and YS have nan values in the back, I need a function that returns the last non-nan value
double lastnonnan(std::array<double,10> XS)
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

double discretize(struct_internalpath* g, SMAsol_struct* sol, int k, int n, std::array<double, 10> XS, std::array<double, 10> YS) 
{
    //todo:Make this an optional algorithm parameter
    double distance_discrete = 0.1;
    std::vector<double> xx,yy;
    double totalpathlength=0;
    //calculate points on arc; reuse function from processpath here
    for(int i=0;i<4;i++)
    {
        //setup
        point startpoint = {g->xCircle[i][0],g->yCircle[i][0]};
        point endpoint = {g->xCircle[i].back(),g->yCircle[i].back()};
        //Logic for appending the straight lines
        point lineVector;
        point runner;
        if(i==0)//first line: use XS(0) instead of endpoint; others: use endpoint from previous arc
        {
            lineVector={startpoint[0]-XS[0],startpoint[1]-YS[0]};
            runner={XS[0],YS[0]};}
        else
        {
            lineVector={startpoint[0]-g->xCircle[i-1].back(),startpoint[1]-g->yCircle[i-1].back()};
            runner={g->xCircle[i-1].back(),g->yCircle[i-1].back()};}
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
        //std::cout<<"start, end, center\n";
        //printpoint(startpoint);
        //printpoint(endpoint);
        //printpoint(centerpoint);
        //Logic for determining whether arc is clockwise or counterclockwise;
        //!angles below 180° are assumed
        point vecCenterStart={centerpoint[0]-startpoint[0],centerpoint[1]-startpoint[1]};
        point vecCenterEnd={centerpoint[0]-endpoint[0],centerpoint[1]-endpoint[1]};
        double cross=vecCenterStart[0]*vecCenterEnd[1]-vecCenterStart[1]*vecCenterEnd[0];
        if(cross<0){cross=1;}else if(cross>0){cross=-1;}
        double arcLength;
        std::vector<point> arcpoints = calculatePointsOnArc(startpoint, centerpoint,endpoint, cross, distance_discrete,0,arcLength);
        totalpathlength+=arcLength;
        // write arc points into xx and yy
        for(int j=0;j<arcpoints.size();j++)
        {
            xx.push_back(arcpoints[j][0]);
            yy.push_back(arcpoints[j][1]);
        }
    }

    //last line from arc to endpoint
    point endpoint = {g->xCircle[k-3].back(),g->yCircle[k-3].back()};
    point lineVector={lastnonnan(XS)-endpoint[0],lastnonnan(YS)-endpoint[1]};
    point runner = {endpoint[0],endpoint[1]};
    double lineVectorLength = sqrt(lineVector[0]*lineVector[0] + lineVector[1]*lineVector[1]);
    lineVector[0]/=lineVectorLength;
    lineVector[1]/=lineVectorLength; //create unit vector
    //runner follows straight line, entries are added to xx,yy
    //std::cout<<"last line: line vector\n";
    //printpoint(lineVector);
    int n_linepoints = floor(lineVectorLength/distance_discrete);
    for(int j=0;j<n_linepoints;j++)
    {
        xx.push_back(runner[0]);
        yy.push_back(runner[1]);
        runner[0]+=lineVector[0]*distance_discrete;
        runner[1]+=lineVector[1]*distance_discrete;
    }

    
    if (!(xx.back() == lastnonnan(XS)) || !(yy.back() == lastnonnan(YS)))
    {
        xx.push_back(lastnonnan(XS));
        yy.push_back(lastnonnan(YS));
    }

    //for debugging: print out points
    //std::cout<<"xx, yy:\n";
    //for(int i=0; i<xx.size();i++)
    //{
    //    std::cout<<"("<<xx[i]<<","<<yy[i]<<")\n";
    //}
    //todo: Handle exceptions or issue warnings if needed
    //if (XS.size() == 2) {
    //    std::array<std::vector<double>,2> temp;
    //    temp = splitLin({XS[0],XS[1]}, {XS.back(),YS.back()}, n);
    //    sol->xx=temp[0];sol->yy=temp[1];
    //}
    sol->xx=xx;sol->yy=yy;
    return totalpathlength;
    }

//get Violation score for rectangles
//todo: While there are round edges in the map svg's corners, these are not implemented here yet
double getViolationRect(struct_internalpath* g, std::vector<double> xobs, std::vector<double> yobs, std::vector<double> xwidth, std::vector<double> ywidth, int nobs, SMAsol_struct* sol, int k, int n,double tolerance) {
    std::for_each(xwidth.begin(), xwidth.end(), [](double& val) { val *= 0.5; }); // modify each value in the two vectors
    std::for_each(ywidth.begin(), ywidth.end(), [](double& val) { val *= 0.5; });
    int amount_points=sol->xx.size();
    int i=0;
    double x,y,obs_left,obs_right,obs_up,obs_down,obs_left_tol,obs_right_tol,obs_up_tol,obs_down_tol;
    double Violation=0;
    double add_collision=100.0/amount_points;//100.0 is important, otherwise we are dividing integers by integers, which returns wrong results
    bool collision_middle_side,collision_corner,collision_x,collision_y,coll_x_tol,coll_y_tol;
    while(i<amount_points) //Iterate over points generated in discretize
    {   
        //std::cout<<"i="<<i<<",";
        x=sol->xx[i];
        y=sol->yy[i];
        for(int j=0;j<nobs;j++) //Iterate over obstacles
        {
            //std::cout<<"j="<<j<<", xobs: "<<xobs[j]<<"yobs: "<<yobs[j]<<"\n";
            if(abs(xobs[j])<=10000000 && abs(yobs[j])<=10000000)//there was a bug before where
            {
            //check if there is overlap in x and y direction without tolerance
            obs_left=xobs[j]-xwidth[j];
            obs_right=xobs[j]+xwidth[j];
            obs_up=yobs[j]+ywidth[j];
            obs_down=yobs[j]-ywidth[j];
            collision_x=(x-obs_left>=0) && (x-obs_right<=0);
            collision_y=(y-obs_down>=0) && (y-obs_up<=0);
            //check if point is on one of the sides (area covered by the tolerance, but not the round corners)
            obs_left_tol=obs_left-tolerance;
            obs_right_tol=obs_right+tolerance;
            obs_up_tol=obs_up+tolerance;
            obs_down_tol=obs_down-tolerance;
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
                //std::cout<<"Collision in the corner\n";
                point rounded_centerp;//center point of the rounded corner the point collides with
                //std::cout<<"calculating sign value of "<<y-yobs[j]<<"\n";
                int upper_lower=sign(y-yobs[j]);
                int right_left=sign(x-xobs[j]);
                //std::cout<<"upper_lower: "<<upper_lower<<";right_left: "<<right_left<<"\n";
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
                if( sqrt(dist[0]*dist[0]+dist[1]*dist[1])<=tolerance ) //if distances small:point in circle => add Violation
                {Violation=Violation+add_collision;}
            }
            }
        }
        i++;
    }
    return Violation;
}

//getViolation for circles
double getViolationCirc(struct_internalpath* g, std::vector<double> xobs, std::vector<double> yobs, std::vector<double> robs, int nobs, SMAsol_struct* sol, double tolerance) {
    int amount_points=sol->xx.size();
    int i=0;
    double Violation=0;
    double add_collision=100.0/amount_points;//100.0 is important, otherwise we are dividing integers by integers, which returns wrong results
    bool collision;
    while(i<amount_points)
    {
        for(int j=0;j<nobs;j++) //Iterate over obstacles
        {
            double dist[2] ={sol->xx[i]-xobs[j],sol->yy[i]-yobs[j]};
            collision= dist[0]*dist[0] + dist[1]*dist[1]<=(robs[j]+tolerance)*(robs[j]+tolerance);
            if(collision){Violation=Violation+add_collision;};
        }
        i++;
    }
    return Violation;
}

//Print out current values for XS and YS
void printxsys(std::array<double, 10> XS, std::array<double, 10> YS){
  std::cout<<"XS, YS: "<<std::endl;
  for(int i=0;i<10;i++){
    std::cout<<XS[i]<<", "<<YS[i]<<std::endl;
  }
  std::cout<<" "<<std::endl;
}

// Function to check if a point (px, py) is inside a rotated rectangle
bool pointInRotatedRectangle(double rectX, double rectY, double rectWidth, double rectHeight, double rotationAngle, double px, double py) {
    // Convert the point coordinates to local coordinates of the rotated rectangle
    rotationAngle*=3.1415926/180;
    double localX = cos(-rotationAngle) * (px - rectX) - sin(-rotationAngle) * (py - rectY);
    double localY = sin(-rotationAngle) * (px - rectX) + cos(-rotationAngle) * (py - rectY);
    // Check if the point is inside the unrotated rectangle
    bool insideUnrotated = (localX >= 0 && localX <= rectWidth && localY >= 0 && localY <= rectHeight);

    return insideUnrotated;
}

//get Violation Score for Charging Station
double getViolationChargingStation(struct_internalpath* g, SMAsol_struct* sol, const model_struct* model, int k, int n) {
    
    //setup
    double Violation=0;
    int amount_points=sol->xx.size();
    double add_collision=100.0/amount_points;//100.0 is important, otherwise we are dividing integers by integers, which returns wrong results
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
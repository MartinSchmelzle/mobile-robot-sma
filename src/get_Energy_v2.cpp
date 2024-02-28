//todo: continue here
#include<array>
#include <iostream>
#include <iomanip>
#include<cmath>
#include<cassert>
#include "../custom_datatypes.h"
#include "../slimemould_func_prototypes.h"

//This function finds the index from which an array should be populated further
int find_i_start(const float array[10000])
{
    int i_start;
    for (i_start = 0; i_start < 10000; ++i_start) {
        //std::cout<<array[i_start]<<"\n";
        if (array[i_start] == 0) {
            break;  // Found the first zero entry
        }
    }
    i_start-=1;
    if(i_start<0){i_start=0;} //this happens at first iteration of first line
    return i_start;
}

//print out result of kinematics calculation for testing / debugging
//If you're planning to use this, I recommend calling the script in an external terminal instead of an IDE-internal one => easier to look at
void printKin(const kinematics_time& kinematics) {
    // Print header
    std::cout << std::setw(10) << "Index" << std::setw(15) << "Time" << std::setw(15) << "Acceleration"
              << std::setw(15) << "Velocity" << std::setw(15) << "x Position"
              << std::setw(15) << "y Position" << std::setw(15) << "distance"
              << std::setw(15) << "omega" <<std::setw(15) << "dw/dt" << std::endl;
    // Print array values
    int i_end = find_i_start(kinematics.a_array) +3; //add +3 to also get a few empty values in
    for (int i = 0; i < i_end; ++i) {
            std::cout << std::setw(10) << i << std::setw(15) << kinematics.time_array[i]
                      << std::setw(15) << kinematics.a_array[i] << std::setw(15) << kinematics.v_array[i]
                      << std::setw(15) << kinematics.x_array[i] << std::setw(15) << kinematics.y_array[i]
                      << std::setw(15) << kinematics.s_array[i] << std::setw(15) << kinematics.w_array[i]
                      << std::setw(15) << kinematics.dwdt_array[i] << std::endl;
    }
}

//This function calculates the amount of distance from the last arc needed for deceleration
//If the last line segment is very short, it cannot cover the entire deceleration process and we also need to break during the last arc
//! assumption: last arc + last line segment are always enough for the entire deceleration process
double get_decel_dist(kinematics_time* kin, double length_last)
{
    double s_break = std::abs(0.5* kin->v_desired*kin->v_desired / kin->a_break);
    //std::cout<<"distance needed for deceleration: "<<s_break<<", length of last segment: "<<length_last<<"\n"; 
    double arc_break_length;
    if(s_break<=length_last) {arc_break_length=0;} //in this case: everything ok, the last line is enough
    else{arc_break_length=s_break - length_last;} //in this case, we need to break during the arc
    return arc_break_length;
    //double v = lineLength * kin->a_break;
}

//This function takes in a straight line denoted by a line Vector and calculates position, velocity and acceleration as a function of time
void get_kinematics_line(kinematics_time* kin, std::array<double,2> lineVector, const double lineLength, const bool breakatend)
{
    //step 1: read out the required data from h, set everything up
    const int time_size=10000;

    //step 1.5: Since we have multiple line segments in sequence, we might not need to start at entry 0, but instead another entry => find that entry
    const int i_start = find_i_start(kin->a_array);
    //std::cout<<"i_start: "<<i_start<<"\n";
    // Now i_start holds the index where you need to continue populating the array (=continue the journey)

    //step 2: create a vector representing time => should be unnecessarily long
    float value=kin->time_array[i_start];
    for (int i = i_start; i < time_size; i++) {
        kin->time_array[i] = value;
        value += kin->delta_t;
    }
    //std::cout<<"Continue populating array at index i="<<i_start<<" at s="<<kin->s_array[i_start]<<" and t="<<kin->time_array[i_start]<<"\n";

    //step 3: calculate a, v and r for straight lines => you need all of them at the same time
    const double s_start = kin->s_array[i_start];
    bool break_cond;
    int i=i_start;
    while(i<time_size){
        double a_current, v_next, x_next, y_next;
        if(breakatend && lineLength - kin->v_desired*kin->v_desired/ (-2 * kin->a_break) - (kin->s_array[i]-s_start)<0.1) // if the robot is far enough that it should break
        {break_cond=1;}
        else{break_cond= 0;}//This sets whether the robot is supposed to decelerate in the current iteration
        //step 3.1: compute value of acceleration for current step
        //todo: Convert this into a proper differential equation with variable mass etc. (long-term)
        if(break_cond)                  {a_current=kin->a_break;} //if robot is decelerating
        else if(kin->v_array[i]>kin->v_max)       {a_current=kin->a_break;} //if robot is going way too fast
        else if(kin->v_array[i]<kin->v_max*0.9)   {a_current=kin->a_acc;} //if robot is accelerating
        else                            {a_current=0.2*(kin->v_array[i]-kin->v_desired)/kin->delta_t;} // If robot is in the right window (90% to 100% of v_max), let it slowly go to the desired value

        //step 3.2: write acceleration into array
        kin->a_array[i]=a_current;

        //step 3.3: calculate v, r of next step
        kin->v_array[i+1] = kin->v_array[i] + a_current * kin->delta_t;
        kin->x_array[i+1] = (kin->v_array[i] * kin->delta_t + 0.5* a_current*kin->delta_t*kin->delta_t)*lineVector[0] + kin->x_array[i];
        kin->y_array[i+1] = (kin->v_array[i] * kin->delta_t + 0.5* a_current*kin->delta_t*kin->delta_t)*lineVector[1] + kin->y_array[i];
        kin->s_array[i+1] = kin->v_array[i] * kin->delta_t + 0.5* a_current*kin->delta_t*kin->delta_t + kin->s_array[i];

        //std::cout<<"i: "<<i<<", time: "<<kin->time_array[i]<<", a: "<<a_current<<", v: "<<kin->v_array[i]<<", x: "<<kin->x_array[i]<<", y: "<<kin->y_array[i]<<"\n";

        //step 3.4: check if the robot finished its journey
        if(break_cond && kin->v_array[i+1]<=0){break;}
        else if(kin->s_array[i] - s_start>=lineLength || (a_current==0 &&kin->v_array[i]==0 && kin->x_array[i]==0)){break;}

        i++;
    }
    //step 4: At the end, the result is between the i-1th and ith entry. Find out the time at which the robot reaches the line's end. This is not true if breakatend=1
    if (!breakatend){
        //use quadratic solution formula for that
        float time_end= (-kin->v_array[i-1]+std::sqrt(kin->v_array[i-1]*kin->v_array[i-1] -2*kin->a_array[i-1]*(kin->s_array[i-1]-s_start-lineLength))) / kin->a_array[i-1];
        time_end += kin->time_array[i-1];
        //set acceleration for end of line with linear interpolation
        float a_current = (time_end-kin->time_array[i-1]) / (kin->time_array[i]-kin->time_array[i-1]) * (kin->a_array[i] - kin->a_array[i-1]) + kin->a_array[i-1];
        kin->time_array[i]=time_end;
        float delta_t_end = time_end - kin->time_array[i-1];
        kin->a_array[i] = a_current;
        //update velocity and position so the robot would end up at the required end position
        kin->v_array[i] = kin->v_array[i-1] + kin->a_array[i-1] * delta_t_end;
        kin->x_array[i] = (kin->v_array[i-1] * delta_t_end + 0.5* kin->a_array[i-1]*delta_t_end*delta_t_end)*lineVector[0] + kin->x_array[i-1];
        kin->y_array[i] = (kin->v_array[i-1] * delta_t_end + 0.5* kin->a_array[i-1]*delta_t_end*delta_t_end)*lineVector[1] + kin->y_array[i-1];
        //std::cout<<"i: "<<i<<", time: "<<kin->time_array[i]<<", a: "<<kin->a_array[i]<<", v: "<<kin->v_array[i]<<", x: "<<kin->x_array[i]<<", y: "<<kin->y_array[i]<<"\n";
    }
}

//this function calculates time-dependent arc kinematics
void get_kinematics_arc(kinematics_time* kin, const double arc_break_length, const double centerPoint[2], const double endpoint[2], const double radius)
{
    //step 10: get i_start to know where to continue
    const int i_start = find_i_start(kin->a_array);

    //step 20: From x,y: derive start point.
    const double startpoint[2] = {kin->x_array[i_start], kin->y_array[i_start]};

    //step 30: With start point and end point, compute arc length  => same as findMid.
    //todo: Make an arc class so that you wouldn't have to calculate this every time
    const double vectorStart[2] = {startpoint[0] - centerPoint[0], startpoint[1] - centerPoint[1]};
    //std::cout<<"vectorStart: ("<<vectorStart[0]<<","<<vectorStart[1]<<")\n";
    const double vectorEnd[2] = {endpoint[0] - centerPoint[0], endpoint[1] - centerPoint[1]};
    //std::cout<<"vectorEnd: ("<<vectorEnd[0]<<","<<vectorEnd[1]<<")\n";

    const double dotProduct = vectorStart[0] * vectorEnd[0] + vectorStart[1] * vectorEnd[1];

    const double magnitudeStart = std::sqrt(vectorStart[0] * vectorStart[0] + vectorStart[1] * vectorStart[1]);
    const double magnitudeEnd = std::sqrt(vectorEnd[0] * vectorEnd[0] + vectorEnd[1] * vectorEnd[1]);

    const double angle = std::acos(dotProduct / (magnitudeStart * magnitudeEnd));
    const double arcLength = radius * angle;
    //std::cout<<"arcLength: "<<arcLength<<"\n";

    //step 35: Compute the arc direction (clockwise or counterclockwise)
    const double crossProduct = vectorStart[0] * vectorEnd[1] - vectorStart[1] * vectorEnd[0];
    const bool isClockwise = (crossProduct < 0); // true if clockwise, false if counterclockwise/ You can determine this based on the sign of the cross product of the vectors

    //step 40: With arc length, call get_kinematics_line to compute a, v, and s
    //std::cout<<"arc_break_length: "<<arc_break_length<<"\n";
    std::array<double,2> vec={arcLength-arc_break_length,0};
    get_kinematics_line(kin,vec,vec[0],0);
    if(arc_break_length>0){ //This section gets called if we need to berak during the (last) arc
        //std::cout<<"breaking during arc.\n";
        vec[0]=arc_break_length;
        vec[1]=0;
        get_kinematics_line(kin,vec,vec[0],1);
    }

    //step 50: with s: calculate correct x and y values
    const int i_end=find_i_start(kin->a_array);
    assert(i_end>=i_start);
    const double s_start = kin->s_array[i_start]; //s value at start of arc
    const double angle_start = std::atan2(vectorStart[1], vectorStart[0]); //angle towards (1,0) at start of arc
    //std::cout<<"angle_start: "<<angle_start<<"\n";

    for(int i=i_start+1;i<=i_end;i++)
    {
        const double delta_angle = (kin->s_array[i]-s_start) / radius;
        const double vector_current_point[2] = {radius * std::cos(angle_start + delta_angle),radius * std::sin(angle_start + delta_angle)};
        kin->x_array[i]=vector_current_point[0]+centerPoint[0];
        kin->y_array[i]=vector_current_point[1]+centerPoint[1];
        //step 60: with v, calculate w
        kin->w_array[i]=kin->v_array[i]/radius; //calculate omega
        //step 70: with tangential acceleration, compute dw/dt
        kin->dwdt_array[i]=kin->a_array[i]/radius; 
        //std::cout<<"i: "<<i<<", time: "<<kin->time_array[i]<<", w: "<<kin->w_array[i]<<", dwdt: "<<kin->dwdt_array[i]<<", x: "<<kin->x_array[i]<<", y: "<<kin->y_array[i]<<"\n";
    }
}

double getEnergy_v2(SMAsol_struct* sol, kinematics_time* kin)
{
    //std::cout<<"getEnergy started\n";
    double lineVector[2];

    //calculate when we have to start decelerating
    double arc_decel_dist = get_decel_dist(kin, euclid_dist({sol->internalpath.arc_endpoints_x[3],sol->internalpath.arc_endpoints_y[3]},{sol->XS[5],sol->YS[5]}));
    //std::cout<<"arc_decel_dist: "<<arc_decel_dist<<"\n";
    for(int i=0;i<4;i++)
    {
        if(i==0)
        {
        kin->x_array[0] = sol->XS[0]; //set start coordinates bcs otherwise the function assumes that we start at (0,0)
        kin->y_array[0] = sol->YS[0];
        } //make sure it's the vector with the actual length, not a unit vector!
        else
        {
        //lineVector[0] =sol->internalpath.arc_startpoints_x[i] - sol->internalpath.arc_endpoints_x[i-1];
        //lineVector[1]=sol->internalpath.arc_startpoints_y[i] - sol->internalpath.arc_endpoints_y[i-1];
        }
        get_kinematics_line(kin, sol->internalpath.lineVector[i], sol->internalpath.lineLength[i], 0);
        //std::cout<<"get_kinematics_line done for i="<<i<<"\n";
        double center[2]={sol->internalpath.m[i][0], sol->internalpath.m[i][1]};
        double end[2]={sol->internalpath.arc_endpoints_x[i],sol->internalpath.arc_endpoints_y[i]};
        //std::cout<<"get_kinematics_arc called for i="<<i<<"\n";
        get_kinematics_arc(kin, (i==3)*arc_decel_dist, center, end, sol->internalpath.r[i]); //in last arc: the robot may have to decelerate
        //std::cout<<"get_kinematics_arc done for i="<<i<<"\n";
    }
    //process last line segment
    lineVector[0] = -sol->internalpath.arc_endpoints_x[3] + sol->XS[5];
    lineVector[1] = -sol->internalpath.arc_endpoints_y[3] + sol->YS[5];
    //std::cout<<"lineVector of last segment: ["<<lineVector[0]<<","<<lineVector[1]<<"]\n";
    get_kinematics_line(kin, sol->internalpath.lineVector[4], sol->internalpath.lineLength[4], 1);
    //std::cout<<"get_kinematics_line done for i=4\n";

    //based on the kinematics, I can now compute the energy demand as a function of time =>HS's work, equation 4.12
    double E_need = 0;
    double mass=kin->mass;
    int i_end=find_i_start(kin->a_array);
    double E_dummy = i_end*kin->P_idle; //P_idle is constant => does not need to be added in each iteration, this way it's faster
    const std::array<double,8> p = {0.244,2.77,1.94,-5.26,0.0934,0.0911,0.00129,-0.00626}; //values for combined drive
    for(int i=0;i<i_end;i++)
    {
        if(std::isnan(kin->v_array[i])){kin->v_array[i]=0;}
        if(std::isnan(kin->a_array[i])){kin->a_array[i]=0;}
        double P_komb =p[0]*kin->v_array[i] * ( mass * (p[1]+p[2] * kin->a_array[i]) + p[3] * std::pow(kin->v_array[i],2))+
                    p[4]*std::abs(kin->w_array[i]) * (mass * (p[5] + p[6] * kin->dwdt_array[i] * sign(kin->w_array[i])) + p[7] *std::pow(kin->w_array[i],2));
        if (std::isnan(P_komb))
        {
            std::cout<<"Data point "<<i<<": P_komb is nan. v: "<<kin->v_array[i]<<", a: "<<kin->a_array[i]<<", w: "<<kin->w_array[i]<<", dwdt: "<<kin->dwdt_array[i]<<"\n";
        }
        E_need += P_komb;
    }
    return std::abs(E_need * kin->delta_t /1000);
}
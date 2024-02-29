#include<iostream>
#include<cassert>
#include <array>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "../slimemould_func_prototypes.h"
#include "../custom_datatypes.h" 

using point = std::array<double, 2>;
using arr4 = std::array<double,4>;

//This function takes in an array of length 4 and returns it, except all entries are unsigned and multiplied by scale
inline std::array<unsigned,4> convert_unsigned_scale(std::array<double,4> input, int scale)
{
  std::array<unsigned,4> res;
  for(int i=0;i<4;i++)
  {
    res[i] = static_cast<unsigned>(input[i]*scale);
  }
  return res;
}

void printPathStruct(const SMAsol_struct* sol) {
    std::cout << "Start Point: (" << sol->XS[0] << ", " << sol->YS[1] << ")\n";
    std::cout << "End Point: (" << sol->XS[5] << ", " << sol->YS[5] << ")\n";
    
    for (int i = 0; i < 4; ++i) {
        std::cout << "Arc " << i+1 << ":\n";
        std::cout << "    Start X: " << sol->internalpath.arc_startpoints_x[i] << "\n";
        std::cout << "    Start Y: " << sol->internalpath.arc_startpoints_y[i] << "\n";
        std::cout << "    End X: " << sol->internalpath.arc_endpoints_x[i] << "\n";
        std::cout << "    End Y: " << sol->internalpath.arc_endpoints_y[i] << "\n";
        std::cout << "    Center X: " << sol->internalpath.m[i][0] << "\n";
        std::cout << "    Center Y: " << sol->internalpath.m[i][1] << "\n";
        std::cout << "    Radius: " << sol->internalpath.r[i] << "\n";
        std::cout << "    Counterclockwise: " << sol->internalpath.arc_direction[i] << "\n\n";
    }
}

//This function takes in a path struct and returns a string array containing the svg lines.
std::array <std::string,4> formatarc_svg(SMAsol_struct* sol, unsigned scale)
{
    std::array <std::string,4> lines;
    //Modify arc data (scale + conversion to unsigned)
    #pragma GCC diagnostic ignored "-Wnarrowing"//Ignore warning about datatype conversion
    std::array<unsigned,4> arc_startpoints_x_u=convert_unsigned_scale(sol->internalpath.arc_startpoints_x,scale);
    std::array<unsigned,4> arc_startpoints_y_u=convert_unsigned_scale(sol->internalpath.arc_startpoints_y,scale);
    std::array<unsigned,4> arc_endpoints_x_u=convert_unsigned_scale(sol->internalpath.arc_endpoints_x,scale);
    std::array<unsigned,4> arc_endpoints_y_u=convert_unsigned_scale(sol->internalpath.arc_endpoints_y,scale);
    #pragma GCC diagnostic warning "-Wnarrowing"

    for(int i=0;i<4;i++)
    {
        //calculate radius
        double arc_radius = sol->internalpath.r[i]*scale;
        //convert flag for clockwise/counterclockwise => in svg, it's 0/1, in the alg it's 1/0/-1
        bool clock=sol->internalpath.arc_direction[i]>=0;

        //write line for arc
        lines[i] = "<path d=\"M " + std::to_string(arc_startpoints_x_u[i]) + "," + std::to_string(arc_startpoints_y_u[i]) +
                          " A " + std::to_string(arc_radius) + "," + std::to_string(arc_radius) + " 0 0 " + std::to_string(clock) + " " +
                          std::to_string(arc_endpoints_x_u[i]) + "," + std::to_string(arc_endpoints_y_u[i]) +
                          "\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n";;
    }
    return lines;
}


//This function draws a given path of type path struct into the svg file "map.svg".
int drawpathintosvg(SMAsol_struct* sol,unsigned scale)
{
    //INPUT FILE
    std::ifstream inputFile("map.svg");

    if (inputFile.is_open()) {
        // Read the content of the file
        std::vector<std::string> content;
        std::string line;
        while (std::getline(inputFile, line)) {
            if (line.find("</svg>") != std::string::npos) {
                break; // Stop reading when </svg> is found so that the line ending the svg would not be read
            }
            content.push_back(line);
        }

        // Close the input file
        inputFile.close();

        //Modify line data (scale + conversion to integer)
        #pragma GCC diagnostic ignored "-Wnarrowing" //ignore warning about loss of data
        std::array<unsigned, 2> startPoint_u={static_cast<int>(sol->XS[0]*scale), static_cast<int>(sol->YS[0]*scale)}; //scale to internal representation, round to integer
        std::array<unsigned, 2> endPoint_u={static_cast<int>(sol->XS[5]*scale), static_cast<int>(sol->YS[5]*scale)};
        std::array<unsigned,4> arc_startpoints_x_u=convert_unsigned_scale(sol->internalpath.arc_startpoints_x,scale);
        std::array<unsigned,4> arc_startpoints_y_u=convert_unsigned_scale(sol->internalpath.arc_startpoints_y,scale);
        std::array<unsigned,4> arc_endpoints_x_u=convert_unsigned_scale(sol->internalpath.arc_endpoints_x,scale);
        std::array<unsigned,4> arc_endpoints_y_u=convert_unsigned_scale(sol->internalpath.arc_endpoints_y,scale);
        #pragma GCC diagnostic warning "-Wnarrowing"

        //STARTPOINT/ENDPOINT
        std::string startpoint_circle = "<circle cx=\"" + std::to_string(startPoint_u[0]) + "\" cy=\"" + std::to_string(startPoint_u[1]) + "\" r=\"4\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n";
        std::string endpoint_circle = "<circle cx=\"" + std::to_string(endPoint_u[0]) + "\" cy=\"" + std::to_string(endPoint_u[1]) + "\" r=\"4\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n";
        content.push_back(startpoint_circle);
        content.push_back(endpoint_circle);

        //STRAIGHT LINES:
        // Construct the line strings. I tried replacing this with a for loop. It didn't work.
        std::array<std::string,5> lines;
        lines[0] = "<line x1=\"" + std::to_string(startPoint_u[0]) + "\" y1=\"" + std::to_string(startPoint_u[1]) +
                   "\" x2=\"" + std::to_string(arc_startpoints_x_u[0]) + "\" y2=\"" + std::to_string(arc_startpoints_y_u[0]) +
                   "\" style=\"stroke:black;stroke-width:2\" />\n";
        lines[1] = "<line x1=\"" + std::to_string(arc_endpoints_x_u[0]) + "\" y1=\"" + std::to_string(arc_endpoints_y_u[0]) +
                   "\" x2=\"" + std::to_string(arc_startpoints_x_u[1]) + "\" y2=\"" + std::to_string(arc_startpoints_y_u[1]) +
                   "\" style=\"stroke:black;stroke-width:2\" />\n";
        lines[2] = "<line x1=\"" + std::to_string(arc_endpoints_x_u[1]) + "\" y1=\"" + std::to_string(arc_endpoints_y_u[1]) +
                   "\" x2=\"" + std::to_string(arc_startpoints_x_u[2]) + "\" y2=\"" + std::to_string(arc_startpoints_y_u[2]) +
                   "\" style=\"stroke:black;stroke-width:2\" />\n";
        lines[3] = "<line x1=\"" + std::to_string(arc_endpoints_x_u[2]) + "\" y1=\"" + std::to_string(arc_endpoints_y_u[2]) +
                   "\" x2=\"" + std::to_string(arc_startpoints_x_u[3]) + "\" y2=\"" + std::to_string(arc_startpoints_y_u[3]) +
                   "\" style=\"stroke:black;stroke-width:2\" />\n";
        lines[4] = "<line x1=\"" + std::to_string(arc_endpoints_x_u[3]) + "\" y1=\"" + std::to_string(arc_endpoints_y_u[3]) +
                   "\" x2=\"" + std::to_string(endPoint_u[0]) + "\" y2=\"" + std::to_string(endPoint_u[1]) +
                   "\" style=\"stroke:black;stroke-width:2\" />\n";

        // Add the constructed lines to the content vector
        for(int i=0;i<5;i++){
        content.push_back(lines[i]);}

        //ARCS: Processed in their own function
        std::array<std::string,4> arcs = formatarc_svg(sol, scale);
        for(int i=0;i<4;i++){
            content.push_back(arcs[i]);
        }

        //Guide Points
        for(int i=0;i<10;i++){
            std::string circleElement = "<circle cx=\"" + std::to_string(static_cast<int>(sol->XS[i] * scale))
                                + "\" cy=\"" + std::to_string(static_cast<int>(sol->YS[i] * scale))
                                + "\" r=\"3\" fill=\"green\" />\n";
            content.push_back(circleElement);
        }
        //for debugging purposes: center points
        for(int i=0;i<4;i++)
           { std::string circleElement = "<circle cx=\"" + std::to_string(static_cast<int>(sol->internalpath.m[i][0] * scale))
                                + "\" cy=\"" + std::to_string(static_cast<int>(sol->internalpath.m[i][1] * scale))
                                + "\" r=\"3\" fill=\"orange\" />\n";
            content.push_back(circleElement);}

        //this line ends the svg file.
        content.push_back("</svg>");

        // OUTPUT FILE
        std::ofstream outputFile("map.svg");

        if (outputFile.is_open()) {
            // Write the modified content back to the file
            for (const auto& line : content) { outputFile << line << std::endl;}

            // Close the output file
            outputFile.close();

            //std::cout << "Successfully modified the SVG file." << std::endl;
        } else {
            std::cout << "Failed to open the output file." << std::endl;
        }
    } else {
        std::cout << "Failed to open the input file." << std::endl;
    }

    return 0;
}


// Function to calculate points along the arc => I'm leaving the legacy version in because discretize calls it
std::vector<point> calculatePointsOnArc_legacy(const point start, const point center, const point end,
                                        int direction, double distance_betw_points,double addtonext, double &arcLength) {
    std::vector<point> arcPoints;

    // Calculate the radius of the circle formed by the arc
    double radius = std::hypot(center[0] - start[0], center[1] - start[1]);

    // Calculate the start and end angles (in radians)
    double startAngle = atan2(start[1] - center[1], start[0] - center[0]);
    double endAngle = atan2(end[1] - center[1], end[0] - center[0]);
    double addnextAngle=addtonext/radius; //angle from adding to next; I hope all of this is in radiants
    //std::cout<<"addnextAngle: "<<addnextAngle<<std::endl;
    
    // Adjust the angles for clockwise direction
    if (direction == -1) {
        if (endAngle < startAngle) {
            endAngle += 2 * 3.1415926;
        }
        startAngle = startAngle + addnextAngle;
    } else {
        // Adjust the end angle if needed for counterclockwise direction
        if (endAngle > startAngle) {
            endAngle -= 2 * 3.1415926;
        }
        startAngle = startAngle - addnextAngle;
    }

    // Calculate the number of points on the arc
    int numberOfPoints = static_cast<int>(std::ceil(arcLength / distance_betw_points));

    // Calculate the increment in angle for each point
    double angleIncrement = (endAngle - startAngle) / numberOfPoints;

    // Generate points along the arc
    for (int i = 0; i <= numberOfPoints; ++i) {
        double currentAngle = startAngle + i * angleIncrement;
        double x = center[0] + radius * cos(currentAngle);
        double y = center[1] + radius * sin(currentAngle);

        arcPoints.push_back({x, y});
    }

    return arcPoints;
}

// Function to calculate points along the arc
std::vector<point> calculatePointsOnArc(const point start, const point center, const point end,
                                        int direction, double distance_betw_points,double addtonext, double &arcLength) {
    std::vector<point> arcPoints;

    // Calculate the radius of the circle formed by the arc
    double radius = euclid_dist(start,center);
    arcLength=getArcLength(start,center,end);

    // Calculate the start and end angles (in radians)
    double startAngle = atan2(start[1] - center[1], start[0] - center[0]);
    if(startAngle<0){startAngle+=2*3.1415926;} //make sure all angles are positive
    double endAngle = atan2(end[1] - center[1], end[0] - center[0]);
    if(endAngle<0){endAngle+=2*3.1415926;}
    //std::cout<<"start angle: "<<startAngle<<", end angle: "<<endAngle<<"\n";
    double addnextAngle=addtonext/radius; //angle from adding to next; I hope all of this is in radiants

    // Adjust the angles for clockwise direction and  endangle<startangle
    if (direction == -1) {
        startAngle = startAngle - addnextAngle;
    } else {
        // Adjust the end angle if needed for counterclockwise direction
        if (endAngle < startAngle) {
            endAngle += 2 * 3.1415926;
        }
        startAngle = startAngle + addnextAngle;
    }
    //std::cout<<"after angle correction - start angle: "<<startAngle<<", end angle: "<<endAngle<<"\n";

    // Calculate the number of points on the arc
    int numberOfPoints = static_cast<int>(std::ceil(arcLength / distance_betw_points));

    // Calculate the increment in angle for each point
    double angleIncrement = (endAngle - startAngle) / numberOfPoints;

    // Generate points along the arc
    for (int i = 0; i <= numberOfPoints; ++i) {
        double currentAngle = startAngle + i * angleIncrement;
        //std::cout<<"currentAngle: "<<currentAngle<<"\n";
        double x = center[0] + radius * cos(currentAngle);
        double y = center[1] + radius * sin(currentAngle);

        arcPoints.push_back({x, y});
    }

    return arcPoints;
}

//This function converts the path_struct into a point cloud
std::vector<point> pathtopointcloud(SMAsol_struct* sol, double distance_betw_points, model_struct model){
    //step 0: Set up vectors for x and y direction, preallocate arrays
//    std::cout<<"Entered step 0\n";
    std::array<point,5> line_start, line_end;
    line_start[0]={sol->XS[0],sol->YS[0]};
    for(int i=0;i<4;i++){
        line_start[i+1]={sol->internalpath.arc_endpoints_x[i],sol->internalpath.arc_endpoints_y[i]};
        line_end[i]={sol->internalpath.arc_startpoints_x[i],sol->internalpath.arc_startpoints_y[i]};
    }
    line_end[4]={sol->XS[5],sol->YS[5]};
    std::array<point,4> arc_start,arc_end;
    for(int i=0;i<4;i++){
        arc_start[i][0]=sol->internalpath.arc_startpoints_x[i];
        arc_start[i][1]=sol->internalpath.arc_startpoints_y[i];
        arc_end[i][0]=sol->internalpath.arc_endpoints_x[i];
        arc_end[i][1]=sol->internalpath.arc_endpoints_y[i];
    }
    std::array<point,5> linevector;
    std::array<double,5> linelength;
    std::array<double,4> arclength;

    //step 1: Get lengths of straight lines and arcs
    //std::cout<<"Entered step 1\n";
    for(int i=0;i<5;i++){
    linevector[i][0] = line_end[i][0] - line_start[i][0];
    linevector[i][1] = line_end[i][1] - line_start[i][1];
    //alternating between straight lines and arcs; even for lines, uneven for arcs
    linelength[i] = sqrt(pow(linevector[i][0], 2) + pow(linevector[i][1], 2));
    //set euclididan norm of linevectors to 1 => transform into unit vectors
    linevector[i][0] /=linelength[i]; //set length to 1
    linevector[i][1] /=linelength[i];
    }
    //get arc lengths
    for(int i=0;i<4;i++){
    arclength[i]=getArcLength(arc_start[i], {sol->internalpath.m[i][0], sol->internalpath.m[i][1]}, arc_end[i]);
    //std::cout<<"arclength["<<i<<"]: "<<arclength[i]<<"\n";
    }
    //step 2 became obsolete
    //step 3: Call functions to calculate point positions
    //std::cout<<"Entered step 3\n";
    std::vector<point> pointcloud;
    std::array<unsigned,5> amountofpoints_line;
    int k=0;
    //step 3.1: let a vector go from start+addtonext to last point and append it to a vector
    point runner;
    std::vector<point> arcpoints;
    double dist_end=0, dist_start=0;
    for(int k=0;k<5;k++){
        amountofpoints_line[k] = floor( (linelength[k]-dist_start) / distance_betw_points);
        if(dist_start<0){break;} //error handling because sometimes dist_start is negative (amountofpoints==0)
        //std::cout<<"k: "<<k<<", amountofpoints_line[k]: "<<amountofpoints_line[k]<<"\n";
        runner = {line_start[k][0]+linevector[k][0]*dist_start -linevector[k][0]*distance_betw_points, line_start[k][1]+linevector[k][1]*dist_start-linevector[k][1]*distance_betw_points};
        if(amountofpoints_line[k]>=1 && amountofpoints_line[k]<=10000){ //<=10000 because when dis_start<0, amountofpoints can go to infinity
            for(int j=0;j<=amountofpoints_line[k];j++) 
            {
                runner = {runner[0]+linevector[k][0]*distance_betw_points, runner[1]+linevector[k][1]*distance_betw_points};
                pointcloud.push_back(runner);
                //std::cout<<"line: pushing back ("<<runner[0]<<","<<runner[1]<<")\n";
                if (runner[0]>=model.ub || runner[1]>=model.ub) {
                    throw std::runtime_error("Point cloud creation failed"); //todo: find better solution here
                }
            }
            dist_end = euclid_dist(runner,line_end[k]); //distance between last point on line and end point of line
            dist_start = distance_betw_points - dist_end;
            assert(dist_start<=distance_betw_points);
            assert(dist_end<=distance_betw_points);
        }
        //set last point of last line to endPoint of model
        if(k==4){pointcloud.pop_back();pointcloud.push_back(line_end[4]);}
        //step 3.2: Take care of arc; add points on arcs, but only from 0 to 3
        if(k<4){
        double uselessdummy; //I also call calculatePointsonArc from discretize, there I need each arc's length, here I don't need it
        //std::cout<<"arc_start[k]: ("<<arc_start[k][0]<<","<<arc_start[k][1]<<"), center: ("<<path.arc_center_points_x[k]<<","<<path.arc_center_points_y[k]<<"), arc_end[k]: ("<<arc_end[k][0]<<","<<arc_end[k][1]<<"), dist_start: "<<dist_start<<"), dist_end: "<<dist_end<<"\n";
        point center={sol->internalpath.m[k][0],sol->internalpath.m[k][1]};
        arcpoints= calculatePointsOnArc(arc_start[k], center,arc_end[k],
                                        sol->internalpath.arc_direction[k], distance_betw_points,dist_start,uselessdummy);
        //std::cout<<"amount of points on arc "<<k<<": "<<arcpoints.size()<<"\n";
        for(int i=0;i<arcpoints.size();i++)
        {
            pointcloud.push_back(arcpoints[i]);
        }
        }
        dist_end = euclid_dist(arcpoints.back(),arc_end[k]); //distance between last point on segment (line/arc) and end point of that seegment
        dist_start = distance_betw_points - dist_end; //distance between start point on segment (line/arc) and first point on that segment
        //std::cout<<"dist_start: "<<dist_start<<", dist_end: "<<dist_end<<"\n";
        
    }
    return pointcloud;
}

int drawpointcloud_svg(std::vector<point> pointcloud, unsigned scale)
{
    //INPUT FILE
    std::ifstream inputFile("map.svg");

    if (inputFile.is_open()) 
    {
        // Read the content of the file
        std::vector<std::string> content;
        std::string line;
        while (std::getline(inputFile, line)) {
            if (line.find("</svg>") != std::string::npos) {
                break; // Stop reading when </svg> is found so that the line ending the svg would not be read
            }
            content.push_back(line);
        }
        // Close the input file
        inputFile.close();

        //put in the points
        std::string point_circle;
        for(int i=0; i<pointcloud.size();i++)
        {
        point_circle = "<circle cx=\"" + std::to_string(pointcloud[i][0]*scale) + "\" cy=\"" + std::to_string(pointcloud[i][1]*scale) + "\" r=\"4\" fill=\"none\" stroke=\"#66BB6A\" stroke-width=\"2\" />\n";
        content.push_back(point_circle);
        }

        //this line ends the svg file.
        content.push_back("</svg>");

        // OUTPUT FILE
        std::ofstream outputFile("map.svg");
        if (outputFile.is_open()) {
            // Write the modified content back to the file
            for (const auto& line : content) {
                outputFile << line << std::endl;
            }
            outputFile.close();
        } else {
            std::cout << "Failed to open the output file." << std::endl;
        }
    }
    return 0;
}

//This function writes the path into aa csv file
int pathtocsv(SMAsol_struct* sol)
{
    //prepare data to be put out
    //step 0: Set up vectors for x and y direction, preallocate arrays
    std::array<double,4> startp_x={sol->internalpath.arc_startpoints_x[0],sol->internalpath.arc_startpoints_x[1],sol->internalpath.arc_startpoints_x[2],sol->internalpath.arc_startpoints_x[3]};
    std::array<double,4> startp_y={sol->internalpath.arc_startpoints_y[0],sol->internalpath.arc_startpoints_y[1],sol->internalpath.arc_startpoints_y[2],sol->internalpath.arc_startpoints_y[3]};
    std::array<double,5> endp_x={sol->XS[0],sol->internalpath.arc_endpoints_x[0],sol->internalpath.arc_endpoints_x[1],sol->internalpath.arc_endpoints_x[2],sol->internalpath.arc_endpoints_x[3]}; //add the StartPoint here because at this point a straight line starts
    std::array<double,5> endp_y={sol->YS[0],sol->internalpath.arc_endpoints_y[0],sol->internalpath.arc_endpoints_y[1],sol->internalpath.arc_endpoints_y[2],sol->internalpath.arc_endpoints_y[3]};

    //open output file
    std::ofstream outputFile("path.csv");
    if (outputFile.is_open()) {
        outputFile << "Type, x Startpoint, y Startpoint, x Endpoint, y Endpoint, x Center Point, y Center Point, Clockwise(-1) / Counterclockwise (1)" << std::endl; // CSV header
        //for loop
        for(int i=0;i<4;i++){
        //write line
        outputFile << "line" << "," << endp_x[i] << "," << endp_y[i] << "," << startp_x[i] << "," << startp_y[i] << std::endl;
        //write arc
        outputFile << "arc" << "," << startp_x[i] << "," << startp_y[i] << "," << endp_x[i+1] << "," << endp_y[i+1] << "," << sol->internalpath.m[i][0] << "," << sol->internalpath.m[i][1] << "," << sol->internalpath.arc_direction[i] << std::endl;
        }
        // Close the file after writing
        outputFile.close();
       // std::cout << "Data has been written to data.csv" << std::endl;
    } else {
        std::cout << "Unable to open CSV file to write into" << std::endl;
    }
    return 0;
}
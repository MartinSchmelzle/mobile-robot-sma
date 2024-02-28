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

void printPathStruct(const path_struct& path) {
    std::cout << "Start Point: (" << path.startPoint[0] << ", " << path.startPoint[1] << ")\n";
    std::cout << "End Point: (" << path.endPoint[0] << ", " << path.endPoint[1] << ")\n";
    
    for (int i = 0; i < 4; ++i) {
        std::cout << "Arc " << i+1 << ":\n";
        std::cout << "    Start X: " << path.arc_startpoints_x[i] << "\n";
        std::cout << "    Start Y: " << path.arc_startpoints_y[i] << "\n";
        std::cout << "    End X: " << path.arc_endpoints_x[i] << "\n";
        std::cout << "    End Y: " << path.arc_endpoints_y[i] << "\n";
        std::cout << "    Center X: " << path.arc_center_points_x[i] << "\n";
        std::cout << "    Center Y: " << path.arc_center_points_y[i] << "\n";
        std::cout << "    Radius: " << path.arc_radii[i] << "\n";
        std::cout << "    Counterclockwise: " << path.arc_counterclockwise[i] << "\n\n";
    }
}

//This function takes in a path struct and returns a string array containing the svg lines.
//Writing to svg requires reformatting the arcs
std::array <std::string,4> formatarc_svg(path_struct path, unsigned scale)
{
    arr4 arc_radius;
    std::array <std::string,4> lines;
    //Modify arc data (scale + conversion to unsigned)
    #pragma GCC diagnostic ignored "-Wnarrowing"//Ignore warning about datatype conversion
    std::array<unsigned, 2> startPoint_u={static_cast<unsigned>(path.startPoint[0]*scale), static_cast<int>(path.startPoint[1]*scale)}; //scale to internal representation, round to integer
    std::array<unsigned, 2> endPoint_u={static_cast<unsigned>(path.endPoint[0]*scale), static_cast<int>(path.endPoint[1]*scale)};
    std::array<unsigned,4> arc_startpoints_x_u=convert_unsigned_scale(path.arc_startpoints_x,scale);
    std::array<unsigned,4> arc_startpoints_y_u=convert_unsigned_scale(path.arc_startpoints_y,scale);
    std::array<unsigned,4> arc_endpoints_x_u=convert_unsigned_scale(path.arc_endpoints_x,scale);
    std::array<unsigned,4> arc_endpoints_y_u=convert_unsigned_scale(path.arc_endpoints_y,scale);
    #pragma GCC diagnostic warning "-Wnarrowing"

    for(int i=0;i<4;i++)
    {
        //calculate radius
        arc_radius[i] = path.arc_radii[i]*scale;
        //convert flag for clockwise/counterclockwise => in svg, it's 0/1, in the alg it's 1/0/-1
        std::array<bool,4> clock;
        if(path.arc_counterclockwise[i]>=0){clock[i]=false;}//includes 180°
        else{clock[i]=true;}
        //TO-DO:calculate how to set the large arc flag
        //point vectorStart = {path.arc_startpoints_x[i] - path.arc_center_points_x[0], path.arc_startpoints_y[i] - path.arc_center_points_y[0]};
        //point vectorEnd = {path.arc_endpoints_x[i] - path.arc_center_points_x[0], path.arc_endpoints_y[i] - path.arc_center_points_y[0]};
        //double crossProduct = vectorStart[0] * vectorEnd[1] - vectorStart[1] * vectorEnd[0];
        point start_i={path.arc_startpoints_x[i],path.arc_startpoints_y[i]};
        point center_i={path.arc_center_points_x[i],path.arc_center_points_y[i]};
        point end_i={path.arc_endpoints_x[i],path.arc_endpoints_y[i]};
        double crossProduct = arc_cross(start_i,center_i,end_i);
        bool large_arc_flag;
        if (clock[i]==false) { //clockwise
            large_arc_flag = (crossProduct < 0) ? 0 : 1;
        } else {
            large_arc_flag = (crossProduct > 0) ? 1 : 0;
        }
        //edit from 05-01: I found out there are errors in the way the large arc flag is set.
        //todo: fix large arc flag. Given that we don't really have arcs with angles>180°, this might not be necessary.
        //provisory:
        large_arc_flag = 0;

        //write line for arc
        lines[i] = "<path d=\"M " + std::to_string(arc_startpoints_x_u[i]) + "," + std::to_string(arc_startpoints_y_u[i]) +
                          " A " + std::to_string(arc_radius[i]) + "," + std::to_string(arc_radius[i]) + " 0 " + std::to_string(large_arc_flag) + " " + std::to_string(clock[i]) + " " +
                          std::to_string(arc_endpoints_x_u[i]) + "," + std::to_string(arc_endpoints_y_u[i]) +
                          "\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n";;
    }
    return lines;
}


//This function draws a given path of type path struct into the svg file "map.svg".
int drawpathintosvg(path_struct path,unsigned scale)
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
        std::array<unsigned, 2> startPoint_u={static_cast<int>(path.startPoint[0]*scale), static_cast<int>(path.startPoint[1]*scale)}; //scale to internal representation, round to integer
        std::array<unsigned, 2> endPoint_u={static_cast<int>(path.endPoint[0]*scale), static_cast<int>(path.endPoint[1]*scale)};
        std::array<unsigned,4> arc_startpoints_x_u=convert_unsigned_scale(path.arc_startpoints_x,scale);
        std::array<unsigned,4> arc_startpoints_y_u=convert_unsigned_scale(path.arc_startpoints_y,scale);
        std::array<unsigned,4> arc_endpoints_x_u=convert_unsigned_scale(path.arc_endpoints_x,scale);
        std::array<unsigned,4> arc_endpoints_y_u=convert_unsigned_scale(path.arc_endpoints_y,scale);
        #pragma GCC diagnostic warning "-Wnarrowing"

        //STARTPOINT/ENDPOINT
        std::string startpoint_circle = "<circle cx=\"" + std::to_string(startPoint_u[0]) + "\" cy=\"" + std::to_string(startPoint_u[1]) + "\" r=\"4\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n";
        std::string endpoint_circle = "<circle cx=\"" + std::to_string(endPoint_u[0]) + "\" cy=\"" + std::to_string(endPoint_u[1]) + "\" r=\"4\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n";
        content.push_back(startpoint_circle);
        content.push_back(endpoint_circle);

        //STRAIGHT LINES:
        // Construct the line strings
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
        std::array<std::string,4> arcs = formatarc_svg(path, scale);
        content.push_back(arcs[0]);
        content.push_back(arcs[1]);
        content.push_back(arcs[2]);
        content.push_back(arcs[3]);

        //Guide Points
        for(int i=0;i<10;i++){
            std::string circleElement = "<circle cx=\"" + std::to_string(static_cast<int>(path.XS[i] * scale))
                                + "\" cy=\"" + std::to_string(static_cast<int>(path.YS[i] * scale))
                                + "\" r=\"3\" fill=\"green\" />\n";
            content.push_back(circleElement);
        }
        //for debugging purposes: center points
        for(int i=0;i<4;i++)
           { std::string circleElement = "<circle cx=\"" + std::to_string(static_cast<int>(path.arc_center_points_x[i] * scale))
                                + "\" cy=\"" + std::to_string(static_cast<int>(path.arc_center_points_y[i] * scale))
                                + "\" r=\"3\" fill=\"orange\" />\n";
            content.push_back(circleElement);}

        //this line ends the svg file.
        content.push_back("</svg>");

        // OUTPUT FILE
        std::ofstream outputFile("map.svg");

        if (outputFile.is_open()) {
            // Write the modified content back to the file
            for (const auto& line : content) {
                outputFile << line << std::endl;
            }

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

    // Calculate the arc length
    //arcLength = radius * std::abs(endAngle - startAngle);
    //std::cout<<"endAngle: "<<endAngle<<"startAngle: "<<startAngle<<std::endl;
    //std::cout<<"arclength: "<<arcLength<<"\n";

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
std::vector<point> pathtopointcloud(path_struct path, double distance_betw_points, model_struct model){
    //step 0: Set up vectors for x and y direction, preallocate arrays
//    std::cout<<"Entered step 0\n";
    std::array<point,5> line_start, line_end;
    line_start[0]=path.startPoint;
    line_start[1]={path.arc_endpoints_x[0],path.arc_endpoints_y[0]};
    line_start[2]={path.arc_endpoints_x[1],path.arc_endpoints_y[1]};
    line_start[3]={path.arc_endpoints_x[2],path.arc_endpoints_y[2]};
    line_start[4]={path.arc_endpoints_x[3],path.arc_endpoints_y[3]};
    line_end[0]={path.arc_startpoints_x[0],path.arc_startpoints_y[0]};
    line_end[1]={path.arc_startpoints_x[1],path.arc_startpoints_y[1]};
    line_end[2]={path.arc_startpoints_x[2],path.arc_startpoints_y[2]};
    line_end[3]={path.arc_startpoints_x[3],path.arc_startpoints_y[3]};
    line_end[4]=path.endPoint;
    std::array<point,4> arc_start,arc_end;
    for(int i=0;i<4;i++){
        arc_start[i][0]=path.arc_startpoints_x[i];
        arc_start[i][1]=path.arc_startpoints_y[i];
        arc_end[i][0]=path.arc_endpoints_x[i];
        arc_end[i][1]=path.arc_endpoints_y[i];
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
    arclength[i]=getArcLength(arc_start[i], {path.arc_center_points_x[i], path.arc_center_points_y[i]}, arc_end[i]);
    //std::cout<<"arclength["<<i<<"]: "<<arclength[i]<<"\n";
    }

    //step 2: Calculate how many points fit on each line
    //std::cout<<"Entered step 2\n";
    std::array<unsigned,5> amountofpoints_line;
    int k=0;

    //step 3: Call functions to calculate point positions
    //std::cout<<"Entered step 3\n";
    std::vector<point> pointcloud;
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
        point center={path.arc_center_points_x[k],path.arc_center_points_y[k]};
        arcpoints= calculatePointsOnArc(arc_start[k], center,arc_end[k],
                                        arc_direction(arc_start[k],center,arc_end[k]), distance_betw_points,dist_start,uselessdummy);
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

//This function converts the path_struct into a point cloud
/*std::vector<point> pathtopointcloud_legacy(path_struct path, double distance_betw_points, model_struct model){
    //step 0: Set up vectors for x and y direction, preallocate arrays
//    std::cout<<"Entered step 0\n";
    std::array<point,5> line_start, line_end;
    line_start[0]=path.startPoint;
    line_start[1]={path.arc_endpoints_x[0],path.arc_endpoints_y[0]};
    line_start[2]={path.arc_endpoints_x[1],path.arc_endpoints_y[1]};
    line_start[3]={path.arc_endpoints_x[2],path.arc_endpoints_y[2]};
    line_start[4]={path.arc_endpoints_x[3],path.arc_endpoints_y[3]};
    line_end[0]={path.arc_startpoints_x[0],path.arc_startpoints_y[0]};
    line_end[1]={path.arc_startpoints_x[1],path.arc_startpoints_y[1]};
    line_end[2]={path.arc_startpoints_x[2],path.arc_startpoints_y[2]};
    line_end[3]={path.arc_startpoints_x[3],path.arc_startpoints_y[3]};
    line_end[4]=path.endPoint;
    std::array<point,4> arc_start,arc_end;
    for(int i=0;i<4;i++){
        arc_start[i][0]=path.arc_startpoints_x[i];
        arc_start[i][1]=path.arc_startpoints_y[i];
        arc_end[i][0]=path.arc_endpoints_x[i];
        arc_end[i][1]=path.arc_endpoints_y[i];
    }
    std::array<point,5> linevector;
    std::array<double,5> linelength;
    std::array<double,4> arclength;

    //step 1: Get lengths of straight lines and arcs
    std::cout<<"Entered step 1\n";
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
    arclength[i]=getArcLength(arc_start[i], {path.arc_center_points_x[i], path.arc_center_points_y[i]}, arc_end[i]);
    std::cout<<"arclength["<<i<<"]: "<<arclength[i]<<"\n";
    }

    //step 1.5:calculate total length, number of points
    std::cout<<"Entered step 1.5\n";
    double totallength=0; //length of all segments together
    for (int i = 0; i < 5; ++i) { //equivalent to sum function
        totallength += linelength[i];
        if(i<4){totallength+=arclength[i];}
    }

    //std::cout<<"total length, "<<totallength<<std::endl;
    unsigned totalpointnumber= ceil(totallength / distance_betw_points);

    //step 2: Calculate how many points fit on each segment, how much is left at the end and from that, how much needs to be added to the next segment
    std::cout<<"Entered step 2\n";
    unsigned i=0;
    std::array<unsigned,5> amountofpoints_line;
    std::array<unsigned,4> amountofpoints_arc;
    std::array<double,9> dis_start{}, dis_end{};
    int k=0;
    amountofpoints_line[0] = floor(linelength[0] / distance_betw_points);
    dis_start[0]=0;
    dis_end[0] = linelength[0] - amountofpoints_line[0]*distance_betw_points; //distance left after last point
//    std::cout<<"dis_end[0]: "<<dis_end[0]<<"\n";
    //assert(dis_end[0]<distance_betw_points);
    for(int i=0;i<4;i++)
    {
        //arc
        dis_start[2*i+1]=distance_betw_points-dis_end[2*i];
        //assert(dis_end[2*i]<distance_betw_points&&dis_end[2*i]>=0);
//        std::cout<<"dis_start[2*i+1]: "<<dis_start[2*i+1]<<"\narclength[i]: "<<arclength[i]<<"\n";
        //assert(dis_start[2*i+1]<distance_betw_points && dis_start[2*i+1]>=0);
        if(arclength[i]>dis_start[2*i+1]){
            amountofpoints_arc[i]=floor((arclength[i]-dis_start[2*i+1])/distance_betw_points);
        }
        else{amountofpoints_arc[i]=0;}
        dis_end[2*i+1] = arclength[i]-amountofpoints_arc[i]*distance_betw_points;
        if(dis_end[2*i+1]>=distance_betw_points){dis_end[2*i+1]-=distance_betw_points;amountofpoints_arc[i]+=1;}
//        std::cout<<"amountofpoints["<<i<<"]: "<<amountofpoints_arc[i]<<"\n";
//        std::cout<<"dis_end["<<2*i+1<<"]: "<<dis_end[2*i+1]<<"\n";
        //line
        dis_start[2*i+2]=distance_betw_points-dis_end[2*i+1];
//        std::cout<<"dis_start["<<2*i+2<<"]: "<<dis_start[2*i+2]<<"\n";
        if(linelength[i+1]>dis_start[2*i+2]){
            amountofpoints_line[i+1]=floor((linelength[i+1]-dis_start[2*i+2])/distance_betw_points);
        }
        else{amountofpoints_line[i+1]=0;}
        dis_end[2*i+2] = linelength[i+1]-amountofpoints_line[i+1]*distance_betw_points;
//        std::cout<<"dis_end["<<2*i+2<<"]: "<<dis_end[2*i+2]<<"\n";
        if(dis_end[2*i+2]>=distance_betw_points){dis_end[2*i+2]-=distance_betw_points;amountofpoints_line[i+1]+=1;}
//        std::cout<<"dis_end["<<2*i+2<<"]: "<<dis_end[2*i+2]<<"\n";
        //assert(dis_end[2*i+1]<distance_betw_points && dis_end[2*i+1]>=0);
        //assert(dis_end[2*i+2]<distance_betw_points && dis_end[2*i+2]>=0);
        //assert(dis_start[2*i+2]<distance_betw_points&&dis_start[2*i+2]>=0);
    }

    //step 3: Call functions to calculate point positions
    std::cout<<"Entered step 3\n";
    std::vector<point> pointcloud;
    //let a vector go from start+addtonext to last point and append it to a vector
    point runner;
    std::vector<point> arcpoints;

    for(int k=0;k<5;k++){
        //std::cout<<"k="<<k<<"\n";
        //std::cout<<"amount of points for lines: "<<amountofpoints_line[k]<<"\ndis_start: "<<dis_start[2*k];
        //std::cout<<"line vector: ("<<linevector[k][0]<<","<<linevector[k][1]<<")\n";
        //std::cout<<"line start: ("<<line_start[k][0]<<","<<line_start[k][1]<<")\n";
        //add points on straight lines
        runner = {line_start[k][0]+linevector[k][0]*dis_start[2*k] , line_start[k][1]+linevector[k][1]*dis_start[2*k]};
        for(int j=0;j<=amountofpoints_line[k]-1;j++) 
        {
            runner = {runner[0]+linevector[k][0]*distance_betw_points, runner[1]+linevector[k][1]*distance_betw_points};
            pointcloud.push_back(runner);
            //std::cout<<"line: pushing back ("<<runner[0]<<","<<runner[1]<<")\n";
            if (runner[0]>=model.ub || runner[1]>=model.ub) {
                throw std::runtime_error("Point cloud creation failed"); //todo: find better solution here
            }
        }
        //set last point of last line to endPoint of model
        if(k==4){pointcloud.pop_back();pointcloud.push_back(line_end[4]);}
        //add points on arcs, but only from 0 to 3
        if(k<4){
        double uselessdummy; //I also call calculatePointsonArc from discretize, there I need each arc's length, here I don't need it
        arcpoints= calculatePointsOnArc(arc_start[k], {path.arc_center_points_x[k],path.arc_center_points_y[k]},
                                        arc_end[k],
                                        path.arc_counterclockwise[k], distance_betw_points,dis_start[2*k+1],uselessdummy);
        for(i=0;i<arcpoints.size();i++)
        {
            pointcloud.push_back(arcpoints[i]);
        }
        }
    }
    return pointcloud;
}*/

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

            // Close the output file
            outputFile.close();

            //std::cout << "Successfully modified the SVG file." << std::endl;
        } else {
            std::cout << "Failed to open the output file." << std::endl;
        }
    }
    return 0;
}

//This function writes the path into aa csv file
int pathtocsv(path_struct path)
{
    //prepare data to be put out
    //step 0: Set up vectors for x and y direction, preallocate arrays
    std::array<double,4> startp_x={path.arc_startpoints_x[0],path.arc_startpoints_x[1],path.arc_startpoints_x[2],path.arc_startpoints_x[3]};
    std::array<double,4> startp_y={path.arc_startpoints_y[0],path.arc_startpoints_y[1],path.arc_startpoints_y[2],path.arc_startpoints_y[3]};
    std::array<double,5> endp_x={path.startPoint[0],path.arc_endpoints_x[0],path.arc_endpoints_x[1],path.arc_endpoints_x[2],path.arc_endpoints_x[3]}; //add the StartPoint here because at this point a straight line starts
    std::array<double,5> endp_y={path.startPoint[1],path.arc_endpoints_y[0],path.arc_endpoints_y[1],path.arc_endpoints_y[2],path.arc_endpoints_y[3]};

    //open output file
    std::ofstream outputFile("path.csv");
    if (outputFile.is_open()) {
        outputFile << "Type, x Startpoint, y Startpoint, x Endpoint, y Endpoint, x Center Point, y Center Point, Clockwise(-1) / Counterclockwise (1)" << std::endl; // CSV header
        //for loop
        for(int i=0;i<4;i++){
        //write line
        outputFile << "line" << "," << endp_x[i] << "," << endp_y[i] << "," << startp_x[i] << "," << startp_y[i] << std::endl;
        //write arc
        outputFile << "arc" << "," << startp_x[i] << "," << startp_y[i] << "," << endp_x[i+1] << "," << endp_y[i+1] << "," << path.arc_center_points_x[i] << "," << path.arc_center_points_y[i] << "," << path.arc_counterclockwise[i] << std::endl;
        }
        // Close the file after writing
        outputFile.close();
       // std::cout << "Data has been written to data.csv" << std::endl;
    } else {
        std::cout << "Unable to open CSV file to write into" << std::endl;
    }

    return 0;
}

/*
int main()
{
    //I will just make a random path now that I can use for testing
    path_struct pathfromalg;
    pathfromalg.startPoint = {2.0, 2.0};
    pathfromalg.endPoint = {88.0, 26.0};
    pathfromalg.arc_startpoints_x={2,55,58,93};
    pathfromalg.arc_startpoints_y={45,48,25,21};
    pathfromalg.arc_endpoints_x={5,58,62,88};
    pathfromalg.arc_endpoints_y={48,45,21,26};
    pathfromalg.arc_center_points_x={5,55,62,93};
    pathfromalg.arc_center_points_y={45,45,25,26};
    pathfromalg.arc_counterclockwise={1,1,-1,-1};
    
    //This should be an alright path. 
    unsigned scale=10;
    drawpathintosvg(pathfromalg, scale);

    //Conversion of path_struct into different formats:
    //point cloud
    double distance_betw_points = 1.6;
    std::vector<point> pointcloud=pathtopointcloud(pathfromalg, distance_betw_points);
    //for(int i=0;i<pointcloud.size();i++)
    //{std::cout<<pointcloud[i][0]<<" "<<pointcloud[i][1]<<std::endl;}
    drawpointcloud_svg(pointcloud, scale);

    //csv-like table with lines and arcs
    pathtocsv(pathfromalg);

    return 0;
}
*/
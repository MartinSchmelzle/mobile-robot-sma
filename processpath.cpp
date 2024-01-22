#include<iostream>
#include <array>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "../07-12-cpp-manual-algorithm/custom_datatypes.h" //makes sure the code works together with alg

using point = std::array<double, 2>;
//This is the data structure that the algorithm puts out
// I will refer to this as representing the path as a struct.
/*struct path_struct { //uncomment for this code to run independently
        point startPoint;
        point endPoint;
        std::array<double, 4> arc_startpoints_x;
        std::array<double, 4> arc_startpoints_y;
        std::array<double, 4> arc_endpoints_x;
        std::array<double, 4> arc_endpoints_y;
        std::array<double, 4> arc_center_points_x;
        std::array<double, 4> arc_center_points_y;
        std::array<int, 4> arc_counterclockwise;

    };*/ 
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
        std::cout << "    Counterclockwise: " << path.arc_counterclockwise[i] << "\n\n";
    }
}

//This function takes in a path struct and returns a string array containing the svg lines.
//Writing to svg requires reformatting the arcs
std::array <std::string,4> formatarc_svg(path_struct path, unsigned scale)
{
    std::array <unsigned,4> arc_radius;
    std::array <std::string,4> lines;
    //Modify arc data (scale + conversion to unsigned)
    #pragma GCC diagnostic ignored "-Wnarrowing"//Ignore warning about datatype conversion
    std::array<unsigned, 2> startPoint_u={static_cast<unsigned>(path.startPoint[0]*scale), static_cast<int>(path.startPoint[1]*scale)}; //scale to internal representation, round to integer
    std::array<unsigned, 2> endPoint_u={static_cast<unsigned>(path.endPoint[0]*scale), static_cast<int>(path.endPoint[1]*scale)};
    std::array<unsigned,4> arc_startpoints_x_u={static_cast<unsigned>(path.arc_startpoints_x[0]*scale), static_cast<unsigned>(path.arc_startpoints_x[1]*scale), static_cast<unsigned>(path.arc_startpoints_x[2]*scale), static_cast<unsigned>(path.arc_startpoints_x[3]*scale)};
    std::array<unsigned,4> arc_startpoints_y_u={static_cast<unsigned>(path.arc_startpoints_y[0]*scale), static_cast<unsigned>(path.arc_startpoints_y[1]*scale), static_cast<unsigned>(path.arc_startpoints_y[2]*scale), static_cast<unsigned>(path.arc_startpoints_y[3]*scale)};
    std::array<unsigned,4> arc_endpoints_x_u={static_cast<unsigned>(path.arc_endpoints_x[0]*scale), static_cast<unsigned>(path.arc_endpoints_x[1]*scale), static_cast<unsigned>(path.arc_endpoints_x[2]*scale), static_cast<unsigned>(path.arc_endpoints_x[3]*scale)};
    std::array<unsigned,4> arc_endpoints_y_u={static_cast<unsigned>(path.arc_endpoints_y[0]*scale), static_cast<unsigned>(path.arc_endpoints_y[1]*scale), static_cast<unsigned>(path.arc_endpoints_y[2]*scale), static_cast<unsigned>(path.arc_endpoints_y[3]*scale)};
    #pragma GCC diagnostic warning "-Wnarrowing"

    for(int i=0;i<4;i++)
    {
        //calculate radius
        arc_radius[i] = sqrt( pow(path.arc_center_points_x[i]-path.arc_startpoints_x[i],2) + pow(path.arc_center_points_y[i]-path.arc_startpoints_y[i],2) ); //Pythagorean theorem, then multiply by scale, then conversion to unsigned
        //scale radius to svg format (int, *scale)
        arc_radius[i] = static_cast<unsigned>(arc_radius[i]*scale);
        //convert flag for clockwise/counterclockwise => in svg, it's 0/1, in the alg it's 1/0/-1
        std::array<bool,4> clock;
        if(path.arc_counterclockwise[i]>=0){clock[i]=false;}//includes 180°
        else{clock[i]=true;}
        //TO-DO:calculate how to set the large arc flag
        point vectorStart = {path.arc_startpoints_x[i] - path.arc_center_points_x[0], path.arc_startpoints_y[i] - path.arc_center_points_y[0]};
        point vectorEnd = {path.arc_endpoints_x[i] - path.arc_center_points_x[0], path.arc_endpoints_y[i] - path.arc_center_points_y[0]};
        double crossProduct = vectorStart[0] * vectorEnd[1] - vectorStart[1] * vectorEnd[0];
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
        std::array<unsigned,4> arc_startpoints_x_u={static_cast<unsigned>(path.arc_startpoints_x[0]*scale), static_cast<unsigned>(path.arc_startpoints_x[1]*scale), static_cast<unsigned>(path.arc_startpoints_x[2]*scale), static_cast<unsigned>(path.arc_startpoints_x[3]*scale)};
        std::array<unsigned,4> arc_startpoints_y_u={static_cast<unsigned>(path.arc_startpoints_y[0]*scale), static_cast<unsigned>(path.arc_startpoints_y[1]*scale), static_cast<unsigned>(path.arc_startpoints_y[2]*scale), static_cast<unsigned>(path.arc_startpoints_y[3]*scale)};
        std::array<unsigned,4> arc_endpoints_x_u={static_cast<unsigned>(path.arc_endpoints_x[0]*scale), static_cast<unsigned>(path.arc_endpoints_x[1]*scale), static_cast<unsigned>(path.arc_endpoints_x[2]*scale), static_cast<unsigned>(path.arc_endpoints_x[3]*scale)};
        std::array<unsigned,4> arc_endpoints_y_u={static_cast<unsigned>(path.arc_endpoints_y[0]*scale), static_cast<unsigned>(path.arc_endpoints_y[1]*scale), static_cast<unsigned>(path.arc_endpoints_y[2]*scale), static_cast<unsigned>(path.arc_endpoints_y[3]*scale)};
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


// Function to calculate the arc length given start, end, and center points
//needed for conversion to pointcloud
double getArcLength(double start_x, double start_y, double end_x, double end_y,
                          double center_x, double center_y, int arc_direction) {
    // Calculate the radius of the arc (distance from center to start or end point)
    double radius = sqrt(pow(center_x - start_x, 2) + pow(center_y - start_y, 2));

    // Calculate angle subtended by the arc using trigonometry
    double angle = atan2(end_y - center_y, end_x - center_x) -
                   atan2(start_y - center_y, start_x - center_x);

    // Adjust angle based on the arc direction (clockwise or counterclockwise)
    switch (arc_direction) {
    case -1: // Clockwise
        if (angle > 0) {
            angle -= 2 * M_PI;
        }
    case 1: // Counterclockwise
        if (angle < 0) {
            angle += 2 * M_PI;
        }
    case 0: //no angle
        angle=0;
    }

    // Calculate the arc length using formula s = r * theta
    double arc_length = fabs(radius * angle);
    return arc_length;
}

// Function to calculate points along the arc
std::vector<point> calculatePointsOnArc(const point start, const point center, const point end,
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
            endAngle += 2 * M_PI;
        }
        startAngle = startAngle + addnextAngle;
    } else {
        // Adjust the end angle if needed for counterclockwise direction
        if (endAngle > startAngle) {
            endAngle -= 2 * M_PI;
        }
        startAngle = startAngle - addnextAngle;
    }

    // Calculate the arc length
    arcLength = radius * std::abs(endAngle - startAngle);
    //std::cout<<"endAngle: "<<endAngle<<"startAngle: "<<startAngle<<std::endl;
    
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

//This function converts the path_struct into a point cloud
std::vector<point> pathtopointcloud(path_struct path, double distance_betw_points){

    //step 0: Set up vectors for x and y direction, preallocate arrays
    std::array<double,4> startp_x={path.arc_startpoints_x[0],path.arc_startpoints_x[1],path.arc_startpoints_x[2],path.arc_startpoints_x[3]};
    std::array<double,4> startp_y={path.arc_startpoints_y[0],path.arc_startpoints_y[1],path.arc_startpoints_y[2],path.arc_startpoints_y[3]};
    std::array<double,5> endp_x={path.startPoint[0],path.arc_endpoints_x[0],path.arc_endpoints_x[1],path.arc_endpoints_x[2],path.arc_endpoints_x[3]}; //add the StartPoint here because at this point a straight line starts
    std::array<double,5> endp_y={path.startPoint[1],path.arc_endpoints_y[0],path.arc_endpoints_y[1],path.arc_endpoints_y[2],path.arc_endpoints_y[3]};

    std::array<point,4> linevector;
    std::array<double,8> linelength;

    //step 1: Get lengths of straight lines and arcs
    for(int i=0;i<4;i++){
    linevector[i] = {startp_x[i] - endp_x[i],startp_y[i] - endp_y[i]};
    //alternating between straight lines and arcs; even for lines, uneven for arcs
    linelength[2*i] = sqrt(pow(linevector[i][0], 2) + pow(linevector[i][1], 2));
    linelength[2*i+1]=getArcLength(startp_x[i], startp_y[i], endp_x[i+1], endp_y[i+1], path.arc_center_points_x[i], path.arc_center_points_y[i], path.arc_counterclockwise[i]);
    //set euclididan norm of linevectors to 1 => transform into unit vectors
    linevector[i][0] /=linelength[2*i]; //set length to 1
    linevector[i][1] /=linelength[2*i];
    }

    //step 1.5:calculate total length, number of points
    double totallength=0; //length of all segments together
    for (int i = 0; i < 8; ++i) { //equivalent to sum function
        totallength += linelength[i];
    }
    //std::cout<<"total length, "<<totallength<<std::endl;
    unsigned totalpointnumber= ceil(totallength / distance_betw_points);

    //step 2: Calculate how many points fit on each segment, how much is left at the end and from that, how much needs to be added to the next segment
    unsigned i=0;
    std::array<unsigned,8> amountofpoints;
    std::array<double,8> dis_leftatend, addtonext;
    do{
    amountofpoints[i] = floor(linelength[i] / distance_betw_points);
    dis_leftatend[i] = linelength[i] - amountofpoints[i]*distance_betw_points; //distance left after last point
    addtonext[i] = distance_betw_points - dis_leftatend[i]; //how much needs to be added to next segment?
    i++;
    }while(i<8);

    //step 3: Call functions to calculate point positions
    std::vector<point> pointcloud;
    //let a vector go from start+addtonext to last point and append it to a vector
    point runner;
    std::vector<point> arcpoints;

    for(int k=0;k<4;k++){
    //add points on straight lines
    runner = {endp_x[k]+linevector[k][0]*addtonext[2*k-1],endp_y[k]+linevector[k][1]*addtonext[2*k-1]};
    pointcloud.push_back(runner);
    for(int j=0;j<amountofpoints[2*k];j++)
    {
        pointcloud.push_back(runner);
        runner = {runner[0]+linevector[k][0]*distance_betw_points,runner[1]+linevector[k][1]*distance_betw_points};
    }
    //add points on arcs
    double uselessdummy; //I also call calculatePointsonArc from discretize, there I need each arc's length, here I don't need it
    arcpoints= calculatePointsOnArc({startp_x[k],startp_y[k]}, 
                                    {path.arc_center_points_x[k],path.arc_center_points_y[k]},
                                    {endp_x[k+1],endp_y[k+1]},
                                    path.arc_counterclockwise[k], distance_betw_points,addtonext[2*k],uselessdummy);
    for(i=0;i<arcpoints.size();i++)
    {
        pointcloud.push_back(arcpoints[i]);
    }
    
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
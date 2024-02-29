#include <iostream>
#include <array>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include "../custom_datatypes.h" //makes sure the code works together with alg

//This function takes a user input color and gets an appropriate RGB value
//15 different colors supported
RGBColor color_mapping(std::string color_input){
    std::unordered_map<col, RGBColor> colorMap = {
        {col::red, {255, 0, 0}},
        {col::green, {0, 255, 0}},
        {col::blue, {0, 0, 255}},
        {col::orange, {255, 165, 0}},
        {col::yellow, {255, 255, 0}},
        {col::purple, {128, 0, 128}},
        {col::cyan, {0, 255, 255}},
        {col::magenta, {255, 0, 255}},
        {col::pink, {255, 192, 203}},
        {col::teal, {0, 128, 128}},
        {col::lime, {0, 255, 0}},
        {col::indigo, {75, 0, 130}},
        {col::gold, {255, 215, 0}},
        {col::silver, {192, 192, 192}},
        {col::maroon, {128, 0, 0}},
        {col::navy, {0, 0, 128}},
        {col::brown, {139, 69, 19}},
        {col::black, {0,0,0}},
        {col::white, {255,255,255}},
        // Add more color mappings...
    };
    RGBColor RGB_color;
    // Check the color input and assign the respective RGB values
    if (color_input == "red") {
        RGB_color = colorMap[col::red];
    } else if (color_input == "green") {
        RGB_color = colorMap[col::green];
    } else if (color_input == "blue") {
        RGB_color = colorMap[col::blue];
    } else if (color_input == "orange") {
        RGB_color = colorMap[col::orange];
    } else if (color_input == "yellow") {
        RGB_color = colorMap[col::yellow];
    } else if (color_input == "purple") {
        RGB_color = colorMap[col::purple];
    } else if (color_input == "cyan") {
        RGB_color = colorMap[col::cyan];
    } else if (color_input == "magenta") {
        RGB_color = colorMap[col::magenta];
    } else if (color_input == "indigo") {
        RGB_color = colorMap[col::indigo];
    } else if (color_input == "pink") {
        RGB_color = colorMap[col::pink];
    } else if (color_input == "teal") {
        RGB_color = colorMap[col::teal];
    } else if (color_input == "lime") {
        RGB_color = colorMap[col::lime];
    } else if (color_input == "gold") {
        RGB_color = colorMap[col::gold];
    } else if (color_input == "silver") {
        RGB_color = colorMap[col::silver];
    } else if (color_input == "maroon") {
        RGB_color = colorMap[col::maroon];
    } else if (color_input == "navy") {
        RGB_color = colorMap[col::navy];
    } else if (color_input == "brown") {
        RGB_color = colorMap[col::brown];
    } else if (color_input == "black") {
        RGB_color = colorMap[col::black];
    } else if (color_input == "white") {
        RGB_color = colorMap[col::white];
    } else {
        throw std::invalid_argument("Unsupported color provided!");
    }
    // Add more elseif blocks for other color inputs...
    //std::cout << RGB_color.red << ", " << RGB_color.green << ", " << RGB_color.blue << std::endl;
    return RGB_color;
}

//class that rectangle, circle etc. inherit from
class obstacle {
    arr2 position; // x and y coordinates of upper left corner
    int scale; //scale: How many pixels (or units in an svg file) are in one meter?
    RGBColor color;
    std::string label;
    //Internally, everything works in svg coordinates, but user wants to enter meters for everything
public://two constructors; one with label and color, another without
    obstacle(arr2 pos, unsigned s){
        // Constructor using initializer list to initialize position
        position[0]=pos[0]*static_cast<int>(s);
        position[1]= pos[1]*static_cast<int>(s);
        scale=s;
        color={0, 0,255}; //Silver by default
        label = std::string("obstacle");
    }

    obstacle(arr2 pos, unsigned s, std::string colortext, std::string Label){
        // Constructor using initializer list to initialize position
        position[0]=pos[0]*s;
        position[1]= pos[1]*s;
        scale=s;
        color=color_mapping(colortext); //Silver by default
        label = Label;
    }

    // Create at least one virtual function to make the thing polymorphic
    virtual ~obstacle() {} 

    //getters and setters
    void setpos(arr2 pos) {
        position[0] = pos[0]*scale;
        position[1] = pos[1]*scale;}
    arr2 getpos() {return {position[0]/scale,position[1]/scale};}
    int getscale() {return scale;}
    RGBColor getcolor() {return color;} //returns RGB value
    void setcolor(std::string color_input) {color=color_mapping(color_input);} //takes in string
    std::string getlabel() {return label;}
    void setlabel(std::string newlabel) {label = newlabel;}
};

class rectangle:public obstacle{ //in rectangle in svg, position refers to upper left corner
    arr2 size;
    public:
    //constructor with default values
    rectangle(arr2 pos, arr2 Size, float s):obstacle(pos,s)
    {
        size[0]=Size[0]*s;
        size[1]=Size[1]*s;
    }
    //second constructor with color, label
    rectangle(arr2 pos, arr2 Size, float s, std::string Colortext, std::string Label):obstacle(pos,s,Colortext, Label)
    {
        size[0]=Size[0]*s;
        size[1]=Size[1]*s;
    }

    arr2 getsize(){return {size[0]/getscale(),size[1]/getscale()};}
    void setsize(std::array<unsigned,2> Size){size[0]=Size[0]; size[1]=Size[1];}
};

class circle:public obstacle{ //For circles in .svg, position refers to the center!
    unsigned radius;
    public:
    circle(arr2 pos, unsigned Radius, float s):obstacle(pos,s)
    {radius=Radius*s;}
    circle(arr2 pos, unsigned Radius, float s,std::string Colortext, std::string Label):obstacle(pos,s,Colortext,Label)
    {radius=Radius*s;}
    unsigned getradius(){return radius/getscale();}
    void setradius(unsigned Radius){radius=Radius*getscale();}
};
//later: maybe also make triangles, circle sections etc. => everything that works in svg

//class for Charging Station obstacle
class ChargingStation:public obstacle{
    double angle,width=2,height=1; //angle relative to horizontal axis, clockwise
    std::array<turned_rectangle,3> blocked_areas;
    public:
    ChargingStation(arr2 pos, double Angle, float s):obstacle(pos,s) //defining two constructors again
    {
        angle=Angle;assert(angle>=0);assert(angle<360);//angle in degrees here
        blocked_areas=calculateBlockedAreas();
    }
    ChargingStation(arr2 pos, double Angle, float s, std::string Colortext, std::string Label):obstacle(pos,s,Colortext,Label)
    {
        angle=Angle;assert(angle>=0);assert(angle<360);
        blocked_areas=calculateBlockedAreas();
    }
    double getangle(){return angle;}
    void setangle(double Angle){angle=Angle;}
    double getwidth(){return width;}
    void setwidth(double w){width=w;}
    double getheight(){return height;}
    void setheight(double h){height=h;}

    //for a given Charging Station, this calculates the three rectangular areas blocking the way
    //in combination, they make up the horseshoe-like shape
    std::array<turned_rectangle,3> calculateBlockedAreas()
    {
        //setup
        std::array<turned_rectangle,3> blocked_areas;
        arr2 center = getpos();
        double angle = getangle();
        double height=getheight();
        double width=getwidth();

        //rectangle 1: below center for angle=0°=> remember here that in svg, y axis goes downward!
        blocked_areas[0].xpos=center[0] - sin(angle/180*3.1415926) * height;
        blocked_areas[0].ypos=center[1] + cos(angle/180*3.1415926)*height;
        blocked_areas[0].angle=angle;
        blocked_areas[0].width=width;
        blocked_areas[0].height=1;

        //rectangle 2: left of center for angle=0°
        blocked_areas[1].width=1;
        blocked_areas[1].xpos=center[0] - cos(angle/180*3.1415926)*(0.75*width);
        blocked_areas[1].ypos=center[1] - sin(angle/180*3.1415926)*(0.75*width);
        blocked_areas[1].angle=angle;
        blocked_areas[1].height=height+2;

        //rectangle 3: on top of center for angle=0°
        blocked_areas[2].xpos=center[0] + sin(angle/180*3.1415926) * height;
        blocked_areas[2].ypos=center[1] - cos(angle/180*3.1415926)*height;
        blocked_areas[2].angle=angle;
        blocked_areas[2].width=width;
        blocked_areas[2].height=1;
        
        //debugging: print it all out
        //std::cout<<"rectangle above center: \n";
        //std::cout<<"pos: ("<<blocked_areas[2].xpos<<","<<blocked_areas[2].ypos<<")\n";
        //std::cout<<"size: ("<<blocked_areas[2].width<<","<<blocked_areas[2].height<<")\n";
        return blocked_areas;
    }
    std::array<turned_rectangle,3> getblocked_areas()
    {return blocked_areas;}
};

class map {
    unsigned tolerance; //space around each obstacle for the robot to avoid
    arr2 size; //length + width; this parameter is also used for rectangles
    int scale; //width of the actual factory floor for scale
    std::vector<rectangle> RectangleVector;
    std::vector<circle> CircleVector;
    std::vector<ChargingStation> ChargerVector;
    //upper left corner is on (0,0) by default
    public:
    //initialize (only one constructor)
    map(arr2 Size, unsigned s, unsigned Tolerance){
        size[0]=Size[0]*s;
        size[1]=Size[1]*s;
        assert(s>0);
        scale = s;
        tolerance = Tolerance*s;
    }

    //getters and setters
    arr2 getsize(){return {size[0]/scale,size[1]/scale};}
    unsigned gettolerance(){return tolerance/scale;}
    unsigned getscale(){return scale;}
    std::vector<rectangle> getRectangleVector(){return RectangleVector;}
    std::vector<circle> getCircleVector(){return CircleVector;}
    std::vector<ChargingStation> getChargerVector(){return ChargerVector;}
    void setscale(unsigned s){scale=s;}
    void setsize(arr2 Size){size={Size[0]*scale,Size[1]*scale};}
    void settolerance(unsigned tol){tolerance=tol*scale;}

    //method that adds an obstacle to our map
    void addobstacle(const obstacle& obs)
    {
        //Differentiate between rectangle and circle
        if (const rectangle* rect = dynamic_cast<const rectangle*>(&obs))
        {
            // Handle rectangle object
            RectangleVector.push_back(static_cast<const rectangle&>(obs));//static cast for conversion to rectangle, push_back = append in Python
        }
        else if (const circle* circ = dynamic_cast<const circle*>(&obs))
        {
            // Handle circle object
            CircleVector.push_back(static_cast<const circle&>(obs));
        }
        else{std::cerr<<"Error during map creation: Added non-viable obstacle.";}
    }



    //method that adds Charging Station to vector
    void addChargingStation(const ChargingStation& newcharger)
    {ChargerVector.push_back(newcharger);}

    //method that prints out map size and list of obstacles => useful for debugging
    void maplist(){
        std::cout<< "map: width "<<size[0]/scale <<", height "<<size[1]/scale<<", tolerance "<< tolerance/scale<<std::endl;
        for (auto& rect : RectangleVector) {
            std::cout << "rectangle: width " << rect.getsize()[0]
              << ", height " << rect.getsize()[1]
              << " at position (" << rect.getpos()[0] << ", " << rect.getpos()[1] << ")"
              << std::endl;}
        for (auto& circ : CircleVector) {
            std::cout << "circle: radius " << circ.getradius()
              << " at position (" << circ.getpos()[0] << ", " << circ.getpos()[1] << ")"
              << std::endl;}
        for (auto& charger : ChargerVector){
            std::cout<< "charging station: angle "<<charger.getangle()<<" at position ("<<charger.getpos()[0]<<", "
            <<charger.getpos()[1]<<")\n";
        }
    }
    
    //create an svg file of the map
    int createsvgmap() {
        //todo: make this function more modular
        if(scale<10){std::cerr << "\033[33mWarning: scale is smaller than 10. This can lead to errors when creating .svg images. Set scale to at least 10 in the map and all obstacles to get correct plots.\033[0m" << std::endl;}
        //open or create .svg file
        std::ofstream svgFile("map.svg");
        if (!svgFile) {
            std::cout << "Error opening file." << std::endl;
            return 1;
        }
        //parameters:
        unsigned fontsize=static_cast<unsigned>(12*scale/10);//Font size scaled to fit resolution chosen by user
        unsigned numGridLines = 9; // Number of grid lines in both coordinate directions
        unsigned strokewidth=static_cast<unsigned>(1*scale/10);//Line thickness. Note that the line is thicker for the frame

        //for modifying the lines writing into the svg file, I highly recommend ChatGPT
        //draw map
        svgFile << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << size[0] << "\" height=\"" << size[1] << "\" version=\"1.1\">\n";
        
        //create frame
        svgFile << "<rect x=\"0\" y=\"0\" width=\"" << size[0] << "\" height=\"" << size[1] << "\" fill=\"none\" stroke=\"black\" stroke-width=\""<<3*strokewidth<<"\" />\n";
        //create numbers for the frame
        arr2 realsize=getsize();
        svgFile << "<text x=\"" << size[0] - 10 << "\" y=\""<<static_cast<unsigned>(20*scale/10)<<"\" font-size=\""<<fontsize<<"\" text-anchor=\"end\">" << realsize[0] << "</text>\n";
        svgFile << "<text x=\"10\" y=\"" << size[1] - 10 << "\" font-size=\""<<fontsize<<"\">" << realsize[1] << "</text>\n";

        //create grid
        for (int i = 1; i <= numGridLines; i++) {
            float y = size[1] * i / (numGridLines + 1);
            svgFile << "<line x1=\"0\" y1=\"" << y << "\" x2=\"" << size[0] << "\" y2=\"" << y << "\" stroke=\"gray\" stroke-width=\""<<strokewidth<<"\" />\n";//horizontal line
            svgFile << "<text x=\"5\" y=\"" << y - 5 << "\" font-size=\""<<fontsize<<"\" fill=\"gray\">" << realsize[1] * y / size[1] << "</text>\n";//numbers horizontal
            float x = size[0] * i / (numGridLines + 1);
            svgFile << "<line x1=\"" << x << "\" y1=\"0\" x2=\"" << x << "\" y2=\"" << size[1] << "\" stroke=\"gray\" stroke-width=\""<<strokewidth<<"\" />\n";//vertical line
            svgFile << "<text x=\"" << x - 5 << "\" y=\""<<static_cast<unsigned>(20*scale/10)<<"\" font-size=\""<<fontsize<<"\" fill=\"gray\" text-anchor=\"end\">" << realsize[0] * x / size[0] << "</text>\n";//numbers vertical
        }

        // draw rectangles into the SVG file with the specified color
        for (auto& rect : RectangleVector) {
            //get color of rectangle, transform into suitable format
            RGBColor rgbval = rect.getcolor();
            std::string svgColor = "rgb(" + std::to_string(rgbval.red) + "," +
                            std::to_string(rgbval.green) + "," +
                            std::to_string(rgbval.blue) + ")";
            //scale up pos and size of rectangle
            arr2 pos = {rect.getpos()[0] * rect.getscale(), rect.getpos()[1] * rect.getscale()};
            arr2 size = {rect.getsize()[0] * rect.getscale(), rect.getsize()[1] * rect.getscale()};
            //draw the rectangles
            svgFile << "<rect x=\"" << pos[0] - tolerance << "\" y=\"" << pos[1] - tolerance
                    << "\" width=\"" << size[0] + 2 * tolerance << "\" height=\"" << size[1] + 2 * tolerance
                    << "\" rx=\"" << tolerance << "\" ry=\"" << tolerance
                    << "\" fill=\"rgba(192,192,192,0.6)\" />\n";//rectangle with tolerance
            svgFile << "<rect x=\"" << pos[0] << "\" y=\"" << pos[1] << "\" width=\"" << size[0]
                    << "\" height=\"" << size[1] << "\" fill=\"" << svgColor << "\" />\n"; // rectangle without tolerance
            
            //write labels; if high rectangle, rotate text by 90°; if dark color, make label white
            std::string textColor;
            if(rgbval.red+rgbval.green+rgbval.blue<400){textColor="rgb(255,255,255)";}else{textColor="rgb(0,0,0)";}
            if (size[0] >= size[1]) {
                svgFile << "<text x=\"" << pos[0] + size[0] / 2 << "\" y=\"" << pos[1] + size[1] / 2 << "\" fill=\"" << textColor << "\" text-anchor=\"middle\" alignment-baseline=\"middle\">" << rect.getlabel() << "</text>\n";
            } else {
                svgFile << "<text x=\"" << pos[0] + size[0] / 2 << "\" y=\"" << pos[1] + size[1] / 2 << "\" fill=\"" << textColor << "\" text-anchor=\"middle\" alignment-baseline=\"middle\" transform=\"rotate(270 " << pos[0] + size[0] / 2 << "," << pos[1] + size[1] / 2 << ")\">" << rect.getlabel() << "</text>\n";
            }
        }

        //draw circles(=obstacles)
        for (auto& circ: CircleVector){
            RGBColor rgbval = circ.getcolor();
            std::string svgColor = "rgb(" + std::to_string(rgbval.red) + "," +
                            std::to_string(rgbval.green) + "," +
                            std::to_string(rgbval.blue) + ")";
            arr2 pos={circ.getpos()[0]*circ.getscale(),circ.getpos()[1]*circ.getscale()};
            unsigned rad=circ.getradius()*circ.getscale();
            svgFile << "<circle cx=\""<<pos[0]<<"\" cy=\""<<pos[1]<<"\" r=\""<<rad+tolerance<<"\" fill=\"rgba(192,192,192,0.6)\" />\n";//circle with tolerance
            svgFile << "<circle cx=\""<<pos[0]<<"\" cy=\""<<pos[1]<<"\" r=\""<<rad<<"\" fill=\"" << svgColor << "\" />\n";//circle without tolerance

            //labels:
            //std::cout<<"label: "<<circ.getlabel()<<std::endl;
            std::string textColor;
            if(rgbval.red+rgbval.green+rgbval.blue<400){textColor="rgb(255,255,255)";}else{textColor="rgb(0,0,0)";}
            svgFile << "<text x=\"" << pos[0] << "\" y=\"" << pos[1] << "\" fill=\"" << textColor << "\" text-anchor=\"middle\" alignment-baseline=\"middle\">" << circ.getlabel() << "</text>\n";
        }

        //draw Charging Stations
        for(auto& charger: ChargerVector){
            RGBColor rgbval = charger.getcolor();
            std::string svgColor = "rgb(" + std::to_string(rgbval.red) + "," +
                            std::to_string(rgbval.green) + "," +
                            std::to_string(rgbval.blue) + ")";
            double width=charger.getwidth()*charger.getscale();
            double height=charger.getheight()*charger.getscale();
            std::array<double,2> pos={charger.getpos()[0]*charger.getscale()-width/2,charger.getpos()[1]*charger.getscale()-height/2};
            double angle=charger.getangle();
            //svgFile << "<rect x=\"" << pos[0] << "\" y=\"" << pos[1] << "\" width=\"" << width
            //    << "\" height=\"" << height << "\" fill=\"" << svgColor
            //    << "\" transform=\"rotate(" << angle << " " << pos[0] + width / 2 << " " << pos[1] + height / 2 << ")\" />\n";
        //draw blocked areas
        for(auto& block: charger.getblocked_areas())
        {
            width=block.width*charger.getscale();
            height=block.height*charger.getscale();
            pos={block.xpos*charger.getscale()-width/2,block.ypos*charger.getscale()-height/2};
            angle=charger.getangle();
            svgFile << "<rect x=\"" << pos[0] << "\" y=\"" << pos[1] << "\" width=\"" << width
                << "\" height=\"" << height << "\" fill=\"" << svgColor
                << "\" transform=\"rotate(" << angle << " " << pos[0] + width / 2 << " " << pos[1] + height / 2 << ")\" />\n";
        }
        }

        //finish svg file
        svgFile << "</svg>\n";
        svgFile.close();

        //std::cout << "SVG file created successfully." << std::endl;

        return 0;
    }
    
};

//class for AGV; This one is supposed to be changed for every query.
class query {
private:
    std::array<double, 2> startPoint, endPoint;
    double startAngle, endAngle;
    //* These three variables currently do nothing in the code.
    std::array<bool,2> Charging_Station;
    double obs_avoid_weight; //facV in Matlab code. Determines how much violation score contributes to fitness value.
    unsigned coll_avoid_resolution; //during collision avoidance, path is diiscretisized into points. This determines how many points the path is divided into.

public:
    // Constructors
    query(const std::array<double, 2>& startP, const std::array<double, 2>& endP,
        double sAngle, double eAngle, std::array<bool,2> Charge, double facV, unsigned nrEl)
        : startPoint(startP), endPoint(endP), startAngle(sAngle), endAngle(eAngle), Charging_Station(Charge), obs_avoid_weight(facV), coll_avoid_resolution(nrEl) {}
    query()
        : startPoint({0.0, 0.0}), endPoint({0.0, 0.0}), startAngle(0.0), endAngle(0.0), Charging_Station({0,1}), obs_avoid_weight(0.05), coll_avoid_resolution(300) {}
    //query()=default;
    // Getters
    std::array<double, 2> getStartPoint() const {return startPoint;}
    std::array<double, 2> getEndPoint() const {return endPoint;}
    double getStartAngle() const {return startAngle;}
    double getEndAngle() const {return endAngle;}
    double getobs_avoid_weight() {return obs_avoid_weight;}
    unsigned getcoll_avoid_resolution(){return coll_avoid_resolution;}
    std::array<bool,2> getCharging_Station() const{return Charging_Station;}

    // Setters
    void setStartPoint(const std::array<double, 2>& startP) {
        startPoint = startP;
        assert(startPoint[0]>=0 && startPoint[1]>=0);}
    //set Charging Station as Start Point
    void setStart(ChargingStation charger){
        arr2 pos=charger.getpos();
        startPoint = {static_cast<double>(pos[0]), static_cast<double>(pos[1])};
        startAngle=charger.getangle();
    }

    //set Charging Station as End Point
    void setEnd(ChargingStation charger){
        arr2 pos=charger.getpos();
        endPoint = {static_cast<double>(pos[0]), static_cast<double>(pos[1])};
        endAngle=charger.getangle();
    }

    void setEndPoint(const std::array<double, 2>& endP) {
        endPoint = endP;
        assert(endPoint[0]>=0 && endPoint[1]>=0);}

    void setStartAngle(double sAngle) {
        startAngle = sAngle;if(startAngle<0){startAngle+=360;}}

    void setEndAngle(double eAngle) {
        endAngle = eAngle;if(endAngle<0){endAngle+=360;}}

    void setStart(std::array<double,2> point,double angle)
    {setStartPoint(point);setStartAngle(angle);}

    void setEnd(std::array<double,2> point, double angle)
    {setEndPoint(point);setEndAngle(angle);}

    void setobs_avoid_weight(double facV){obs_avoid_weight=facV;assert(facV>0);}
    void setcoll_avoid_resolution(unsigned nrEl){coll_avoid_resolution=nrEl;}

    void setCharging_Station(std::array<bool,2> charge) {Charging_Station=charge; assert(Charging_Station[0] == false or Charging_Station[1] == false);}
};

//function that prints out the model struct => useful for debugging
void printModel(const model_struct& Model) {
    std::cout << "AisChargingStation: " << Model.AisChargingStation << std::endl;
    std::cout << "BisChargingStation: " << Model.BisChargingStation << std::endl;
    std::cout << "SA: " << Model.SA << std::endl;
    std::cout << "EA: " << Model.EA << std::endl;
    std::cout << "tolerance: " << Model.tolerance << std::endl;
    std::cout << "xs: " << Model.xs << std::endl;
    std::cout << "ys: " << Model.ys << std::endl;
    std::cout << "xt: " << Model.xt << std::endl;
    std::cout << "yt: " << Model.yt << std::endl;

    std::cout << "xobs_rect: ";
    for (const auto& x : Model.xobs_rect) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    std::cout << "yobs_rect: ";
    for (const auto& y : Model.yobs_rect) {
        std::cout << y << " ";
    }
    std::cout << std::endl;

    std::cout << "xobs_circ: ";
    for (const auto& x : Model.xobs_circ) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    std::cout << "yobs_circ: ";
    for (const auto& y : Model.yobs_circ) {
        std::cout << y << " ";
    }
    std::cout << std::endl;

    std::cout << "xwidth: ";
    for (const auto& xw : Model.xwidth) {
        std::cout << xw << " ";
    }
    std::cout << std::endl;

    std::cout << "ywidth: ";
    for (const auto& yw : Model.ywidth) {
        std::cout << yw << " ";
    }
    std::cout << std::endl;

    std::cout << "robs_circ: ";
    for (const auto& yw : Model.robs_circ) {
        std::cout << yw << " ";
    }
    std::cout << std::endl;

    std::cout << "n: " << Model.n << std::endl;

    std::cout << "lb: "<<Model.lb<<std::endl;
    std::cout << "ub: "<<Model.ub<<std::endl;
    
    //std::cout << "facV: "<<Model.facV<< std::endl;
    std::cout<<"nrEl: "<<Model.nrEl<<std::endl;

    std::cout << std::endl;
}

//function that transforms map object into a model the algorithm can work with.
model_struct maptomodel(map Map, query Query)
{
    //setup
    std::vector<rectangle> RectangleVector=Map.getRectangleVector();
    std::vector<circle> CircleVector=Map.getCircleVector();
    unsigned n_rect = RectangleVector.size();
    unsigned n_circ = CircleVector.size();
    struct model_struct Model;
    struct model_struct *p_Model = &Model;

    //While this step should theoretically be redundant, it avoids a segfault, so it's good to have
    Model.xobs_rect.resize(n_rect); 
    Model.yobs_rect.resize(n_rect);
    Model.xwidth.resize(n_rect);
    Model.ywidth.resize(n_rect);

    //start filling up model struct
    Model.n=n_rect+n_circ;
    assert(Model.xobs_rect.size()==n_rect);

    //Rectangles
    double tolerance = Map.gettolerance();
    std::vector<arr2> rect_pos(n_rect);
    std::vector<arr2> rect_center_pos(n_rect);
    std::vector<arr2> rect_size(n_rect);
    arr2 mapsize = Map.getsize();
    
    for(int i=0;i<n_rect;i++) //extract pos and size of rectangles, get rectangles into right shape
    {
        rect_pos[i] = RectangleVector[i].getpos();
        rect_size[i] = RectangleVector[i].getsize();
        arr2 rect_size_05={rect_size[i][0] /2,rect_size[i][1] /2};
        rect_center_pos[i][0] = rect_pos[i][0] + rect_size_05[0];
        rect_center_pos[i][1] = rect_pos[i][1] + rect_size_05[1];

        p_Model->xobs_rect[i] = rect_center_pos[i][0]; 
        p_Model->yobs_rect[i] = rect_center_pos[i][1]; 
        p_Model->xwidth[i] = rect_size[i][0];
        Model.ywidth[i] = rect_size[i][1];
    }
    assert(Model.xobs_rect[0]==rect_center_pos[0][0]);
    assert(Model.ywidth[1]==rect_size[1][1]);

    //set obstacle to account for non-squared map size
    //!Since the algorithm only works with quadratic maps, we add an obstacle to make sure it uses the actual map size
    //It uses the larger size as the square length, then adds an obs at the end
    //If the map is actually square-shaped, the code jumps (yes, with a goto-statement!) to the end of this section.
    p_Model->lb=0;
    double max_value = std::max(mapsize[0], mapsize[1]);
    Model.ub=max_value;
    int index = (mapsize[0] == max_value) ? 0 : 1;
    double xwidth, ywidth, xpos, ypos;
    
    if (mapsize[0]==mapsize[1]) { // Both equally as long
        Model.ub = mapsize[0];
        goto skip_non_quad_compensation;
    } else if (index == 0) { // If size in x direction larger
        // Shorten Model in y direction
        //std::cout<<"Horizontal format!"<<std::endl;
        xwidth = mapsize[0];
        ywidth = mapsize[0] - mapsize[1]; // Adjust for tolerance here
        xpos = 0.5 * mapsize[0];
        ypos = mapsize[1] + 0.5 * ywidth;
    } else { // If size in y direction larger
        // Shorten Model in x direction
        //std::cout<<"Vertical format!"<<std::endl;
        ywidth = mapsize[1] + 2 * Model.tolerance;
        xwidth = mapsize[1] - mapsize[0] + 2 * Model.tolerance; // Adjust for tolerance here
        xpos = mapsize[0] + 0.5 * xwidth + Model.tolerance;
        ypos = 0.5 * mapsize[1] + Model.tolerance;}

    //Model.n += 1; // Adjust number of obstacles
    Model.xobs_rect.push_back(xpos);
    Model.yobs_rect.push_back(ypos);
    Model.xwidth.push_back(xwidth);
    Model.ywidth.push_back(ywidth);
    Model.n += 1; // Adjust number of obstacles
    assert(Model.xobs_rect.size()==n_rect+1);
    //std::cout<<"xpos: "<<xpos<<" ypos: "<<ypos<<" xwidth: "<<xwidth<<" ywidth: "<<ywidth<<std::endl;
    skip_non_quad_compensation:

    //Circles
    for(int i=0;i<n_circ;i++)
    {
        arr2 pos=CircleVector[i].getpos();
        Model.xobs_circ.push_back( pos[0] );
        Model.yobs_circ.push_back( pos[1] );
        double radius = CircleVector[i].getradius();
        Model.robs_circ.push_back(radius);
    }

    //Values of start, end etc. => extract from query class
    std::array<double,2> point_temp= Query.getStartPoint();
    Model.xs=point_temp[0]; Model.ys=point_temp[1];
    point_temp= Query.getEndPoint();
    Model.xt=point_temp[0]; Model.yt=point_temp[1];
    Model.SA=Query.getStartAngle();
    Model.EA=Query.getEndAngle();
    Model.tolerance=Map.gettolerance();
    std::array<bool,2> charge =Query.getCharging_Station();
    Model.AisChargingStation=charge[0];
    Model.BisChargingStation=charge[1];
    //Model.facV=Query.getobs_avoid_weight(); //outdated: We now write algorithm data into the config file
    Model.nrEl=Query.getcoll_avoid_resolution();

    //Charging Stations
    for(auto& charger: Map.getChargerVector())
    {
        Model.blocked_area_down.push_back(charger.getblocked_areas()[0]);
        Model.blocked_area_left.push_back(charger.getblocked_areas()[1]);
        Model.blocked_area_up.push_back(charger.getblocked_areas()[2]);
    }

    return Model;
}

/*
int main()
{
    //create map
    unsigned scale = 10; //How many pixels (svg units) are in one meter? 
    map mymap(arr2 {100,50},scale,2);
    mymap.addobstacle(rectangle(arr2{10,10}, arr2{30,5}, scale));
    mymap.addobstacle(rectangle(arr2{50,10}, arr2{40,5}, scale));
    mymap.addobstacle(rectangle(arr2{10,25}, arr2{5,20}, scale));
    mymap.addobstacle(rectangle(arr2{20,25}, arr2{5,20}, scale));
    mymap.addobstacle(rectangle(arr2{30,25}, arr2{5,20}, scale));
    mymap.addobstacle(rectangle(arr2{40,25}, arr2{5,20}, scale));
    mymap.addobstacle(rectangle(arr2{50,25}, arr2{5,20}, scale));
    mymap.addobstacle(circle(arr2{75,35}, 12, scale));
    mymap.maplist();
    mymap.createsvgmap();

    //set AGV start and end point with start and end angles.
    query query1;
    query1.setStart({2.0, 2.0}, 270); //start point and angle at start point
    query1.setEnd({88.0, 26.0}, 270); //same for end point & angle
    query1.setCharging_Station({0,1}); //optional; sets start or end as a Charging Station
    query1.setcoll_avoid_resolution(400); //optional: sets resolution of collision avoidance
    query1.setobs_avoid_weight(0.01); //optional: sets how much the algorithm is punished for collision

    model Model= maptomodel(mymap, query1);
    printModel(Model);
}
*/
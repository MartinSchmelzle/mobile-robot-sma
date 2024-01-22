#include <cmath>
#include "slimemould_func_prototypes.h"
#include "mapbuilder.cpp"
#include "processpath.cpp"
#include<string>

int main()
{
  //mapbuilder
  unsigned scale = 100; //How many pixels (svg units) are in one meter? 
  map mymap(arr2 {100,50},scale,1);
  mymap.addobstacle(rectangle(arr2{10,10}, arr2{30,5}, scale, "white", "CNC machines"));
  mymap.addobstacle(rectangle(arr2{50,10}, arr2{40,5}, scale, "teal", "3D printers"));
  mymap.addobstacle(rectangle(arr2{10,25}, arr2{5,20}, scale, "red", "Casting"));
  mymap.addobstacle(rectangle(arr2{30,25}, arr2{5,20}, scale));
  mymap.addobstacle(rectangle(arr2{40,25}, arr2{5,20}, scale));
  mymap.addobstacle(rectangle(arr2{50,25}, arr2{5,20}, scale));
  mymap.addobstacle(circle(arr2{75,35}, 12, scale, "maroon", "Nuclear Reactor"));
  //add charging stations
  ChargingStation Station_1(arr2{2,2},90,scale,"orange","station_1");
  ChargingStation Station_2(arr2{5,2},90,scale,"orange","station_2");
  ChargingStation Station_3(arr2{20,23},135,scale,"orange","station_3");
  mymap.addChargingStation(Station_1);
  mymap.addChargingStation(Station_2);
  mymap.addChargingStation(Station_3);

  //print out map, create svg
  //mymap.maplist();
  mymap.createsvgmap();

  //set AGV start and end point with start and end angles.
  query query1;
  //query1.setStart({2.0, 2.0}, 90); //start point and angle at start point
  query1.setStart(Station_1);
  //query1.setEnd(Station_3);
  query1.setEnd({88.0, 26.0}, 270); //same for end point & angle
  query1.setCharging_Station({0,1}); //optional; sets start or end as a Charging Station
  query1.setcoll_avoid_resolution(400); //optional: sets resolution of collision avoidance
  query1.setobs_avoid_weight(0.01); //optional: sets how much the algorithm is punished for collision

  model_struct model= maptomodel(mymap, query1);
  //printModel(model);

  //declare empty struct for path
  path_struct pathfromalg;
  //run algorithm
  bool k=main_alg(&model,&pathfromalg);//call main function
  //if(k==0){std::cout<<"path computed successfully!"<<std::endl;}else{std::cout<<"error with path calculation!"<<std::endl;}

  //path processing
    //printPathStruct(pathfromalg);
    drawpathintosvg(pathfromalg, scale);
    double distance_betw_points = 1.6;
    std::vector<point> pointcloud=pathtopointcloud(pathfromalg, distance_betw_points);
    drawpointcloud_svg(pointcloud, scale); //todo This function should be reviewed before handing it in => prone to bugs
    pathtocsv(pathfromalg);
    std::cout<<"Path processed successfully!"<<std::endl;
    return 0;
}
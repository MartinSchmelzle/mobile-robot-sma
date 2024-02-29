#include <cmath>
#include <string>
#include "../slimemould_func_prototypes.h"
#include "../src/mapbuilder.cpp"
#include<string>
#include <chrono>

int main()
{
  //mapbuilder
  unsigned scale = 10; //How many pixels (svg units) are in one meter? 
  map mymap(arr2 {100,50},scale,1);
  //mymap.addobstacle(rectangle(arr2{0,-15}, arr2{100,16}, scale));
  mymap.addobstacle(rectangle(arr2{10,10}, arr2{30,5}, scale, "white", "CNC machines"));
  mymap.addobstacle(rectangle(arr2{50,10}, arr2{43,5}, scale, "teal", "3D printers"));
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

  //measure performance from here bcs this part changes every query
  auto start = std::chrono::high_resolution_clock::now();
  //set AGV start and end point with start and end angles.
  query query1;
  //query1.setStart({2.0, 2.0}, 90); //start point and angle at start point
  query1.setStart(Station_1);
  //query1.setEnd(Station_3);
  query1.setEnd({88.0, 45.0}, 270); //same for end point & angle

  model_struct model= maptomodel(mymap, query1);
  //printModel(model);

  //declare empty struct for path
  SMAsol_struct* sol = new SMAsol_struct;
  //run algorithm
  bool k=mr_sma_alg(&model,sol);//call main function
  //measure performance until here bcs the plot is not executed every time
  auto end = std::chrono::high_resolution_clock::now(); 
  std::chrono::duration<double> duration = end - start;
  if(k==0){
  //path processing
    //printPathStruct(sol);
    drawpathintosvg(sol, scale);
    double distance_betw_points = 1.6;
    std::vector<point> pointcloud=pathtopointcloud(sol, distance_betw_points, model);
    drawpointcloud_svg(pointcloud, scale);
    pathtocsv(sol);
    //std::cout<<"Path processed successfully!"<<std::endl;
  } 
  delete sol;
  std::cout << "Algorithm performance: " << duration.count() << " s\n";
    return 0;
}
#include <cmath>
#include <string>
#include "../slimemould_func_prototypes.h"
#include "../src/mapbuilder.cpp"
#include <string>
#include <chrono>
#include <fstream>
#include <sstream>

int main()
{
  //mapbuilder
  unsigned scale = 10; //How many pixels (svg units) are in one meter? 
  map mymap(arr2 {132,171},scale,1);

    // Open the CSV file
    std::ifstream file("examples/example-2-obs.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    std::vector<std::array<int, 4>> data; // Vector of arrays to store the values

    std::string line;
    while (std::getline(file, line)) { // Read each line from the file
        std::istringstream iss(line);
        std::array<int, 4> values; // Array to store the values of each line

        for (int i = 0; i < 4; ++i) {
            std::string value_str;
            if (!std::getline(iss, value_str, ',')) {
                std::cerr << "Error parsing line: " << line << std::endl;
                return 1;
            }
            
            // Convert string to integer
            char* end;
            int value = std::strtol(value_str.c_str(), &end, 10);
            if (*end != '\0') {
                std::cerr << "Error converting string to integer: " << value_str << std::endl;
                return 1;
            }
            values[i] = value;
        }

        data.push_back(values); // Add the array to the vector
    }

  for(const auto& arr: data){
    mymap.addobstacle(rectangle(arr2{arr[0],arr[1]}, arr2{arr[2],arr[3]}, scale));
  }

  //add charging stations
  ChargingStation Station_1(arr2{2,2},90,scale,"orange","station_1");
  ChargingStation Station_2(arr2{5,2},90,scale,"orange","station_2");
  mymap.addChargingStation(Station_1);
  mymap.addChargingStation(Station_2);

  //print out map, create svg
  //mymap.maplist();
  mymap.createsvgmap();

  //measure performance from here bcs this part changes every query
  auto start = std::chrono::high_resolution_clock::now();
  //set AGV start and end point with start and end angles.
  query query1;
  query1.setStart({2,2}, 90); //start point and angle at start point
  query1.setEnd({66,17.1}, 180); //same for end point & angle

  model_struct model= maptomodel(mymap, query1);
  printModel(model);

  //declare empty struct for path
  path_struct pathfromalg;
  //run algorithm
  bool k=mr_sma_alg(&model,&pathfromalg);//call main function
  //measure performance until here bcs the plot is not executed every time
  auto end = std::chrono::high_resolution_clock::now(); 
  std::chrono::duration<double> duration = end - start;
  if(k==0){
  //path processing
    //printPathStruct(pathfromalg);
    drawpathintosvg(pathfromalg, scale);
    double distance_betw_points = 1.6;
    std::vector<point> pointcloud=pathtopointcloud(pathfromalg, distance_betw_points, model);
    drawpointcloud_svg(pointcloud, scale);
    pathtocsv(pathfromalg);
    //std::cout<<"Path processed successfully!"<<std::endl;
  } 
  std::cout << "Algorithm performance: " << duration.count() << " s\n";
    return 0;
}
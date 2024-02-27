A path planning algorithm designed for mobile robots that optimizes energy efficiency.
![image](https://github.com/MartinSchmelzle/mobile-robot-sma/assets/120244663/d0be3993-e10a-46ab-acf8-cb04fd9e269f)

It provides a simple way to model a map layout, compute an optimal path and process the results.
This file describes how to get it to run and what it can do for you.

# Installation
1. Install compiler and build system: g++, CMake and Ninja / Make. The code comes with a CMakeLists.txt.
2. git clone the repository into a directory of your choosing.
3. It is recommended to delete CMakeCache.txt before running the program for the first time.
4. To run the code, I recommend using the terminal: After `cd` ing into the build directory, type `cmake ..; cmake --build; ninja` to build the project
5. An executable `slime_alg_exe.exe` is generated in the build directory. Run it. The example `CALL-ALG.cpp` should be executed automatically.

# Usage
See the script `CALL-ALG.cpp` as an example.
As for includes, you will need `#include "slimemould_func_prototypes.h"`, `#include "mapbuilder.cpp"` and `#include "processpath.cpp"`.
The first step is map creation. Create a rectangular map `map mymap(arr2 {100,50},scale,tol);` with a size of your choosing (in this example, 100x50).
If your map is not rectangular, that is not a problem; in this case, create a rectangular map with the maximum dimensions and create obstacles later to emulate your desired map shape.
Whether you prefer metric or imperial units is entirely up to you. The variable `scale` defines the precision of the map. It should be set to a value 10 or higher.
`tol` sets a space around each obstacle which the robot may not enter. This helps with collision avoidance. Please set this value in accordance with your robots geometry.

After these values are decided, you can start adding obstacles using `map.addobstacle`. There are currently two kinds of obstacles supported: Rectangles and circles.
Rectangles have two values for position (referring to the upper left corner) and two values for size (referring to length and width).
You can optionally add a color and a label to each obstacle: `mymap.addobstacle(rectangle(pos, size, scale, color, label));`
For circles, the input variables are almost the same: `size` only has one entry instead of two (the radius). Additionally, the position here refers to the circle's center.
Future versions may include triangles, rotated rectangles and / or boolean operations.

After adding obstacles, you likely want to add charging stations using `map.addChargingStation`. 
They have a fixed place and orientation, moreover they include a U-shaped obstacle, meaning that robots can only enter and exit from one side.
Charging stations have a position (x and y), an angle in degrees (relative to the positive x axis), a scale and, optionally, a color and a label.

After adding obstacles and charging stations, you may want to visalize and check your map. The function `map.createsvgmap()` creates an svg file called map.svg based on the map you built.

As a next step, a query is made: Using `query.setStart` and `query.setEnd`, start and end points are set. Instead of a 2D position, you can also set a Charging Station as an input.
Before the algorithm is called, the map is turned into a model using `maptomodel(map,query)`, a path is declared and the maximum amount of iterations (we recommend 150-300 iterations) is set.
The algorithm itself is called using `  bool k=main_alg(&model,&path,max_iter);`.

The algorithm will attempt to find a solution to the query five times. If it cannot find a solution after the fifth attempt, it will return -1, if it is successfull, it will return 0.

After the algorithm runs, the path can be processed in different ways:
- `printPathStruct(path);` prints the resulting path into the terminal.
- `drawpathintosvg(pathfromalg, scale);` draws the resulting path into the map's svg file.
- `double dist = 1.6; std::vector<point> pointcloud=pathtopointcloud(pathfromalg, dist, model);` generates a 2D point cloud from the path.
- `drawpointcloud_svg(pointcloud, scale);` draws the resulting point cloud into the svg file.
- `pathtocsv(pathfromalg);` writes the path struct into a csv file.

# Community
If you have difficulty using the algorithm, would like to suggest improvements or would like to work together with us, feel free to message any of us:
https://de.linkedin.com/in/martin-schmelzle






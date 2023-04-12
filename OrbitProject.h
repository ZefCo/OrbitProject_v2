#pragma once
#include <tuple>
#include <map>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include "json/json.h"
#include "RungeKutta.h"

// https://stackoverflow.com/questions/72604248/cannot-overload-functions-distinguished-by-return-type-alone-but-the-function

// const int X = 0, Y = 1, RN = 0, VN = 1;

namespace fs = std::filesystem;


// Imports the data from a JSON file. The cool part of this is that the JSON file can be easily modified with more or less masses
// to simulate different systems. The objects are not hard coded in, expcept for how the object titled "Sun" is treated.
std::tuple<std::vector<std::string>, std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> create_system(fs::path filename);


// Because the acceleration is -G*sum(Mj [vec]) I find it easier to just multiply all the
// masses by -G, meaning that the only multiplication that needs to be done over and over is
// the Mj [vec] part, which will change every iteration anyway.
std::vector<double> gmass(std::vector<double> masses, double G);


// Calls the RK4 function, finds the time step, and adapts it as needed. The only object hard coded in (which is actually hard coded into the RungeKutta.cpp file)
// is one titled "Sun". An object title "Sun" is set to have 0 acc, 0 vel, and 0 position: it doesn't move. This is just to save on some computation time
// and looking for the adaptive time step. In earlier iterations I did let the sun move, but it just took up resrouces to do nearly nothing. It would be fun to see
// the sun wobble from Jupiter, but it wasn't enough to really see on a graph.
void orbits(double time_steps, double h, std::vector<std::string> names, std::vector<double> masses, std::vector<std::vector<double>> rn_data, std::vector<std::vector<double>> vn_data, double delta);


// Finds the euclidean distance between ri and rj by row. Any value of 0 (or very close to it) is set to 16 (2**4). This is used when looking for rho.
std::vector<double> euc_dist(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);


// Finds the worst possible value for rho and adapts the step size to it.
double wrost_rho(std::vector<double> a, double delta, double h);
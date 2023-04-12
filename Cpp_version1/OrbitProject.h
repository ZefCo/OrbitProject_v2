#pragma once
// #include <tuple>
#include <map>
#include <vector>
#include <filesystem>
#include <fstream>
#include "json/json.h"
#include "Wanderer.h"

// https://stackoverflow.com/questions/72604248/cannot-overload-functions-distinguished-by-return-type-alone-but-the-function


namespace fs = std::filesystem;


// Loads the JSON file with all the orbital body data and returns it as a tuple for easy unpacking.
std::tuple<std::vector<double>,
           std::vector<double>, 
           std::vector<double>, 
           std::vector<double>, 
           std::vector<double>, 
           std::vector<std::string>> loadJSON_old(fs::path filename);


// Loads the JSON data into the wanderer class for each object. Each object is loaded into a map as
// key = object name
// value = Wanderer class
// std::map<std::string, Wanderer> solar_system(fs::path filename);
std::map<std::string, Wanderer> create_system(fs::path filename);


//
void write_csv(fs::path output_path, std::map<std::string, Wanderer> planet_data, int rows);


std::tuple<double, double, double, double> import_settings(fs::path settings_file);
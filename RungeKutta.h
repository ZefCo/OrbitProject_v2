#pragma once
#include <tuple>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <cmath>
#include <iostream>

const int X = 0, Y = 1, RN = 0, VN = 1;

// rn and gm have the same size: Nx2.
std::vector<std::vector<double>> acc(std::vector<std::string> names, std::vector<double> masses, std::vector<std::vector<double>> rn, std::vector<std::vector<double>> dv, double h);


std::vector<std::vector<double>> vel(std::vector<std::vector<double>> vn, std::vector<std::vector<double>> dr, double h);


std::vector<std::vector<double>> nnext(std::vector<std::vector<double>> a, std::vector<std::vector<double>> d1, std::vector<std::vector<double>> d2, std::vector<std::vector<double>> d3, std::vector<std::vector<double>> d4, double h);


std::tuple< std::vector<std::vector<double>>, std::vector<std::vector<double>> > RK4(std::vector<std::string> names, std::vector<double> masses, std::vector<std::vector<double>> rn, std::vector<std::vector<double>> vn, double h);
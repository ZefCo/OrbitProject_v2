#pragma once
#include <map>
#include <fstream>
#include <iostream>
#include <filesystem>
#include "Wanderer.h"

namespace fs = std::filesystem;

// const double G = 6.67408E-11;  // changing G from [m**3 / (kg s**2)] to [AU**3 / (e d**2)] where e = mass of the earth and d = day
const double G = 880E-12;  // changing G from [m**3 / (kg s**2)] to [AU**3 / (e d**2)] where e = mass of the earth and d = day
// const double G = 1;  // changing G from [m**3 / (kg s**2)] to [AU**3 / (e d**2)] where e = mass of the earth and d = day
const int X = 0, Y = 1;


class NBodies {
    public:
        NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h);
        ~NBodies();

        void orbits();
        void time_step(int t);

        // returns the system
        std::map<std::string, Wanderer> get_system();

        // velocity
        // std::array<double, 2> acc(Wanderer Iander, std::vector<Wanderer> Jander, std::array<double, 2> kappa, double h);
        // std::array<double, 2> acc(std::array<double, 2> xy, std::vector<Wanderer> Jander, std::array<double, 2> zeta_in, double h);  
        std::array<double, 2> acc(std::array<double, 2> xy, std::vector<Wanderer> Jander, int t);
        
        // acceleration/force, = G * sum[ Mj * ((ri + hk) - rj) / ((ri + hk) - rj)**3 ]
        // takes a 2D vector in of the mass and position: mass is index 0 and position is index 1
        // assumes that the i index is not in the j index, so mr_j should be N - 1 in length
        // std::array<double, 2> vel(Wanderer Iander, std::array<double, 2> eta_in, double h);  
        std::array<double, 2> vel(std::array<double, 2> v_in, std::array<double, 2> eta_in, double h);

        // Do not use: keeping here just in case I have an issue with the other write_csv... again
        void write_csv(fs::path output_path);

        void init_table();
        std::vector<std::vector<double>> blank_table;


    private:
        int tsteps;
        double h;
        int bodies = psystem.size();

        std::map<std::string, Wanderer> psystem;  // I'm going to be accessing this a lot so it really doesn't make sense to make this private
};
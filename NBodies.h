#pragma once
#include <map>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <tuple>
#include "Wanderer.h"

namespace fs = std::filesystem;



class NBodies{
    public:
        NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h, double G);
        ~NBodies();

        // starts the time iteration using the RK4
        void orbits();

        // multiplies each M by G
        void init_system();

        // the actual RK4 step
        void time_step(int t);

        // evaluates v = v + kh
        std::tuple<double, double> vel(Wanderer Iander, int t, int k);

        // evaluates a = -G Sum(M * [(ri + ki*h) - (rj + kj*h)] / [(ri + ki*h) - (rj + kj*h)]**3)
        std::tuple<double, double> acc(Wanderer Iander, std::vector<Wanderer> Jander, int t, int k);

        std::map<std::string, Wanderer> get_system();

    private:
        const int tsteps;
        const double h;
        const int bodies;
        const double G;

        std::map<std::string, Wanderer> psystem;
};
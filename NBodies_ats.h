#pragma once
#include <map>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <tuple>
#include "Wanderer.h"

namespace fs = std::filesystem;

// const double G = 6.67408E-11;  // changing G from [m**3 / (kg s**2)] to [AU**3 / (e d**2)] where e = mass of the earth and d = day
const int X = 0, Y = 1;


class NBodies {
    public:
        NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h, double delta, double G);
        ~NBodies();

        // Oribts the bodies for time steps
        void orbits();
        
        // Does one time step. Requires the time step currently being preformed to pull the correct position and velocity data
        void adaptive_time_step(int t);
        // A single step
        std::tuple<double, double, double, double> step(Wanderer body, int t, double h);

        // returns the system
        std::map<std::string, Wanderer> get_system();

        // acceleration/force, = G * sum[ Mj * ((ri + hk) - rj) / ((ri + hk) - rj)**3 ]
        // takes a 2D vector in of the mass and position: mass is index 0 and position is index 1
        // assumes that the i index is not in the j index, so mr_j should be N - 1 in length
        std::array<double, 2> acc(std::array<double, 2> xy, std::vector<Wanderer> Jander, int t);

        // velocity
        std::array<double, 2> vel(std::array<double, 2> v_in, std::array<double, 2> eta_in, double h);

        // Do not use: keeping here just in case I have an issue with the other write_csv... again
        void write_csv(fs::path output_path);

        // Creates a blank table. This way I don't have to create a new empty table every single iteration, just create it once
        // and reuse it.
        void init_table();
        std::vector<std::vector<double>> blank_table;

        // Returns the abs between point 1 and point 2
        double euc_error(std::array<double, 2> r1, std::array<double, 2> r2, std::array<double, 2> v1, std::array<double, 2> v2);

        // A single step but with an adaptive step size
        // Finds two time steps, one at h and one at 2h, then determines rho = 30*h*delta/sqrt(e_x**2 + e_y**2). If rho > 1 we keep
        // the current step but increase h to 2h. If rho < 1 then repeat the calulation with a smaller h.
        // bool adative_step(int t);

        // the local tsize. Can and will be modified through, while tsize0 is the original step size
        double tsize;

        void adjust_tsize(double rho);



    private:
        const int tsteps;
        int bodies;
        const double delta;
        const double tsize0;
        const double G;

        std::map<std::string, Wanderer> psystem;
};
#pragma once
#include <cmath>
#include <fstream>
#include <vector>
#include <array>
#include <iostream>
// #include <math.h>

class Wanderer {
    public:
        // inits the Wanderer class. Requires a mass, and inital x and y position, an inital velocity, and a (unique) name.
        // the uniquness of the name is highly recomended because these are held in a map
        Wanderer(double mass, double x, double y, double vx, double vy, std::string name);
        // destructor
        ~Wanderer();

        // returns the mass, current r position, and name respectivly
        double get_mass();
        double get_gass();
        double get_r();
        std::array<double, 2> get_xy();
        std::array<double, 2> get_vxy();
        std::string get_name();
        // sets the size based off an input size and pushes in the inital values to the first position.
        void set_vec_size(int tsteps);

        // Inserts something into the nth position
        void nth(int n);

        // returns the nth position (index 0 and 1), and veloctiy (index 2 and 3)
        std::array<double, 4> get_nth(int n);

        // stores values for later pushing
        void storage(double xu, double yu, double vxu, double vyu);

        // stores the x, y, vx, vy positions in a vector
        std::vector<double> xn;
        std::vector<double> yn;
        std::vector<double> vxn;
        std::vector<double> vyn;

        void Gmass(double G);


    private:
        const double mass;
        double gass;
        std::string name;

        double x;
        double y;
        double vx;
        double vy;

        double  x_store;
        double  y_store;
        double vx_store;
        double vy_store;

};
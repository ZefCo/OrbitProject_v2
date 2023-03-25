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

        // returns the mass
        double get_mass();
        // returns the G*mass: figured this would be computationally easier then continually calling G
        // The other way is to factor out the G from the sum and do it once in NBodies, but I thought this
        // was more fun
        double get_gass();
        // returns sqrt(x**2 + y**2)
        double get_r();
        // Get the initial x and y values
        std::array<double, 2> get_xy();
        // Get the initial vx and vy values
        std::array<double, 2> get_vxy();
        // returns the name
        std::string get_name();
        // sets the size based off an input size and pushes in the inital values to the first position.
        void set_vec_size(int tsteps);

        // Inserts something into the nth position
        void nth(int n);

        // returns the nth position (index 0 and 1), and veloctiy (index 2 and 3)
        std::array<double, 4> get_nth(int n);

        // stores values for later pushing
        void storage(double xu, double yu, double vxu, double vyu);
        // stores the values for later comparison to determine the step size
        void storage1(double xu, double yu, double vxu, double vyu);
        // stores the values for later comparison to determine the step size
        void storage2(double xu, double yu, double vxu, double vyu);

        // stores the Xn positions in a vector
        std::vector<double> xn;
        // stores the Yn positions in a vector
        std::vector<double> yn;
        // stores the VXn positions in a vector
        std::vector<double> vxn;
        // stores the VYn positions in a vector
        std::vector<double> vyn;

        // sets the G*mass value of the Wanderer
        void Gmass(double G);

        void set_time(double h);

        double get_time();


    private:
        const double mass;
        double gass;
        std::string name;

        // OK I realize this might be Faux Pas in C++ to have these Private and then have a bunch of getters and setters for them
        // but I'm trying to limit the access to these variables, that way I don't have accidental values being shoved into them,
        // and be able to access them when needed.

        double x;
        double y;
        double vx;
        double vy;

        double  x_store;
        double  y_store;
        double vx_store;
        double vy_store;

        double  x1_store;
        double  y1_store;
        double vx1_store;
        double vy1_store;

        double  x2_store;
        double  y2_store;
        double vx2_store;
        double vy2_store;

        double tsize;

};
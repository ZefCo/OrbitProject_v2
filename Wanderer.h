#pragma once
#include <cmath>
#include <fstream>
#include <vector>
#include <array>
#include <iostream>
#include <tuple>

// not because I'm laxy, but because I have to deal with so many x and y
// this is going to be easier to index them
const int X = 0, Y = 1;


class Wanderer {
    public:
        // constructor
        Wanderer(std::string name, double mass, double x0, double y0, double vx0, double vy0);
        // destructor
        ~Wanderer();

        // get mass
        double get_mass();

        // set GMass: mass times a constant. Useful for not having to do the G*M over and over again
        // if there's a constant that has to be multiplied by the mass over and over again you can set
        // it here and save on some computations
        void set_gass(double G);

        // get GMass
        double get_gass();

        // push nth states
        void add_nth(double x, double y, double vx, double vy);

        // get nth states
        std::tuple<double, double, double, double> get_nth(int n);

        // get dv
        std::tuple<double, double> get_dv(int k);

        // set dv
        void update_dv(std::array<double, 2> dxy, int k);

        // clears dv: sets all values to 0
        void clear_dv();

        // get dr
        std::tuple<double, double> get_dr(int k);

        // set dr
        void update_dr(std::array<double, 2> dxy, int k);

        // clears dr: sets all values to 0
        void clear_dr();

        // stores the dv values for later use
        void store_dv(std::array<double, 2> dxy);
        // stores the dr values for later use
        void store_dr(std::array<double, 2> dxy);

        std::tuple<double, double> get_dr_store();
        std::tuple<double, double> get_dv_store();

    private:
        // Yes it's kind of Faux Pas to have all these private variables and then have
        // ways to easily grab them, but I want protections against weird numbers being put
        // into these variables. This way they'll change only really if I want them to
        // change

        // holds the mass
        const double mass;
        // holds the Mass * constant value
        double gass;
        // the *unique* name of the object (it should be unique, else there will be problems
        // down the road)
        const std::string name;

        // vector of the x positions
        std::vector<double> xn;
        // vector of the y positions
        std::vector<double> yn;
        // vector of the vx velocity - maybe this could be a single value as I don't really need
        // to know it's velocity at every point in time
        std::vector<double> vxn;
        // vector of the vy velocity
        std::vector<double> vyn;

        // k vector for the velocity. Note it's a 2D vector
        std::array<std::array<double, 4>, 2> dv;
        // k vector for the position. Note it's 2D: 0 = X and 1 = Y
        std::array<std::array<double, 4>, 2> dr;

        std::array<double, 2> dummy_dv;
        std::array<double, 2> dummy_dr;
};
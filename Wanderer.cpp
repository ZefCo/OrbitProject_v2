// #include <math.h>
// #include <numeric>
#include "Wanderer.h"
// #include <cmath>

Wanderer::Wanderer(double mass, 
                   double x, 
                   double y, 
                   double vx, 
                   double vy, 
                   std::string name):
mass(mass), x(x), y(y), vx(vx), vy(vy), name(name) {}


Wanderer::~Wanderer() {}

double Wanderer::get_mass() {return mass;}

double Wanderer::get_gass() {return gass;}

std::string Wanderer::get_name() {return name;}

double Wanderer::get_r() {return sqrt(pow(x, 2) + pow(y, 2));}

std::array<double, 2> Wanderer::get_xy() {return {x, y};}

std::array<double, 2> Wanderer::get_vxy() {return {vx, vy};}

std::array<double, 4> Wanderer::get_nth(int n) {return {xn[n], yn[n], vxn[n], vyn[n]};}

void Wanderer::nth(int n) {
    // this ->  x =  x_store;
    // this ->  y =  y_store;
    // this -> vx = vx_store;
    // this -> vy = vy_store;

    // std::cout << "Stored values for " << name << " are: " << x_store << "\t" << y_store << "\t" << vx_store << "\t" << vy_store << std::endl;

    xn[n] = x_store;
    yn[n] = y_store;
    vxn[n] = vx_store;
    vyn[n] = vy_store;
    // std::cout << "Size of " << name << " Xn = " << xn.size() << std::endl; 
}

void Wanderer::storage(double x, double y, double vx, double vy){
    // this ->  x =  x_store;
    // this ->  y =  y_store;
    // this -> vx = vx_store;
    // this -> vy = vy_store;
    x_store = x;
    y_store = y;
    vx_store = vx;
    vy_store = vy;

    // std::cout << "Stored X, Y, Vx, Vy = " << x_store << " " << y_store << " " << vx_store << " " << vy_store << std::endl;
}

void Wanderer::Gmass(double G) {gass = G * mass;}

void Wanderer::set_vec_size(int tsteps) {
    xn.resize(tsteps + 1), yn.resize(tsteps + 1), vxn.resize(tsteps + 1), vyn.resize(tsteps + 1);
    // std::cout << "x = " << x << " y = " << y << " vx = " << vx << " vy = " << vy << std::endl;
    xn[0] = x, yn[0] = y, vxn[0] = vx, vyn[0] = vy;

}

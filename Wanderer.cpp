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
mass(mass), x(x), y(y), vx(vx), vy(vy), name(name){}


Wanderer::~Wanderer() {}


double Wanderer::get_mass() {return mass;}


double Wanderer::get_gass() {return gass;}


std::string Wanderer::get_name() {return name;}


double Wanderer::get_r() {return sqrt(pow(x, 2) + pow(y, 2));}


std::array<double, 2> Wanderer::get_xy() {return {x, y};}


std::array<double, 2> Wanderer::get_vxy() {return {vx, vy};}


std::array<double, 4> Wanderer::get_nth(int n) {return {xn[n], yn[n], vxn[n], vyn[n]};}


void Wanderer::nth(int n) {
    xn[n] = x_store;
    yn[n] = y_store;
    vxn[n] = vx_store;
    vyn[n] = vy_store;
}


void Wanderer::storage(double x, double y, double vx, double vy){
    x_store = x;
    y_store = y;
    vx_store = vx;
    vy_store = vy;
}


void Wanderer::storage1(double x, double y, double vx, double vy){
    x1_store = x;
    y1_store = y;
    vx1_store = vx;
    vy1_store = vy;
}


void Wanderer::storage2(double x, double y, double vx, double vy){
    x2_store = x;
    y2_store = y;
    vx2_store = vx;
    vy2_store = vy;
}


void Wanderer::Gmass(double G) {gass = G * mass;}


void Wanderer::set_vec_size(int tsteps) {
    xn.resize(tsteps + 1, 0), yn.resize(tsteps + 1, 0), vxn.resize(tsteps + 1, 0), vyn.resize(tsteps + 1, 0);
    // std::cout << "x = " << x << " y = " << y << " vx = " << vx << " vy = " << vy << std::endl;
    xn[0] = x, yn[0] = y, vxn[0] = vx, vyn[0] = vy;

}


void Wanderer::set_time(double h) {tsize = h;}


double Wanderer::get_time() {return tsize;}
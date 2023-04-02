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
mass(mass), x(x), y(y), vx(vx), vy(vy), name(name){update_position(x, y, vx, vy);}


Wanderer::~Wanderer() {}


double Wanderer::get_mass() {return mass;}


double Wanderer::get_gass() {return gass;}


std::string Wanderer::get_name() {return name;}


double Wanderer::get_r() {return sqrt(pow(x, 2) + pow(y, 2));}


std::array<double, 2> Wanderer::get_xy() {return {x, y};}


std::array<double, 2> Wanderer::get_vxy() {return {vx, vy};}

std::array<std::array<double, 4>, 2> Wanderer::get_dr() {return dr;}
std::array<std::array<double, 4>, 2> Wanderer::get_dv() {return dv;}

std::tuple<double, double, double, double> Wanderer::get_nth(int n) {
    double x, y, vx, vy;
    x = xn[n]; y = yn[n]; vx = vxn[n]; vy = vyn[n];

    return {x, y, vx, vy};
}

void Wanderer::update_dr(int k) {dr[0][k] = dr_store[0]; dr[1][k] = dr_store[1];}
void Wanderer::update_dv(int k) {dv[0][k] = dv_store[0]; dr[1][k] = dv_store[1];}

void Wanderer::store_dr(std::array<double, 2> store) {dr_store = store;}
void Wanderer::store_dv(std::array<double, 2> store) {dv_store = store;}


void Wanderer::nth(int n) {
    xn[n] = x_store;
    yn[n] = y_store;
    vxn[n] = vx_store;
    vyn[n] = vy_store;
}


void Wanderer::update_position(double x, double y, double vx, double vy) {
    xn.push_back(x);
    yn.push_back(y);
    vxn.push_back(vx);
    vyn.push_back(vy);
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
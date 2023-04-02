#include "Wanderer.h"


Wanderer::Wanderer(std::string name,
                   double mass,
                   double x0,
                   double y0,
                   double vx0,
                   double vy0):
    mass(mass), name(name) {add_nth(x0, y0, vx0, vy0);}



Wanderer::~Wanderer() {}



void Wanderer::add_nth(double x, double y, double vx, double vy) 
{xn.push_back(x); yn.push_back(y); vxn.push_back(vx); vyn.push_back(vy);}



double Wanderer::get_mass() 
{return mass;}



void Wanderer::set_gass(double G) 
{gass = G*mass;
std::cout << name << " GM = " << gass << std::endl;}



double Wanderer::get_gass() 
{return gass;}



std::tuple<double, double, double, double> Wanderer::get_nth(int n)
{return {xn[n], yn[n], vxn[n], vyn[n]};}



void Wanderer::clear_dr()
{   for (int k = 0; k < 4; k++)
    {dr[X][k] = 0; dr[Y][k] = 0;}}



void Wanderer::clear_dv()
{   for (int k = 0; k < 4; k++)
    {dv[X][k] = 0; dv[Y][k] = 0;}}



void Wanderer::update_dr(std::array<double, 2> dxy, int k)
{dr[X][k] = dxy[X]; dr[Y][k] = dxy[Y];}



void Wanderer::update_dv(std::array<double, 2> dxy, int k)
{dv[X][k] = dxy[X]; dv[Y][k] = dxy[Y];}



std::tuple<double, double> Wanderer::get_dr(int k)
{return {dr[X][k], dr[Y][k]};}



std::tuple<double, double> Wanderer::get_dv(int k)
{return {dv[X][k], dv[Y][k]};}



void Wanderer::store_dv(std::array<double, 2> dxy) 
{dummy_dv = dxy;}



void Wanderer::store_dr(std::array<double, 2> dxy) 
{dummy_dr = dxy;}



std::tuple<double, double> Wanderer::get_dr_store()
{return {dummy_dr[X], dummy_dr[Y]};}



std::tuple<double, double> Wanderer::get_dv_store()
{return {dummy_dv[X], dummy_dv[Y]};}

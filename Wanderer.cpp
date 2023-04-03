#include "Wanderer.h"

// Trying using this -> because it looks like the value is getting overridden a lot and then gets
// turned into garbage
// https://stackoverflow.com/questions/24643388/variables-printing-as-nan-or-inf-instead-of-actual-value


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



void Wanderer::update_dr(int k)
{dr[X][k] = dummy_dr[X]; dr[Y][k] = dummy_dr[Y];}



void Wanderer::update_dv(int k)
{dv[X][k] = dummy_dv[X]; dv[Y][k] = dummy_dv[Y];}



std::tuple<double, double> Wanderer::get_dr(int k)
{double drx = dr[X][k]; double dry = dr[Y][k];
 if (drx != drx) {std::cout << "\t~!~!~! error in getting dr[X] !~!~!~" << std::endl;}
 if (dry != dry) {std::cout << "\t~!~!~! error in getting dr[Y] !~!~!~" << std::endl;}
 return {drx, dry};}



std::tuple<double, double> Wanderer::get_dv(int k)
{double dvx = dv[X][k]; double dvy = dv[Y][k];
 if (dvx != dvx) {std::cout << "\t~!~!~! error in getting dr[X] !~!~!~" << std::endl;}
 if (dvy != dvy) {std::cout << "\t~!~!~! error in getting dr[Y] !~!~!~" << std::endl;}
 return {dvx, dvy};}
 


void Wanderer::store_dv(std::array<double, 2> dxy) 
{this -> dummy_dv[X] = dxy[X]; this -> dummy_dv[Y] = dxy[Y];
if (dummy_dv[X] != dummy_dv[X]) {std::cout << "\t\t~~!!~~!! error in setting dummy dv X !!~~!!~~" << std::endl;}
if (dummy_dv[Y] != dummy_dv[Y]) {std::cout << "\t\t~~!!~~!! error in setting dummy dv Y !!~~!!~~" << std::endl;}}



void Wanderer::store_dr(std::array<double, 2> dxy) 
{this -> dummy_dr[X] = dxy[X]; this -> dummy_dr[Y] = dxy[Y];
if (dummy_dr[X] != dummy_dr[X]) {std::cout << "\t\t~~!!~~!! error in setting dummy dr X !!~~!!~~" << std::endl;}
if (dummy_dr[Y] != dummy_dr[Y]) {std::cout << "\t\t~~!!~~!! error in setting dummy dr Y !!~~!!~~" << std::endl;}}



std::tuple<double, double> Wanderer::get_dr_store()
{return {dummy_dr[X], dummy_dr[Y]};}



std::tuple<double, double> Wanderer::get_dv_store()
{return {dummy_dv[X], dummy_dv[Y]};}


std::string Wanderer::get_name()
{return name;}



void Wanderer::set_dr(std::array<double, 2> dxy, int k)
{dr[X][k] = dxy[X]; dr[Y][k] = dxy[Y];}



void Wanderer::set_dv(std::array<double, 2> dxy, int k)
{dv[X][k] = dxy[X]; dv[Y][k] = dxy[Y];}
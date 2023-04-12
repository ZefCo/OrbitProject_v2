#include <iostream>
#include <chrono>
#include <tuple>
#include <math.h>
#include "Wanderer.h"
#include "OrbitProject.h"
// #include "NBodies_ats.h"
#include "NBodies.h"

// g++ main.cpp jsoncpp.cpp OrbitProject.cpp Wanderer.cpp NBodies.cpp -o OP.exe
// g++ -Wall -Wextra -g main.cpp jsoncpp.cpp OrbitProject.cpp Wanderer.cpp NBodies.cpp -o OP.exe

int main() {
    int time_steps;
    double h;
    double G;
    double delta;
    int years;
    double scale;

    // time_steps = 5;
    // std::cout << "Input a value for h (double): ";
    // std::cout << "Step size set to 60*60*24 seconds (seconds in a day)" << std::endl << std::endl;

    std::tie(G, delta, h, scale) = import_settings(fs::current_path() / "Settings.json");

    std::cout << "Input the number of years to run the simulation (int): ";
    // time_steps = 50;
    // std::cout << "Value set to " << time_steps << std::endl;
    std::cin >> years;

    time_steps = years * 365 * (int)scale;
    h = h / scale;

    std::cout << "G = " << G << std::endl;
    std::cout << "h = " << h << std::endl;

    auto ts = std::chrono::high_resolution_clock::now();
    std::map<std::string, Wanderer> system;

    fs::path ao_filename = fs::current_path() / "AstronomicalObjects.json";
    fs::path ao_outfile = fs::current_path() / "Orbit_Table.csv";

    system = create_system(ao_filename);

    // NBodies motion(system, time_steps, h, 10E-6, G); // adaptive time step
    NBodies motion(system, time_steps, h, G); // non adaptive time step

    motion.init_system();

    motion.orbits();

    system = motion.get_system();

    auto te = std::chrono::high_resolution_clock::now();
    auto td = std::chrono::duration_cast<std::chrono::milliseconds>(te - ts).count();
    std::cout << "\nFinished Code in: " << td << " miliseconds\n*not including wrting to csv" << std::endl;

    // motion.write_csv(ao_outfile);

    write_csv(ao_outfile, system, time_steps);

    // double sanity = pow(2, 2);

}
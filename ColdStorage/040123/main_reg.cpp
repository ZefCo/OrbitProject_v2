#include <iostream>
#include <chrono>
#include <tuple>
// #include <math.h>
#include "Wanderer.h"
#include "OrbitProject.h"
// #include "NBodies_ats.h"
#include "NBodies.h"

// Windows complier commands
//          RG4
//      g++ main_reg.cpp jsoncpp.cpp OrbitProject.cpp Wanderer.cpp NBodies.cpp -o OP_reg.exe
//          RG4 with adaptive step
//      g++ main_ats.cpp jsoncpp.cpp OrbitProject.cpp Wanderer.cpp NBodies_ats.cpp -o OP_ats.exe
// https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html

int main() {
    int time_steps;
    double h;
    double G;
    double delta;

    // time_steps = 5;
    // std::cout << "Input a value for h (double): ";
    // std::cout << "Step size set to 60*60*24 seconds (seconds in a day)" << std::endl << std::endl;

    std::tie(G, delta, h) = import_settings(fs::current_path() / "Settings.json");
    // std::cin >> h;
    // h = 60 * 60 * 24;
    std::cout << "G = " << G << std::endl;
    // std::cout << "delta = " << delta << std::endl;  // delta not needed for this one
    std::cout << "h = " << h << std::endl;

    std::cout << "Input the number of time steps - days (int): ";
    // std::cout << "Value set to 2" << std::endl;
    std::cin >> time_steps;

    auto ts = std::chrono::high_resolution_clock::now();
    std::map<std::string, Wanderer> system;

    fs::path ao_filename = fs::current_path() / "AstronomicalObjects.json";
    fs::path ao_outfile = fs::current_path() / "Orbit_Table.csv";

    system = create_system(ao_filename);

    // NBodies motion(system, time_steps, h, 10E-6, G); // adaptive time step
    NBodies motion(system, time_steps, h, G); // non adaptive time step

    motion.orbits();

    system = motion.get_system();

    auto te = std::chrono::high_resolution_clock::now();

    // motion.write_csv(ao_outfile);

    write_csv(ao_outfile, system, time_steps);

    // double sanity = pow(2, 2);
    auto td = std::chrono::duration_cast<std::chrono::milliseconds>(te - ts).count();
    std::cout << "\nFinished Code in: " << td << " miliseconds\n*not including wrting to csv" << std::endl;

}
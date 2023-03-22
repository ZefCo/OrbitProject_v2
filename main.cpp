#include <iostream>
// #include <math.h>
#include "Wanderer.h"
#include "OrbitProject.h"
#include "NBodies.h"

// g++ main.cpp jsoncpp.cpp OrbitProject.cpp Wanderer.cpp NBodies.cpp -o OP.exe

int main() {
    double h;
    int time_steps;

    std::cout << "Input the number of time steps (int): ";
    // std::cout << "Value set to 2" << std::endl;
    std::cin >> time_steps;
    // time_steps = 5;
    std::cout << "Input a value for h (double): ";
    // std::cout << "Value set to 0.1" << std::endl << std::endl;
    std::cin >> h;
    // h = 0.1;


    std::map<std::string, Wanderer> system;

    fs::path ao_filename = fs::current_path() / "AstronomicalObjects.json";
    fs::path ao_outfile = fs::current_path() / "Orbit_Table.csv";

    system = create_system(ao_filename);

    NBodies motion(system, time_steps, h);

    motion.orbits();

    system = motion.get_system();

    // motion.write_csv(ao_outfile);

    write_csv(ao_outfile, system, time_steps);

    // double sanity = pow(2, 2);
    std::cout << "Finished" << std::endl;

}
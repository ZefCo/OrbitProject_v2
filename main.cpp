#include "OrbitProject.h"
#include "RungeKutta.h"

// Don't forget the json.cpp file and folder
// Windows
// g++ main.cpp RungeKutta.cpp OrbitProject.cpp jsoncpp.cpp -o OP.exe

const double G = 6.67408e-11;
const double scale = 1;
double dseconds = 60*60*24/scale;
const double yseconds = dseconds*365*scale;
// const int X = 0, Y = 1, RN = 0, VN = 1;

int main() {

    std::vector<double> masses; std::vector<std::string> names;
    std::vector<std::vector<double>> rn, vn;
    double years;
    std::cout << "Number of years to simulate: ";
    std::cin >> years;
    double delta = 1E-2;
    fs::path ao_filename = fs::current_path() / "AstronomicalObjects.json";

    std::cout << "years in seconds = " << yseconds << std::endl;
    std::cout << "years = " << years << std::endl;

    double time = yseconds * years;
    std::cout << "Max time = " << time << std::endl;

    std::tie(names, masses, rn, vn) = create_system(ao_filename);

    masses = gmass(masses, G);
    
    std::cout << std::endl;

    orbits(time, dseconds, names, masses, rn, vn, delta);



}
#include "OrbitProject.h"
#include "Wanderer.h"
// #include <vector>
// #include <map>

std::tuple<std::vector<double>, 
           std::vector<double>, 
           std::vector<double>, 
           std::vector<double>, 
           std::vector<double>, 
           std::vector<std::string>> loadJSON_old(fs::path filename) {
    
    std::tuple<std::vector<double>, 
               std::vector<double>, 
               std::vector<double>, 
               std::vector<double>, 
               std::vector<double>, 
               std::vector<std::string>> return_data;
    
    std::vector<double> x, y, m;
    std::vector<double> vx, vy;
    std::vector<std::string> names;

    Json::Value root;
    std::ifstream file(filename);
    file >> root;

    for (Json::Value::iterator ob = root["System"].begin(); ob != root["System"].end(); ++ob) {
        x.push_back((*ob)["x0"].asDouble());
        y.push_back((*ob)["y0"].asDouble());
        m.push_back((*ob)["m"].asDouble());
        vx.push_back((*ob)["vx0"].asDouble());
        vy.push_back((*ob)["vy0"].asDouble());
        names.push_back((*ob)["name"].asString());
    }

    return_data = {x, y, m, vx, vy, names};

    return return_data;

}


std::map<std::string, Wanderer> create_system(fs::path filename) {
    Json::Value root;
    std::ifstream file(filename);
    file >> root;
    std::map<std::string, Wanderer> system;

    for (Json::Value::iterator ob = root["System"].begin(); ob != root["System"].end(); ++ob) {
        Wanderer local_object((*ob)["name"].asString(),
                              (*ob)["m"].asDouble(), 
                              (*ob)["x0"].asDouble(), 
                              (*ob)["y0"].asDouble(), 
                              (*ob)["vx0"].asDouble(), 
                              (*ob)["vy0"].asDouble());

        // system[(*ob)["name"].asString()] = local_object;
        system.emplace((*ob)["name"].asString(), local_object);
    }

    file.close();

    return system;

}



void write_csv(fs::path output_path, std::map<std::string, Wanderer> planet_data, int rows) {

    double xn, yn, vxn, vyn;

    std::cout << "Writing to file:\n" << output_path << std::endl;

    // int rows = table.size(), cols = table[0].size();
    // int cols = planet_data.size();

    std::ofstream fileout(output_path);

    // std::cout << std::endl;

    std::string header = "t";
    for (const auto& [name, body]: planet_data) {header = header + "," + name + "_x," + name + "_y";}
    header = header + "\n";

    fileout << header;

    // std::cout << header;

    for (int t = 0; t < rows + 1; t++) {

        std::string row_data;
        fileout << t;

        for (auto& [name, body]: planet_data) {
            std::tie(xn, yn, vxn, vyn) = body.get_nth(t);
            fileout << "," << xn << "," << yn;
        }

        fileout << "\n";

        fileout << row_data;
    }

    fileout.close();

}


std::tuple<double, double, double, double> import_settings(fs::path settings_file) {
    Json::Value root;
    std::ifstream file(settings_file);
    file >> root;
    // std::tuple<double, double, double> settings;

    double G = root["G"].asDouble();
    double delta = root["delta"].asDouble();
    double time_step = root["time_step"].asDouble();
    double scale = root["scale"].asDouble();

    file.close();

    return {G, delta, time_step, scale};



}
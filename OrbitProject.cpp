#include "OrbitProject.h"


std::tuple<std::vector<std::string>, std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> create_system(fs::path filename) {
    Json::Value root;
    std::ifstream file(filename);
    file >> root;
    std::vector<std::string> names;
    std::vector<double> masses;
    std::vector<std::vector<double>> rn;
    std::vector<std::vector<double>> vn;

    for (Json::Value::iterator ob = root["System"].begin(); ob != root["System"].end(); ++ob) {
        // std::cout << (*ob)["name"] << " loaded" << std::endl;
        names.push_back((*ob)["name"].asString()); masses.push_back((*ob)["m"].asDouble());
        rn.push_back({(*ob)["x0"].asDouble(), (*ob)["y0"].asDouble()});
        vn.push_back({(*ob)["vx0"].asDouble(), (*ob)["vy0"].asDouble()});
    }

    file.close();

    return {names, masses, rn, vn};

}


std::vector<double> gmass(std::vector<double> masses, double G) {
    for (auto& mass: masses) {
        mass *= (-1)*G;
    }

    return masses;

}



void orbits(double time_steps, 
            double h,
            std::vector<std::string> names,
            std::vector<double> masses, 
            std::vector<std::vector<double>> rn, 
            std::vector<std::vector<double>> vn, 
            double delta) {

    double t = 0;
    int n = 0;

    fs::path output_path = fs::current_path() / "Orbit_Table.csv";

    std::cout << "Writing data to file as its generated:\n" << output_path << std::endl;
    std::ofstream fileout(output_path);
    fileout << "N";
    for (const auto& name: names) {
        fileout << "," << name << "_x," << name << "_y";}
    fileout << "\n";

    fileout << n;
    for (int nn = 0; nn < rn.size(); nn++) {
        fileout << "," << rn[nn][X] << "," << rn[nn][Y];
    }
    fileout << "\n";

    while (t < time_steps) {
        double dt;
        std::vector<std::vector<double>> rn0, vn0, rn1, vn1, rn2, vn2;

        double rho = 0;
        // std::cout << "rho start = " << rho << std::endl;
        // std::cout << "h start = " << h << std::endl;
        int safty = 0;

        while (rho < 1) {
            std::vector<double> diff;
            double hrho;

            dt = 2*h;

            std::tie(rn0, vn0) = RK4(names, masses, rn, vn, h);
            std::tie(rn1, vn1) = RK4(names, masses, rn0, vn0, h);
            std::tie(rn2, vn2) = RK4(names, masses, rn, vn, 2*h);

            // for (int n = 0; n < rn1.size(); n++) {
            //     std::cout << "x = " << rn[n][X] << " vs " << rn0[n][X] << " vs " << rn1[n][X] << " vs " << rn2[n][X] << std::endl;
            //     std::cout << "y = " << rn[n][Y] << " vs " << rn0[n][Y] << " vs " << rn1[n][Y] << " vs " << rn2[n][Y] << std::endl;
            // }

            diff = euc_dist(rn1, rn2);
            
            rho = wrost_rho(diff, delta, h);

            hrho = h * pow(rho, 1./4.);
            // std::cout << "\th * rho = " << hrho << " h = " << h  << " rho = " << rho << " rho**1/4 " << pow(rho, 1./4.) << std::endl;
            if (hrho < (2*h)) {h = hrho;}
            else {h = 2*h;}

            if (safty > 200) {std::cout << "Hit safty limit; t = " << t << std::endl; break;}

            safty += 1;            
        }
        // std::cout << " rho end = " << rho << std::endl;
        // std::cout << " h end = " << h << std::endl;
    
        rn = rn1; vn = vn1;

        // then write these to a file... that way I don't keep holding this data forever
        n += 1;
        fileout << n;
        for (int nn = 0; nn < rn.size(); nn++) {
            fileout << "," << rn[nn][X] << "," << rn[nn][Y];
        }
        fileout << "\n";

        // std::cout << "t = " << t << " dt = " << dt << std::endl;
        t += dt;

    }
    std::cout << "Iterations = " << n << std::endl;
    fileout.close();
}



std::vector<double> euc_dist(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b) {
    std::vector<double> c;

    c.resize(a.size());

    for (int i = 0; i < a.size(); i++) {
        double abx = pow(a[i][X] - b[i][X], 2);
        double aby = pow(a[i][Y] - b[i][Y], 2);

        double p = sqrt(abx + aby);

        if (abs(p) > 0) {c[i] = p;}
        else {c[i] = 16;}
    }

    return c;
}


double wrost_rho(std::vector<double> a, double delta, double h) {
    double rho;

    rho = h * delta / a[0];
   
    for (int i = 1; i < a.size(); i ++) {
        double iho = h * delta / a[i];
        if (rho > iho) {rho = iho;}
    }

    return rho;
}



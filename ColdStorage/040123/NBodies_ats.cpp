#include "NBodies_ats.h"


NBodies::NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h, double delta, double G):
    psystem(psystem), tsteps(tsteps), tsize0(h), delta(delta), G(G) {tsize = tsize0; bodies = psystem.size();}


NBodies::~NBodies() {}


void NBodies::orbits() {
    safties = 0;

    init_table();

    for (auto& [body_name, body]: psystem) {
        body.Gmass(G);
        body.set_vec_size(tsteps);
        body.set_time(tsize0);
    }

    for (int t = 0; t < tsteps; t++) {
        // bool insufficent_accuracy = true;

        // while (insufficent_accuracy){insufficent_accuracy = adative_step(t);}

        adaptive_time_step(t);

        for (auto& [body_name, body]: psystem) {
            body.nth(t + 1);  // note to self: the reason for storing the body xnn etc and then pushing later: I want the system to move all at once
            // and not be influenced by where things are seperetly.
            // std::cout << body_name << std::endl;
            // std::cout << "\t" << "x = " << body.xn[t + 1] << std::endl;
            // std::cout << "\t" << "y = " << body.yn[t + 1] << std::endl;
            // std::cout << "\t" << "vx = " << body.vxn[t + 1] << std::endl;
            // std::cout << "\t" << "vy = " << body.vyn[t + 1] << std::endl;
            // std::cout << "\t" << "h = " << body.get_time() << std::endl;

        }

        if ((t % 10000) == 0) {std::cout << "Finished step " << t << std::endl;}

    }

    std::cout << "Total number of safty incidents = " << safties << std::endl;
}


void NBodies::adaptive_time_step(int t) {

    for (auto& [body_name, body]: psystem) {
        double x1, y1, vx1, vy1;
        double x2, y2, vx2, vy2;
        double rho;

        bool insufficent_accuracy = true;

        int saftey = 0;

        while (insufficent_accuracy) {

            // std::cout << "Initial Time = " << body.get_time() << std::endl;

            std::tie(x1, y1, vx1, vy1) = step(body, t, body.get_time());
            std::tie(x2, y2, vx2, vy2) = step(body, t, 2*body.get_time());

            double err = euc_error_4D({x1, y1}, {x2, y2}, {vx1, vy1}, {vx2, vy2});
            // double err = euc_error_2D({x1, y1}, {x2, y2});
            // std::cout << "Err = " << err << std::endl;

            if (err == 0) {rho = 2*2*2*2;}
            else {rho = body.get_time() * delta / err;}
            std::cout << body.get_name() << ": " << body.get_time() << ", " << rho << std::endl;

            if (rho < 1.0) {

                // std::cout << "Rho = " << rho << std::endl; 
                // std::cout << "current time = " << body.get_time() << std::endl;
                double new_time = body.get_time() * pow(rho, 0.25);
                
                // std::cout << "New time = " << new_time << std::endl;
                body.set_time(new_time);
                
                saftey += 1;

                if (saftey > 100) {
                    body.storage(x1, y1, vx1, vy1);
                    // std::cout << "Hit max saftey" << std::endl;

                    // std::cout << "Rho = " << rho << std::endl;

                    // if (rho > 1) {std::cout << "What the hell? It's greater then 1" << std::endl;}
                    // else if (rho == 1) {std::cout << "What the actual fuck?" << std::endl;}
                    // else if (rho < 1) {std::cout << "What the shit? It's less then 1" << std::endl;}
                    // else {std::cout << "OK c++ is just mocking me" << std::endl;}

                    // std::exit(1);
                    safties += 1;

                    break;
                } // just in case I missed something and the while loop will not exit
            } // too much error: make h smaller
            else {

                body.set_time(2*body.get_time());
                body.storage(x1, y1, vx1, vy1);
                
                // insufficent_accuracy = false;
                // std::cout << "Rho = " << rho << std::endl;
                // std::cout << "breaking out" << std::endl;
                insufficent_accuracy = false;
                break;
            } // keep x1 and make h = 2*h
        }
        // body.storage(xnn, ynn, vxnn, vynn);
        // return {xnn, ynn, vxnn, vynn};
        // if (storage == 1) {body.storage1(xnn, ynn, vxnn, vynn);}
        // else {body.storage2(xnn, ynn, vxnn, vynn);}
    }
}



std::map<std::string, Wanderer> NBodies::get_system() {return NBodies::psystem;}



std::array<double, 2> NBodies::acc(std::array<double, 2> xy, std::vector<Wanderer> Jander, int t) {
    std::array<double, 2> ai = {0, 0};

    int length = Jander.size();

    for (int j = 0; j < length; j++) {
        std::array<double, 4> jxy = Jander[j].get_nth(t);
        double r = sqrt(((jxy[X] - xy[X])*(jxy[X] - xy[X])) + ((jxy[Y] - xy[Y])*(jxy[Y] - xy[Y])));

        double a = Jander[j].get_gass() / (r * r * r);  // note: G is moved to G*Mj, meaning it's tied to the mass

        ai[X] += a * (jxy[X] - xy[X]);
        ai[Y] += a * (jxy[Y] - xy[Y]);

    }

    return ai;

}



std::array<double, 2> NBodies::vel(std::array<double, 2> v_in, std::array<double, 2> zeta_in, double h) {
    std::array<double, 2> v_out;

    v_out[X] = v_in[X] + (h * zeta_in[X]);
    v_out[Y] = v_in[Y] + (h * zeta_in[Y]);

    return v_out;

}


void NBodies::write_csv(fs::path output_path) {

    int cols = psystem.size();

    // std::ofstream fileout(output_path);

    // std::cout << std::endl;

    std::string header = "t";
    for (const auto& [name, body]: psystem) {header = header + "," + name + "_x," + name + "_y";}
    header = header + "\n";

    // fileout << header;

    // std::cout << header;

    for (int t = 0; t < tsteps + 1; t++) {

        // std::string row_data;
        // fileout << std::to_string(r);
        // row_data = t;
        // std::cout << t;
        
        for (auto& [name, body]: psystem) {
            // row_data = row_data + "," + std::to_string(body.xn[t]) + "," + std::to_string(body.yn[t]);
            // std::cout << "\t" << body.xn[t] << "\t" << body.yn[t];
        }
        
        // row_data = row_data + "\n";
        // std::cout << "\n";

        // std::cout << row_data;

        // fileout << row_data;
    }

    // fileout.close();

}


void NBodies::init_table() {
    blank_table.resize(2);
    blank_table[0].resize(4), blank_table[1].resize(4);
}


double NBodies::euc_error_4D(std::array<double, 2> r1, 
                             std::array<double, 2> r2, 
                             std::array<double, 2> v1, 
                             std::array<double, 2> v2) {
    
    // auto f = [](float a, float b) -> float {return sqrt(a*a + b*b);};

    double R1, R2;
    double V1, V2;
    double err;

    err = (1/30) * sqrt(((r1[X] - r2[X])*(r1[X] - r2[X])) + ((r1[Y] - r2[Y])*(r1[Y] - r2[Y])) + ((v1[X] - v2[X])*(v1[X] - v1[X])) + ((v1[Y] - v2[Y])*(v1[Y] - v1[Y])));

    return err;
}


double NBodies::euc_error_2D(std::array<double, 2> r1, 
                             std::array<double, 2> r2) {
    

    double R1, R2;
    double err;

    err = sqrt(((r1[X] - r2[X])*(r1[X] - r2[X])) + ((r1[Y] - r2[Y])*(r1[Y] - r2[Y])));

    return err;
}



std::tuple<double, double, double, double> NBodies::step(Wanderer body, int t, double h) {
        std::vector<std::vector<double>> eta_table = blank_table, zeta_table = blank_table;
        std::vector<Wanderer> Jander;
        std::array<double, 2> local_eta, local_zeta;
        std::array<double, 4> xyv;
        double xn, yn, vxn, vyn;
        double xnn, ynn, vxnn, vynn;
        std::string body_name;

        xyv = body.get_nth(t);
        body_name = body.get_name();
        xn = xyv[X], yn = xyv[Y], vxn = xyv[2], vyn = xyv[3];

        for (const auto& [j_name, j_body]: psystem) {
            if (body_name != j_name) {
                Jander.push_back(j_body);
                // rm_j[0].push_back(j_body.)
            }
            // else {continue;}
        }
        // k1
        local_eta = acc({xn, yn}, Jander, t);
        eta_table[X][0] =  local_eta[X];
        eta_table[Y][0] =  local_eta[Y];

        local_zeta = {vxn, vyn};
        zeta_table[X][0] = local_zeta[X];
        zeta_table[Y][0] = local_zeta[Y];

        // k2
        local_eta = acc({xn + (zeta_table[X][0] * h/2), yn + (zeta_table[Y][0] * h/2)}, Jander, t);
        eta_table[X][1] =  local_eta[X];
        eta_table[Y][1] =  local_eta[Y];

        local_zeta = {vxn + (eta_table[X][0] * h/2), vyn + (eta_table[Y][0] * h/2)};
        zeta_table[X][1] = local_zeta[X];
        zeta_table[Y][1] = local_zeta[Y];

        // k3
        local_eta = acc({xn + (zeta_table[X][1] * h/2), yn + (zeta_table[Y][1] * h/2)}, Jander, t);
        eta_table[X][2] =  local_eta[X];
        eta_table[Y][2] =  local_eta[Y];

        local_zeta = {vxn + (eta_table[X][1] * h/2), vyn + (eta_table[Y][1] * h/2)};
        zeta_table[X][2] = local_zeta[X];
        zeta_table[Y][2] = local_zeta[Y];

        // k4
        local_eta = acc({xn + (zeta_table[X][2] * h), yn + (zeta_table[Y][2] * h)}, Jander, t);
        eta_table[X][3] =  local_eta[X];
        eta_table[Y][3] =  local_eta[Y];

        local_zeta = {vxn + (eta_table[X][2] * h), vyn + (eta_table[Y][2] * h)};
        zeta_table[X][3] = local_zeta[X];
        zeta_table[Y][3] = local_zeta[Y];

        vxnn = vxn + (h/6)*(  eta_table[X][0] + (2 *  eta_table[X][1]) + (2 *  eta_table[X][2]) +  eta_table[X][3]);
        vynn = vyn + (h/6)*(  eta_table[Y][0] + (2 *  eta_table[Y][1]) + (2 *  eta_table[Y][2]) +  eta_table[Y][3]);
        xnn =  xn + (h/6)*( zeta_table[X][0] + (2 * zeta_table[X][1]) + (2 * zeta_table[X][2]) + zeta_table[X][3]);
        ynn =  yn + (h/6)*( zeta_table[Y][0] + (2 * zeta_table[Y][1]) + (2 * zeta_table[Y][2]) + zeta_table[Y][3]);

        return {xnn, ynn, vxnn, vynn};
    
}



void NBodies::adjust_tsize(double rho) {
    tsize = tsize * pow(rho, .25);
}
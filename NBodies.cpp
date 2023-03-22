#include "NBodies.h"


NBodies::NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h):
    psystem(psystem), tsteps(tsteps), h(h) {}


NBodies::~NBodies() {}


void NBodies::orbits() {

    for (auto& [body_name, body]: psystem) {
        body.Gmass(G);
        // std::array<double, 2> xy;
        // xy = body.get_xy();
        // std::cout << body_name << " " << xy[X] << " " << xy[Y] << std::endl;
    }

    for (int t = 0; t < NBodies::tsteps; t++) {
        int bodies = NBodies::psystem.size();
        NBodies::time_step(t);

        for (auto& [body_name, body]: psystem) {
            body.nth();
            // std::cout << body_name << std::endl;
            // for (int n = 0; n < body.xn.size(); n++) {
            //     std::cout << "\t" << body.xn[n] << ", " << body.yn[n] << "  " << body.vxn[n] << ", " << body.vyn[n] << std::endl;
            // }
            // std::array<double, 4> xyv = body.get_nth(t + 1);
            // std::cout << body_name << " x = " << xyv[X] << " y = " << xyv[Y] << " vx = " << xyv[2] << " vy = " << xyv[3] << std::endl << std::endl;
        }

        if ((t % 1000) == 0) {std::cout << "Finishing step " << t << std::endl;}
    }

    // // for (int t = 0; t < NBodies::tsteps; t++) {
    //     for (auto& [body_name, body]: psystem) {
    //         std::cout << body_name << std::endl;
    //         for (int n = 0; n < body.xn.size(); n++) {
    //             std::cout << "\t" << body.xn[n] << ", " << body.yn[n] << "  " << body.vxn[n] << ", " << body.vyn[n] << std::endl;
    //         }
    //     }
    // // }
}


void NBodies::time_step(int t) {

    // int bodies = psystem.size();

    std::vector<std::vector<double>> blank_table;

    blank_table.resize(2);
    blank_table[0].resize(4), blank_table[1].resize(4);

    // int n = 0;
    for (auto& [body_name, body]: psystem) {
        std::vector<std::vector<double>> eta_table = blank_table, zeta_table = blank_table;
        std::vector<Wanderer> Jander;
        std::array<double, 2> local_eta, local_zeta, prev_eta, prev_zeta;
        std::array<double, 4> xyv;
        double xn, yn, vxn, vyn;
        double xnn, ynn, vxnn, vynn;

        xyv = body.get_nth(t);
        xn = xyv[X], yn = xyv[Y], vxn = xyv[2], vyn = xyv[3];
        // std::array<double, 2> xy, vxy;
        // xy = body.get_xy(), vxy = body.get_vxy();
        // xy = body.get_xy(), vxy = body.get_vxy();

        for (const auto& [j_name, j_body]: psystem) {
            if (body_name != j_name) {
                Jander.push_back(j_body);
                // rm_j[0].push_back(j_body.)
            }
            else {continue;}
        }
        // k1
        // std::cout << "K1" << std::endl;
        local_eta = NBodies::acc({xn, yn}, Jander, {0, 0}, h);
        local_zeta = NBodies::vel({vxn, vyn}, {0, 0}, h);
        prev_eta = local_eta, prev_zeta = local_zeta;

         eta_table[X][0] =  local_eta[X];
         eta_table[Y][0] =  local_eta[Y];
        zeta_table[X][0] = local_zeta[X];
        zeta_table[Y][0] = local_zeta[Y];

        // k2
        // std::cout << "K2" << std::endl;
        local_eta = NBodies::acc({xn, yn}, Jander, prev_zeta, h/2);
        local_zeta = NBodies::vel({vxn, vyn}, prev_eta, h/2);
        prev_eta = local_eta, prev_zeta = local_zeta;

         eta_table[X][1] =  local_eta[X];
         eta_table[Y][1] =  local_eta[Y];
        zeta_table[X][1] = local_zeta[X];
        zeta_table[Y][1] = local_zeta[Y];

        // k3
        // std::cout << "K3" << std::endl;
        local_eta = NBodies::acc({xn, yn}, Jander, prev_zeta, h/2);
        local_zeta = NBodies::vel({vxn, vyn}, prev_eta, h/2);
        prev_eta = local_eta, prev_zeta = local_zeta;

         eta_table[X][2] =  local_eta[X];
         eta_table[Y][2] =  local_eta[Y];
        zeta_table[X][2] = local_zeta[X];
        zeta_table[Y][2] = local_zeta[Y];

        // k4
        // std::cout << "K4" << std::endl;
        local_eta = NBodies::acc({xn, yn}, Jander, prev_zeta, h);
        local_zeta = NBodies::vel({vxn, vyn}, prev_eta, h);

         eta_table[X][3] =  local_eta[X];
         eta_table[Y][3] =  local_eta[Y];
        zeta_table[X][3] = local_zeta[X];
        zeta_table[Y][3] = local_zeta[Y];

        // std::cout << std::endl;

        vxnn = vxn + (h/6)*(  eta_table[X][0] + (2 *  eta_table[X][1]) + (2 *  eta_table[X][2]) +  eta_table[X][3]);
        vynn = vyn + (h/6)*(  eta_table[Y][0] + (2 *  eta_table[Y][1]) + (2 *  eta_table[Y][2]) +  eta_table[Y][3]);
         xnn =  xn + (h/6)*( zeta_table[X][0] + (2 * zeta_table[X][1]) + (2 * zeta_table[X][2]) + zeta_table[X][3]);
         ynn =  yn + (h/6)*( zeta_table[Y][0] + (2 * zeta_table[Y][1]) + (2 * zeta_table[Y][2]) + zeta_table[Y][3]);

        // std::cout << "eta Table" << std::endl;
        // std::cout << eta_table[X][0] << "\t" << eta_table[X][1] << "\t" << eta_table[X][2] << "\t" << eta_table[X][3] << std::endl;
        // std::cout << eta_table[Y][0] << "\t" << eta_table[Y][1] << "\t" << eta_table[Y][2] << "\t" << eta_table[Y][3] << std::endl;

        // std::cout << "zeta Table" << std::endl;
        // std::cout << zeta_table[X][0] << "\t" << zeta_table[X][1] << "\t" << zeta_table[X][2] << "\t" << zeta_table[X][3] << std::endl;
        // std::cout << zeta_table[Y][0] << "\t" << zeta_table[Y][1] << "\t" << zeta_table[Y][2] << "\t" << zeta_table[Y][3] << std::endl;

        // std::cout << "Xnn = " << xnn << " Ynn = " << ynn << " Vxnn = " << vxnn << " Vynn = " << vynn << std::endl; 

        body.storage(xnn, ynn, vxnn, vynn);

        // std::cout << std::endl;
    }
}


std::map<std::string, Wanderer> NBodies::get_system() {return NBodies::psystem;}


// std::array<double, 2> NBodies::acc(Wanderer Iander, std::vector<Wanderer> Jander, std::array<double, 2> eta, double h) {
std::array<double, 2> NBodies::acc(std::array<double, 2> xy, std::vector<Wanderer> Jander, std::array<double, 2> eta, double h) {
    std::array<double, 2> ai = {0, 0};
    std::array<double, 2> rn;
    // std::array<double, 2> xy = Iander.get_xy();

    int length = Jander.size();

    rn[X] = xy[X] + (eta[X] * h);
    rn[Y] = xy[Y] + (eta[Y] * h);
    // double ri = sqrt(pow(rn[X], 2) + pow(rn[Y], 2));

    // std::cout << "Ri = " << r << std::endl; 

    for (int j = 0; j < length; j++) {
        std::array<double, 2> jxy = Jander[j].get_xy();
        double r = pow((rn[X] - jxy[X]), 2) + pow((rn[Y] - jxy[Y]), 2);

        // std::cout << "Mj = " << Jander[j].get_gass() << " Rj = " << Jander[j].get_r() << " r = " << pow(r - Jander[j].get_r(), 3) << " {" << jxy[X] << ", " << jxy[Y] << "} ";
        
        double a = (-1) * Jander[j].get_gass() / pow(r, 3);  // note: G is moved to G*Mj, meaning it's tied to the mass

        // std::cout << " a = " << a << std::endl;
        
        ai[X] += a * (rn[X] - jxy[X]);
        ai[Y] += a * (rn[Y] - jxy[Y]);

        // std::cout << "Ax = " << ai[X] << " Ay = " << ai[Y] << std::endl;
    }

    return ai;

}


// std::array<double, 2> NBodies::vel(Wanderer Iander, std::array<double, 2> zeta_in, double h) {
std::array<double, 2> NBodies::vel(std::array<double, 2> v_in, std::array<double, 2> zeta_in, double h) {
    std::array<double, 2> v_out;

    // std::cout << "h = " << h << std::endl;
    // std::cout << "Vx_in = " << v_in[X] << "\tVy_in = " << v_in[Y] << std::endl;
    // std::cout << "nx = " << zeta_in[X] << "\tny = " << zeta_in[Y] << std::endl;

    v_out[X] = v_in[X] + (h * zeta_in[X]);
    v_out[Y] = v_in[Y] + (h * zeta_in[Y]);

    // std::cout << "Vx_out = " << v_out[X] << "\tVy_out = " << v_out[Y] << std::endl; 

    // return v_in + (h * zeta_in);

    return v_out;

}


void NBodies::write_csv(fs::path output_path) {

    // int rows = table.size(), cols = table[0].size();
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
#include "NBodies.h"


NBodies::NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h, double G):
    psystem(psystem), tsteps(tsteps), h(h), G(G) {}


NBodies::~NBodies() {}


void NBodies::orbits() {

    init_table();

    for (auto& [body_name, body]: psystem) {
        body.Gmass(G);
        // body.set_vec_size(tsteps);
    }

    for (int t = 0; t < NBodies::tsteps; t++) {
        time_step(t);

        // for (auto& [body_name, body]: psystem) {
        //     body.nth(t + 1);  // note to self: the reason for storing the body xnn etc and then pushing later: I want the system to move all at once
        //     // and not be influenced by where things are seperetly.

        //     // std::cout << body_name << std::endl;
        //     // std::cout << body.xn[t] << ", " << body.yn[t] << "\t" << body.vxn[t] << ", " << body.vyn[t] << std::endl; 
        //     // std::cout << body.xn[t + 1] << ", " << body.yn[t + 1] << "\t" << body.vxn[t + 1] << ", " << body.vyn[t + 1] << std::endl; 

        //     // for (int n = 0; n < body.xn.size(); n++) {
        //     //     std::cout << "\t" << body.xn[n] << ", " << body.yn[n] << "  " << body.vxn[n] << ", " << body.vyn[n] << std::endl;
        //     // }
        //     // std::array<double, 4> xyv = body.get_nth(t + 1);
        //     // std::cout << body_name << " x = " << xyv[X] << " y = " << xyv[Y] << " vx = " << xyv[2] << " vy = " << xyv[3] << std::endl << std::endl;
        // }

        if ((t % 10000) == 0) {std::cout << "Finished step " << t << std::endl;}
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

    // k1
    for (auto& [body_name, body]: psystem) {
        std::vector<std::vector<double>> dr_table = blank_table, dv_table = blank_table;
        std::vector<Wanderer> Jander;
        std::array<double, 2> dr_local, dv_local;
        std::array<double, 4> xyv;
        double xn, yn, vxn, vyn;
        double xnn, ynn, vxnn, vynn;
        double h;

        std::tie(xn, yn, vxn, vyn) = body.get_nth(t);
        h = body.get_time();

        for (const auto& [j_name, j_body]: psystem) {
            if (body_name != j_name) {
                Jander.push_back(j_body);
                // rm_j[0].push_back(j_body.)
            }
            // else {continue;}
        }

        dr_local = vel(body, t, 0);
        dv_local = acc(body, Jander, t, 0);

        body.store_dr(dr_local);
        body.store_dv(dv_local);

    }

    for (auto& [body_name, body]: psystem) {body.update_dr(0); body.update_dv(0);}

    // k2-4
    for (int k = 1; k < 4; k++) {
        for (auto& [body_name, body]: psystem) {
            std::vector<std::vector<double>> dr_table = blank_table, dv_table = blank_table;
            std::vector<Wanderer> Jander;
            std::array<double, 2> dr_local, dv_local;
            std::array<double, 4> xyv;
            double xn, yn, vxn, vyn;
            double xnn, ynn, vxnn, vynn;
            double h;

            std::tie(xn, yn, vxn, vyn) = body.get_nth(t);
            h = body.get_time();

            for (const auto& [j_name, j_body]: psystem) {
                if (body_name != j_name) {
                    Jander.push_back(j_body);
                    // rm_j[0].push_back(j_body.)
                }
                // else {continue;}
            }

            dr_local = vel(body, t, k - 1);
            dv_local = acc(body, Jander, t, k - 1);

            body.store_dr(dr_local);
            body.store_dv(dv_local);

        }

        for (auto& [body_name, body]: psystem) {body.update_dr(k); body.update_dv(k);}
    }

    for (auto& [body_name, body]: psystem) {
        std::array<std::array<double, 4>, 2> dr;
        std::array<std::array<double, 4>, 2> dv;
        std::array<double, 4> xvy;
        double xn, yn, vxn, vyn;
        double xnn, ynn, vxnn, vynn;
        double h;

        h = body.get_time();

        dr = body.get_dr(); dv = body.get_dv();

        std::tie(xn, yn, vxn, vyn) = body.get_nth(t);

        xnn =   xn + (h / 6.0) * (dr[X][0] + 2.0*dr[X][1] + 2.0*dr[X][2] + dr[X][3]);
        ynn =   yn + (h / 6.0) * (dr[Y][0] + 2.0*dr[Y][1] + 2.0*dr[Y][2] + dr[Y][3]);
        vxnn = vxn + (h / 6.0) * (dv[X][0] + 2.0*dv[X][1] + 2.0*dv[X][2] + dv[X][3]);
        vynn = vyn + (h / 6.0) * (dv[Y][0] + 2.0*dv[Y][1] + 2.0*dv[Y][2] + dv[Y][3]);

        body.update_position(xnn, ynn, vxnn, vynn);

    }
}


std::map<std::string, Wanderer> NBodies::get_system() {return NBodies::psystem;}


// // std::array<double, 2> NBodies::acc(Wanderer Iander, std::vector<Wanderer> Jander, std::array<double, 2> eta, double h) {
// std::array<double, 2> NBodies::acc(std::array<double, 2> xy, std::vector<Wanderer> Jander, std::array<double, 2> zeta_in, double h) {
//     std::array<double, 2> ai = {0, 0};
//     std::array<double, 2> rn;
//     // std::array<double, 2> xy = Iander.get_xy();

//     int length = Jander.size();

//     rn[X] = xy[X] + (zeta_in[X] * h);
//     rn[Y] = xy[Y] + (zeta_in[Y] * h);
//     // double ri = sqrt(pow(rn[X], 2) + pow(rn[Y], 2));

//     // std::cout << "Ri = " << r << std::endl; 

//     for (int j = 0; j < length; j++) {
//         std::array<double, 2> jxy = Jander[j].get_xy();
//         // double r = pow((jxy[X] - rn[X]), 2) + pow((jxy[Y] - rn[Y]), 2);
//         double r = sqrt(((jxy[X] - rn[X])*(jxy[X] - rn[X])) + ((jxy[Y] - rn[Y])*(jxy[Y] - rn[Y])));

//         // std::cout << "Mj = " << Jander[j].get_gass() << " Rj = " << Jander[j].get_r() << " r = " << pow(r - Jander[j].get_r(), 3) << " {" << jxy[X] << ", " << jxy[Y] << "} ";
        
//         // double a = Jander[j].get_gass() / pow(r, 3);  // note: G is moved to G*Mj, meaning it's tied to the mass
//         double a = Jander[j].get_gass() / (r * r * r);  // note: G is moved to G*Mj, meaning it's tied to the mass

//         // std::cout << " a = " << a << std::endl;
        
//         ai[X] += a * (jxy[X] - rn[X]);
//         ai[Y] += a * (jxy[Y] - rn[Y]);

//         // std::cout << "Ax = " << ai[X] << " Ay = " << ai[Y] << std::endl;
//     }

//     return ai;

// }


std::array<double, 2> NBodies::acc(Wanderer Iander, std::vector<Wanderer> Jander, int t, int ki) {
    std::array<double, 2> ai = {0, 0};
    double xn, yn, vxn, vyn;


    std::tie(xn, yn, vxn, vyn) = Iander.get_nth(t);
    std::array<std::array<double, 4>, 2> DR = Iander.get_dr();
    double drx, dry;
    if (ki > 0) {drx = DR[X][ki]; dry = DR[Y][ki];}
    else {drx = 0; dry = 0;}

    double Ix = (xn + drx*h);
    double Iy = (yn + dry*h);

    int length = Jander.size();

    for (int j = 0; j < length; j++) {
        double jx, jy, jvx, jvy;
        std::tie(jx, jy, jvx, jvy) = Jander[j].get_nth(t);
        std::array<std::array<double, 4>, 2> JDR = Jander[j].get_dr();
        double jh = Jander[j].get_time();
        double Jx, Jy;

        if (ki > 0 ) {Jx = (jx + JDR[X][ki]*jh); Jy = jy + JDR[Y][ki]*jh;}
        else {Jx = jx; Jy = jy;}

        double r = sqrt( ((Jx - Ix)*(Jx - Ix)) + ((Jy - Iy)*(Jy - Iy)) );

        double a = Jander[j].get_gass() / (r * r * r);  // note: G is moved to G*Mj, meaning it's tied to the mass

        ai[X] += a * (Jx - Ix);
        ai[Y] += a * (Jy - Ix);

    }

    return ai;

}



// // std::array<double, 2> NBodies::vel(Wanderer Iander, std::array<double, 2> zeta_in, double h) {
// std::array<double, 2> NBodies::vel(std::array<double, 2> v_in, std::array<double, 2> zeta_in, double h) {
//     std::array<double, 2> v_out;

//     // std::cout << "h = " << h << std::endl;
//     // std::cout << "Vx_in = " << v_in[X] << "\tVy_in = " << v_in[Y] << std::endl;
//     // std::cout << "nx = " << zeta_in[X] << "\tny = " << zeta_in[Y] << std::endl;

//     v_out[X] = v_in[X] + (h * zeta_in[X]);
//     v_out[Y] = v_in[Y] + (h * zeta_in[Y]);

//     // std::cout << "Vx_out = " << v_out[X] << "\tVy_out = " << v_out[Y] << std::endl; 

//     // return v_in + (h * zeta_in);

//     return v_out;

// }


// std::array<double, 2> NBodies::vel(Wanderer Iander, std::array<double, 2> zeta_in, double h) {
std::array<double, 2> NBodies::vel(Wanderer Iander, int ki, int t) {
    std::array<double, 2> v_out;
    double xn, yn, vxn, vyn;
    std::tie(xn, yn, vxn, vyn) = Iander.get_nth(t);
    double h = Iander.get_time();
    std::array<std::array<double, 4>, 2> DV = Iander.get_dv();
    double dvx, dvy;

    if (ki > 0) {dvx = DV[X][ki]; dvy = DV[Y][ki];}
    else {dvx = 0; dvy = 0;}

    // std::cout << "h = " << h << std::endl;
    // std::cout << "Vx_in = " << v_in[X] << "\tVy_in = " << v_in[Y] << std::endl;
    // std::cout << "nx = " << zeta_in[X] << "\tny = " << zeta_in[Y] << std::endl;

    v_out[X] = vxn + (h * dvx);
    v_out[Y] = vyn + (h * dvy);

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


void NBodies::init_table() {
    blank_table.resize(2);
    blank_table[0].resize(4), blank_table[1].resize(4);
}
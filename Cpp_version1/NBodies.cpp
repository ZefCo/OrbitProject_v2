#include "NBodies.h"

// something is odd with the k2 values... e-300? It might not be storing the correct values and indexing the wrong area.

NBodies::NBodies(std::map<std::string, Wanderer> psystem, int tsteps, double h, double G):
psystem(psystem), tsteps(tsteps), h(h), G(G), bodies(psystem.size()) {}

NBodies::~NBodies() {}



void NBodies::orbits() 
{
    for (int t = 0; t < tsteps; t++) 
    {time_step(t); 
    n_plusone(t); 
    if ((t % 10000) == 0) {std::cout << "Finished step " << t << std::endl;}}}



void NBodies::init_system()
{for (auto& [body_name, body]: psystem) {body.set_gass(G); body.clear_dr(); body.clear_dv();
double x0, y0, vx0, vy0;
std::tie(x0, y0, vx0, vy0) = body.get_nth(0);

std::cout << body_name << ": x0 = " << x0 << " y0 = " << y0 << " vx0 = " << vx0 << " vy0 = " << vy0 << std::endl;}}



void NBodies::time_step(int t)
{
    std::cout << "#### " << t << " ####" << std::endl;

    for (int kn = 0; kn < 4; kn++) {
        int k = kn - 1;
        std::map<std::string, std::array<double, 2>> dr_storage;
        std::map<std::string, std::array<double, 2>> dv_storage;
        if (k < 0) {k = 0;} // the k values for a given body are set to 0, so indexing with kn = 0 and k = 0 means the return values are 0.
                            // it's a little bit of work to be super lazy: I can use the same function over and over with no alteration, just
                            // need to make sure that the tables are set to 0 at the start of the loop
        
        for (auto& [body_name, body]: psystem) {
            std::vector<Wanderer> Jander;
            double dvx, dvy, drx, dry;

            for (auto& [jody_name, jody]: psystem) {
                if (body_name != jody_name) {Jander.push_back(jody);}
            }

            std::tie(drx, dry) = vel(body, t, k);  // returns dr/dt -> uses dv/dt to find the k value
            std::tie(dvx, dvy) = acc(body, Jander, t, k);  // returns dv/dt -> uses dr/dt to find the k vaule

            // std::cout << k << " " << drx << " " << dry << " " << dvx << " " << dvy << std::endl;

            // body.store_dr({drx, dry});
            // body.store_dv({dvx, dvy});

            dr_storage[body_name] = {drx, dry};
            dv_storage[body_name] = {dvx, dvy};

            if ((drx != drx) || (dry != dry) || (dvx != dvx) || (dvy != dvy)) {
                std::cout << "!!!! Problem at " << kn << " " << body_name << std::endl;

                if (drx != drx) {std::cout << "\tdrx!!!!" << std::endl;}
                if (dry != dry) {std::cout << "\tdry!!!!" << std::endl;}
                if (dvx != dvx) {std::cout << "\tdvx!!!!" << std::endl;}
                if (dvy != dvy) {std::cout << "\tdvy!!!!" << std::endl;}}
        }

        for (auto& [body_name, body]: psystem) {
            double dvx, dvy, drx, dry;

            for(auto const& [name, drr]: dr_storage) {
                if (name == body_name) {
                    body.set_dr(drr, kn);
                    body.set_dv(dv_storage[name], kn);
                    break;
                }
            }
            // body.update_dr(kn);
            // body.update_dv(kn);
        }
    }

}

void NBodies::n_plusone(int t){  
    for (auto& [body_name, body]: psystem) {
        double xnn, ynn, vxnn, vynn;

        // instead of grabbing vectors and such I'm just going to grab each value individually
        // helps keep my sanity
        double xn, yn, vxn, vyn;
        double drx1, drx2, drx3, drx4;
        double dry1, dry2, dry3, dry4;
        double dvx1, dvx2, dvx3, dvx4;
        double dvy1, dvy2, dvy3, dvy4;

        std::tie(xn, yn, vxn, vyn) = body.get_nth(t);
        std::tie(drx1, dry1) = body.get_dr(0); std::tie(drx2, dry2) = body.get_dr(1); std::tie(drx3, dry3) = body.get_dr(2); std::tie(drx4, dry4) = body.get_dr(3);
        std::tie(dvx1, dvy1) = body.get_dv(0); std::tie(dvx2, dvy2) = body.get_dv(1); std::tie(dvx3, dvy3) = body.get_dv(2); std::tie(dvx4, dvy4) = body.get_dv(3);

        // std::cout << "~~~Body " << body_name << " ~~~" << std::endl;
        // std::cout << "h = " << h << std::endl;
        // std::cout << "xn = " << xn << " yn = " << yn << " vxn = " << vxn << " vyn = " << vyn << std::endl;

        // std::cout << "drx1 = " << drx1 << " drx2 = " << drx2 << " drx3 = " << drx3 << " drx4 = " << drx4 << std::endl; 
        // std::cout << "dry1 = " << dry1 << " dry2 = " << dry2 << " dry3 = " << dry3 << " dry4 = " << dry4 << std::endl; 
        // std::cout << "dvx1 = " << dvx1 << " dvx2 = " << dvx2 << " dvx3 = " << dvx3 << " dvx4 = " << dvx4 << std::endl; 
        // std::cout << "dvy1 = " << dvy1 << " dvy2 = " << dvy2 << " dvy3 = " << dvy3 << " dvy4 = " << dvy4 << std::endl; 

        xnn = xn + (h / 6.0) * (drx1 + 2*drx2 + 2*drx3 + drx4);
        ynn = yn + (h / 6.0) * (dry1 + 2*dry2 + 2*dry3 + dry4);
        
        vxnn = vxn + (h / 6.0) * (dvx1 + 2*dvx2 + 2*dvx3 + dvx4);
        vynn = vyn + (h / 6.0) * (dvy1 + 2*dvy2 + 2*dvy3 + dvy4);

        // std::cout << "xnn = " << xnn << " ynn = " << ynn << " vxnn = " << vxnn << " vynn = " << vynn << std::endl;

        if ((xnn != xnn) || (ynn != ynn) || (vxnn != vxnn) || (vynn != vynn)) {
            std::cout << "!!!! Issues @ " << body_name << " !!!!" << std::endl;
            abort();}


        body.add_nth(xnn, ynn, vxnn, vynn);
        // std::cout << "clearning dr" << std::endl;
        body.clear_dr(); 
        // std::cout << "clearing dv" << std::endl;
        body.clear_dv();
    }
}




std::tuple<double, double> NBodies::vel(Wanderer Iander, int t, int k)
{
    double dummy1, dummy2;
    double vx, vy;
    double kx, ky;
    double kxn, kyn;
    double hl;

    if ((k == 0) || (k == 3)) {hl = h;}
    else {hl = h/2;}

    std::tie(dummy1, dummy2, vx, vy) = Iander.get_nth(t);
    std::tie(kx, ky) = Iander.get_dv(k);

    kxn = vx + kx*hl;
    kyn = vy + ky*hl;

    if ((vx != vx) || (vy != vy) || (kx != kx) || (ky != ky) || (hl != hl) || (k != k) || (t != t) || (h != h))
    {std::cout << "something went wrong in vel function for " << Iander.get_name() << std::endl;
    if (vx != vx) {std::cout << "~~~ vx ~~~" << std::endl;}
    if (vy != vy) {std::cout << "~~~ vy ~~~" << std::endl;}
    if (kx != kx) {std::cout << "~~~ kx ~~~" << std::endl;}
    if (ky != ky) {std::cout << "~~~ ky ~~~" << std::endl;}
    if (hl != hl) {std::cout << "~~~ hl ~~~" << std::endl;}
    if (k != k) {std::cout << "~~~ k ~~~" << std::endl;}
    if (t != t) {std::cout << "~~~ t ~~~" << std::endl;}
    if (h != h) {std::cout << "~~~ h ~~~" << std::endl;}
    }

    return {kxn, kyn};
}



std::tuple<double, double> NBodies::acc(Wanderer Iander, std::vector<Wanderer> Jander, int n, int k)
{
    double kxn = 0, kyn = 0;

    double dummy1, dummy2;
    double xi, yi;
    double xi0, yi0;
    double kix, kiy;
    double hl;
    
    if ((k == 0) || (k == 3)) {hl = h;}
    else {hl = h/2;}

    int length = Jander.size();

    std::tie(xi0, yi0, dummy1, dummy2) = Iander.get_nth(n);
    std::tie(kix, kiy) = Iander.get_dr(k);

    xi = xi0 + kix*hl;
    yi = yi0 + kiy*hl;

    for (int j = 0; j < length; j++) {
        double xj, yj;
        double xj0, yj0;
        double kjx, kjy;

        std::tie(xj0, yj0, dummy1, dummy2) = Jander[j].get_nth(n);
        std::tie(kjx, kjy) = Jander[j].get_dr(n);

        xj = xj0 + kjx*hl;
        yj = yj0 + kjy*hl;

        double rx = xi - xj; 
        double ry = yi - yj;

        double r = sqrt(rx*rx + ry*ry);

        double a = (-1) * Jander[j].get_gass() / (r*r*r);

        kxn += a * rx; 
        kyn += a * ry;
    }

    return {kxn, kyn};
}



std::map<std::string, Wanderer> NBodies::get_system() 
{return psystem;}
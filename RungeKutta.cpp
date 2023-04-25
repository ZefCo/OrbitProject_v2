#include "RungeKutta.h"


std::tuple< std::vector<std::vector<double>>, std::vector<std::vector<double>> > RK4(std::vector<std::string> names, 
                                                                                     std::vector<double> masses, 
                                                                                     std::vector<std::vector<double>> rn, 
                                                                                     std::vector<std::vector<double>> vn, 
                                                                                     double h) {
    std::vector<std::vector<double>> d0;
    std::vector<std::vector<double>> dr1, dr2, dr3, dr4;
    std::vector<std::vector<double>> dv1, dv2, dv3, dv4;
    std::vector<std::vector<double>> rnn, vnn;

    d0.resize(names.size());
    dr1.resize(names.size()); dr2.resize(names.size()); dr3.resize(names.size()); dr4.resize(names.size());
    dv1.resize(names.size()); dv2.resize(names.size()); dv3.resize(names.size()); dv4.resize(names.size());
    for (int n = 0; n < names.size(); n++) {
        d0[n].resize(2); 
        dr1[n].resize(2); dr2[n].resize(2); dr3[n].resize(2); dr4[n].resize(2);
        dv1[n].resize(2); dv2[n].resize(2); dv3[n].resize(2); dv4[n].resize(2);
    }

    for (int i = 0; i < names.size(); i++) {d0[i][X] = 0; d0[i][Y] = 0;}

    dr1 = vel(vn, d0, h);
    dv1 = acc(names, masses, rn, d0, h);

    dr2 = vel(vn, dv1, h/2);
    dv2 = acc(names, masses, rn, dr1, h/2);

    dr3 = vel(vn, dv2, h/2);
    dv3 = acc(names, masses, rn, dr2, h/2);

    dr4 = vel(vn, dv3, h);
    dv4 = acc(names, masses, rn, dr3, h);

    rnn = nnext(rn, dr1, dr2, dr3, dr4, h);     
    vnn = nnext(vn, dv1, dv2, dv3, dv4, h);

    return {rnn, vnn};

}



std::vector<std::vector<double>> nnext(std::vector<std::vector<double>> a, 
                                       std::vector<std::vector<double>> d1,
                                       std::vector<std::vector<double>> d2,
                                       std::vector<std::vector<double>> d3,
                                       std::vector<std::vector<double>> d4, 
                                       double h) {
    std::vector<std::vector<double>> an;
    an.resize(a.size());
    for (int n = 0; n < a.size(); n++) {an[n].resize(2);}

    for (int i = 0; i < a.size(); i++) {
        an[i][X] = a[i][X] + (h / 6.0)*(d1[i][X] + 2*d2[i][X] + 2*d3[i][X] + d4[i][X]);
        an[i][Y] = a[i][Y] + (h / 6.0)*(d1[i][Y] + 2*d2[i][Y] + 2*d3[i][Y] + d4[i][Y]);
    }

    return an;
}





std::vector<std::vector<double>> acc(std::vector<std::string> names, 
                                     std::vector<double> masses, 
                                     std::vector<std::vector<double>> rn,
                                     std::vector<std::vector<double>> dr, 
                                     double h){
    std::vector<std::vector<double>> an;
    int i = 0;

    an.resize(names.size());
    for (int n = 0; n < names.size(); n++) {an[n].resize(2);}

    for (auto& name : names) {
        // std::array<double, 2> ri;
        double ax = 0, ay = 0;
        if (name == "Sun") {an[i][X] = ax; an[i][Y] = ay; i += 1; continue;}

        double xi = rn[i][X] + (dr[i][X] * h);
        double yi = rn[i][Y] + (dr[i][Y] * h);

        int j = 0;
        for (auto& jame : names) {
            if (name == jame) {j += 1; continue;}

            double xj = rn[j][X] + (dr[j][X] * h);
            double yj = rn[j][Y] + (dr[j][Y] * h);

            double r = sqrt(pow(xi - xj, 2) + pow(yi - yj, 2) + 10000.0);  // adds a soft distance of 1000 meters. This way, if by some oddity the bodies
                                                                           // are on top of each other, there is no divide by 0 error

            // The original way of doing it, wich adds in an if statement.
            // if (r < 1000) {r = 1000;} // this makes r act as if it's at least 1000 meters away. It probably adds to the error, but this way
            //                           // I don't have to deal with commets or asteroids going into planets and then r ~ 0 and there's a divide
            //                           // by zero error that shows up.

            double a = masses[j] / (r*r*r);
            ax += a * (xi - xj);
            ay += a * (yi - yj);

            j += 1;
        }
        an[i][X] = ax; 
        an[i][Y] = ay;

        i += 1;
    }

    return an;
}


std::vector<std::vector<double>> vel(std::vector<std::vector<double>> vn, 
                                     std::vector<std::vector<double>> dv, 
                                     double h) {

    std::vector<std::vector<double>> vout;

    vout.resize(vn.size());
    for (int n = 0; n < vn.size(); n++) {vout[n].resize(2);}

    for (int i = 0; i < vn.size(); i++) {
        vout[i][X] = vn[i][X] + (dv[i][X] * h);
        vout[i][Y] = vn[i][Y] + (dv[i][Y] * h);
    }

    return vout;
}




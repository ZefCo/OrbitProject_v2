#include <iostream>
#include <math.h>

int main() {
    double wtf = 1.04675e-08;
    double wth = 86400;

    std::cout << "WTF = " << wtf << std::endl;
    std::cout << "WTH = " << wth << std::endl;
    std::cout << "WTF**.25 = " << pow(wtf, 0.25) << std::endl;
    std::cout << "WTF*WTH = " << wtf*wth << std::endl;
    std::cout << "WTH * WTF**.25 = " << wth * pow(wtf, 0.25) << std::endl;

}
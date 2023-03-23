#include <iostream>

int main() {
    long double x0 = 10000;
    long double y1 = 1;
    long double y2 = 0.1;
    long double y3 = 0.01;
    long double y4 = 0.001;
    long double y5 = 0.0001;
    long double y6 = 0.00001;
    long double y7 = 0.000001;
    long double y8 = 0.0000001;
    long double y9 = 0.00000001;
    long double y10 = 0.000000001;

    std::cout << "x0 size = " << sizeof(x0) << " bytes" << std::endl;
    std::cout << "y1 size = " << sizeof(y1) << " bytes" << std::endl;
    std::cout << "y2 size = " << sizeof(y2) << " bytes" << std::endl;
    std::cout << "y3 size = " << sizeof(y3) << " bytes" << std::endl;
    std::cout << "y4 size = " << sizeof(y4) << " bytes" << std::endl;
    std::cout << "y5 size = " << sizeof(y5) << " bytes" << std::endl;
    std::cout << "y6 size = " << sizeof(y6) << " bytes" << std::endl;
    std::cout << "y7 size = " << sizeof(y7) << " bytes" << std::endl;
    std::cout << "y8 size = " << sizeof(y8) << " bytes" << std::endl;
    std::cout << "y9 size = " << sizeof(y9) << " bytes" << std::endl;
    std::cout << "y10 size = " << sizeof(y10) << " bytes" << std::endl;

    std::cout << "x0 - y1 = " << x0 - y1 << std::endl;
    std::cout << "x0 - y2 = " << x0 - y2 << std::endl;
    std::cout << "x0 - y3 = " << x0 - y3 << std::endl;
    std::cout << "x0 - y4 = " << x0 - y4 << std::endl;
    std::cout << "x0 - y5 = " << x0 - y5 << std::endl;
    std::cout << "x0 - y6 = " << x0 - y6 << std::endl;
    std::cout << "x0 - y7 = " << x0 - y7 << std::endl;
    std::cout << "x0 - y8 = " << x0 - y8 << std::endl;
    std::cout << "x0 - y9 = " << x0 - y9 << std::endl;
    std::cout << "x0 - y10 = " << x0 - y10 << std::endl;
}
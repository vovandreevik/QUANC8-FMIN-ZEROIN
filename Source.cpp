#include <iostream>
#include <iomanip>
#include <cmath>
#include "Forsythe.h"
#include "quanc8.h"

int T = 0;
const long double PI = 3.141592653589793238462643383279502884;

double f(double x) {
	return 64.77 / (pow(T, 4) * pow(x, 5) * (exp(1.432 / (T * x)) - 1));
}

double zeroin(double x) {
    return 2 * sqrt(x) - cos(x * PI / 2);
}

double fmin(double z) {
    return exp(z) * (2 * pow(z, 2) - 4) + pow((2 * pow(z, 2)) - 1, 2) + exp(2 * z) - 3 * pow(z, 4);
}

int main() {
    const double UPPER_LIMIT_COEFFICIENT = 31.66675 * pow(10, -5), LOWER_LIMIT_COEFFICIENT  = 3.039830 * pow(10, -5), ABSERR = 0.0, RELERR = 1.0e-10;
    const double FMIN_INTERVAL_START = -2, FMIN_INTERVAL_END = -1;

    const double LOWER_LIMIT = -LOWER_LIMIT_COEFFICIENT * FMin(fmin, FMIN_INTERVAL_START, FMIN_INTERVAL_END, RELERR);
    const double UPPER_LIMIT = UPPER_LIMIT_COEFFICIENT *  Zeroin(zeroin, 0, 1, RELERR);

    std::cout << "FMIN: " << FMin(fmin, FMIN_INTERVAL_START, FMIN_INTERVAL_END, RELERR) << "\n";
    std::cout << "ZEROIN: " << Zeroin(zeroin, 0, 1, RELERR) << "\n";

    std::cout << "LOWER_LIMIT: " << LOWER_LIMIT << "\n";
    std::cout << "UPPER_LIMIT: " << UPPER_LIMIT << "\n";

    std::cout << "\nQUANC8\n";
    const int BEGIN = 1000, END = 9000, STEP = 1000;
    double result = 0, errest = 0, posn = 0; int nofun = 0, flag = 0;

    std::cout << std::setw(1) << "T"
        << std::setw(16) << "result"
        << std::setw(17) << "errest"
        << std::setw(19) << "nofun"
        << std::setw(10) << "posn"
        << std::setw(10) << "flag"
        << std::endl;
    for (T = BEGIN; T <= END; T += STEP) {
        quanc8(f, LOWER_LIMIT, UPPER_LIMIT, ABSERR, RELERR, &result, &errest, &nofun, &posn, &flag);
        std::cout << std::setw(3) << T
            << std::setw(20) << std::fixed << std::setprecision(10) << result
            << std::scientific << std::setw(20) << std::setprecision(10) << errest
            << std::setw(6) << std::fixed << nofun
            << std::setw(13) << std::setprecision(2) << posn
            << std::setw(7) << std::fixed << flag
            << std::endl;
    }

    std::cout << "\nResearch on the sustainability of solutions\n";
    double l1 = LOWER_LIMIT * 1.01;
    double l2 = UPPER_LIMIT * 1.01;
    double l3 = LOWER_LIMIT * 0.99;
    double l4 = UPPER_LIMIT * 0.99;
    double r1, r2, r3, r4, r5, r6, r7, r8;
    std::cout << std::setw(25) << "r1: LOWER_LIMIT + 1%, UPPER_LIMIT\n"
        << std::setw(20) << "r2: LOWER_LIMIT + 1%, UPPER_LIMIT + 1%\n"
        << std::setw(20) << "r3: LOWER_LIMIT + 1%, UPPER_LIMIT - 1%\n"
        << std::setw(20) << "r4: LOWER_LIMIT     , UPPER_LIMIT + 1%\n"
        << std::setw(20) << "r5: LOWER_LIMIT - 1%, UPPER_LIMIT\n"
        << std::setw(20) << "r6: LOWER_LIMIT - 1%, UPPER_LIMIT + 1%\n"
        << std::setw(20) << "r7: LOWER_LIMIT - 1%, UPPER_LIMIT - 1%\n"
        << std::setw(20) << "r8: LOWER_LIMIT     , UPPER_LIMIT - 1%\n"
        << std::endl;

    std::cout << std::setw(1) << "T"
        << std::setw(25) << "|result - r1|"
        << std::setw(20) << "|result - r2|"
        << std::setw(20) << "|result - r3|"
        << std::setw(20) << "|result - r4|"
        << std::setw(20) << "|result - r5|"
        << std::setw(20) << "|result - r6|"
        << std::setw(20) << "|result - r7|"
        << std::setw(20) << "|result - r8|"
        << std::endl;

    for (T = BEGIN; T <= END; T += STEP) {
        quanc8(f, LOWER_LIMIT, UPPER_LIMIT, ABSERR, RELERR, &result, &errest, &nofun, &posn, &flag);

        quanc8(f, l1, UPPER_LIMIT, ABSERR, RELERR, &r1, &errest, &nofun, &posn, &flag);
        quanc8(f, l1, l2, ABSERR, RELERR, &r2, &errest, &nofun, &posn, &flag);
        quanc8(f, l1, l4, ABSERR, RELERR, &r3, &errest, &nofun, &posn, &flag);
        quanc8(f, LOWER_LIMIT, l2, ABSERR, RELERR, &r4, &errest, &nofun, &posn, &flag);

        quanc8(f, l3, UPPER_LIMIT, ABSERR, RELERR, &r5, &errest, &nofun, &posn, &flag);
        quanc8(f, l3, l2, ABSERR, RELERR, &r6, &errest, &nofun, &posn, &flag);
        quanc8(f, l3, l4, ABSERR, RELERR, &r7, &errest, &nofun, &posn, &flag);
        quanc8(f, LOWER_LIMIT, l4, ABSERR, RELERR, &r8, &errest, &nofun, &posn, &flag);

        std::cout << std::setw(3) << T
            << std::setw(20) << std::setprecision(6) << abs(r1 - result)
            << std::setw(20) << std::setprecision(6) << abs(r2 - result)
            << std::setw(20) << std::setprecision(6) << abs(r3 - result)
            << std::setw(20) << std::setprecision(6) << abs(r4 - result)
            << std::setw(20) << std::setprecision(6) << abs(r5 - result)
            << std::setw(20) << std::setprecision(6) << abs(r6 - result)
            << std::setw(20) << std::setprecision(6) << abs(r7 - result)
            << std::setw(20) << std::setprecision(6) << abs(r8 - result)
            << std::endl;          
    }

    return 0;
}
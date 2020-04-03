#include "boost/math/special_functions/jacobi_elliptic.hpp"
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <vector>
#define PI 3.141592654

struct EllipticVariables
{
    double sn;
    double cn;
    double dn;
    void destroy()
    {
        sn = 6969;
        cn = 6969;
        dn = 6969;
    }
    void print()
    {
        std::cout << sn << " " << cn << " " << dn << "\n";
    }
};

void printArray(std::vector<double> &array)
{
    for (int id = 0; id < array.size(); ++id)
    {
        if (id == array.size() - 1)
        {
            std::cout << array.at(id);
        }
        else
        {
            std::cout << array.at(id) << " , ";
        }
    }
    std::cout << "\n";
}

EllipticVariables jacobiAmplitude(double q, double k)
{
    auto result = new EllipticVariables;
    //the variable q is on the second spot
    //the complex number k  (constant) is the first argument
    result->sn = boost::math::jacobi_sn(k, q);
    result->cn = boost::math::jacobi_cn(k, q);
    result->dn = boost::math::jacobi_dn(k, q);
    if (!isnan(result->sn) && !isnan(result->cn) && !isnan(result->dn))
        return *result;
    result->destroy();
    return *result;
}

double A_Fct(double spin, double oddSpin, double a1, double a2, double theta)
{
    auto j2_Component = oddSpin * sin(theta * PI / 180.0);
    auto result = a2 * (1.0 - j2_Component / spin) - a1;
    if (!isnan(result))
        return result;
    return 6969;
}

double u_Fct(double spin, double oddSpin, double a1, double a2, double a3, double theta)
{
    auto result = static_cast<double>((a3 - a1) / A_Fct(spin, oddSpin, a1, a2, theta));
    if (!isnan(result))
        return result;
    return 6969;
}

double v0_Fct(double spin, double oddSpin, double a1, double a2, double theta)
{
    auto j1_Component = oddSpin * cos(theta * PI / 180.0);
    auto result = -1.0 * static_cast<double>((a1 * j1_Component) / A_Fct(spin, oddSpin, a1, a2, theta));
    if (!isnan(result))
        return result;
    return 6969;
}

double k_Fct(double spin, double oddSpin, double a1, double a2, double a3, double theta)
{
    auto result = sqrt(u_Fct(spin, oddSpin, a1, a2, a3, theta));
    if (!isnan(result))
        return result;
    return 6969;
}

double V_Rotor(double q, double spin, double oddSpin, double a1, double a2, double a3, double theta)
{
    auto I = spin;
    auto k = k_Fct(I, oddSpin, a1, a2, a3, theta);
    if (k == 6969)
        return 6969;
    auto v0 = v0_Fct(I, oddSpin, a1, a2, theta);
    if (v0 == 6969)
        return 6969;
    auto JA = jacobiAmplitude(q, k);
    auto RotorPotential = (I * (I + 1.0) * pow(k, 2) + pow(v0, 2)) * pow(JA.sn, 2) + (2.0 * I + 1.0) * v0 * JA.cn * JA.dn;
    if (!isnan(RotorPotential))
        return RotorPotential;
    return 6969;
}

int main()
{
    std::cout << "OK"
              << "\n";
    std::vector<double> vTable;
    std::vector<double> qTable;
    //set of parameters
    const double i1 = 20.0;
    const double i2 = 100.0;
    const double i3 = 40.0;
    auto a1 = static_cast<double>(1.0 / (2.0 * i1));
    auto a2 = static_cast<double>(1.0 / (2.0 * i2));
    auto a3 = static_cast<double>(1.0 / (2.0 * i3));
    const double I = 45.0 / 2.0;
    const double j = 6.5;
    const double theta = 30.0;
    auto startTime = std::chrono::high_resolution_clock::now();
    //using u instead of q
    for (double u = -8.0; u < 8; u += 0.1)
    {
        //the first value is the parameter m (or k^2)
        //the second value represents the variable q
        qTable.emplace_back(u);
        auto vRot = V_Rotor(u, I, j, a1, a2, a3, theta);
        if (vRot != 6969)
            vTable.emplace_back(vRot);
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    if (qTable.size() == vTable.size())
    {
        printArray(qTable);
        printArray(vTable);
    }
    std::cout << "Process took..." << static_cast<double>(duration / 1000) << " s\n";
}
#include "boost/math/special_functions/jacobi_elliptic.hpp"
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <vector>

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

int main()
{
    std::cout << "OK"
              << "\n";
    std::vector<double> snValues;
    double k = 4.0;
    auto startTime = std::chrono::high_resolution_clock::now();
    for (double u = -8.0; u < 8.1; u += 0.1)
    {
        //the first value is the parameter m (or k^2)
        //the second value represents the variable q
        auto currentValue = boost::math::jacobi_sn(k, u);
        snValues.emplace_back(currentValue);
    }
    auto vars = jacobiAmplitude(2.0, 3.0);
    vars.print();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    // printArray(snValues);
    std::cout << "Process took..." << static_cast<double>(duration / 1000) << " s\n";
}
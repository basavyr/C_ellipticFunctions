# Calculating the Elliptic Functions in `C++`

Constructing the Jacobi elliptic functions using the `boost`[^1] library within a C++ algorithm.

The algorithm successfully calculates the potential energy of a triaxial rigid nucleus.
The expression of the potential $V$ is expressed in terms of the general coordinate $q$.  

* $q$ is determined with the help of the `Jacobi Amplitude` -> the elliptic variables $sn$, $cn$ and $dn$ are required.
* `boost` library provides a direct way of computing the mentioned functions.


[^1]:https://www.boost.org/
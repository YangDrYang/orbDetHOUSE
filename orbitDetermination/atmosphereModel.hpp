
#include <Eigen/Dense>
#include <vector>
#include "satRefSys.hpp"

// #include "coordTrans.hpp"

// Use C atmospheric model header
extern "C" {
    #include <nrlmsise-00.h>
}



// Calculates drag Force (in Newtons) from nrlmsise-00 atmospheric model
// Given:
// rSat (3d position vector [m])
// vSat (3d velocity vector [m/s]) 
// IERS (IERS Struct)
// Area (frontal area [m^2])
// Cd (coefficient of drag [no unit])
// t (modified Julian date)
// Returns:
// dragForce (3d force vector [N])
Eigen::Vector3d calculateDragForce(const Eigen::Vector3d& rSat, 
                            const Eigen::Vector3d& vSat,
                            IERS iers,
                            double Area,
                            double Cd,
                            double t);
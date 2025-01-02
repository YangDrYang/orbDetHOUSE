
#include "satRefSys.hpp"
#include <Eigen/Dense>
#include <vector>

// #include "coordTrans.hpp"

// Use C atmospheric model header
extern "C"
{
#include <nrlmsise-00/nrlmsise-00.h>
}

using namespace Eigen;
using namespace std;

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
Vector3d calculateDragForce(const Vector3d &rSat,
                            const Vector3d &vSat,
                            IERS iers,
                            double Area,
                            double Cd,
                            double t);
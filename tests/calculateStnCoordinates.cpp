#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXd readCSV(const string &filename, int headerLinesToSkip)
{
    ifstream file(filename);
    vector<vector<double>> data;

    string line;
    int linesSkipped = 0;

    while (getline(file, line))
    {
        if (linesSkipped < headerLinesToSkip)
        {
            linesSkipped++;
            continue; // Skip header lines
        }

        istringstream iss(line);
        string value;
        vector<double> row;

        while (getline(iss, value, ','))
        {
            try
            {
                row.push_back(stod(value));
            }
            catch (const invalid_argument &e)
            {
                // Handle conversion errors, e.g., non-numeric values
                cerr << "Invalid argument: " << e.what() << endl;
                // You can choose to skip this value or handle it differently
                row.push_back(0.0); // Default value if conversion fails
            }
        }

        data.push_back(row);
    }

    int rows = static_cast<int>(data.size());
    int cols = static_cast<int>(data[0].size());

    MatrixXd matrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            matrix(i, j) = data[i][j];
        }
    }

    return matrix;
}

#include <iostream>
#include <vector>
#include <cmath>

#include <iostream>
#include <vector>
#include <cmath>

int findClosestIndex(const VectorXd &tSec, double target)
{
    double minDiff = std::numeric_limits<double>::max();
    int closestIndex = -1;

    for (int i = 0; i < tSec.size(); i++)
    {
        double diff = std::abs(tSec[i] - target);
        if (diff < minDiff)
        {
            minDiff = diff;
            closestIndex = i;
        }
    }

    return closestIndex;
}

// Function to interpolate the given data points using Lagrange's formula
// xi corresponds to the new data point whose value is to be obtained
// n represents the number of known data points
double interpolate(const MatrixXd &f, double xi, int n)
{
    double result = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double term = f(i, 1);
        for (int j = 0; j < n; ++j)
        {
            if (j != i)
                term *= (xi - f(j, 0)) / (f(i, 0) - f(j, 0));
        }
        result += term;
    }

    return result;
}

// Function to compute Greenwich Mean Sidereal Time (GMST) in degrees
double jd2Gmst(double jd)
{
    // Find the Julian Date of the previous midnight, jd0
    double jd0 = std::floor(jd - 0.5);
    double h = (jd - jd0) * 24;  // Time in hours past previous midnight
    double d = jd - 2451545.0;   // Compute the number of days since J2000
    double d0 = jd0 - 2451545.0; // Compute the number of days since J2000
    double t = d / 36525;        // Compute the number of centuries since J2000

    // Calculate GMST in hours (0h to 24h) ... then convert to degrees
    double gmst = fmod(6.697374558 + 0.06570982441908 * d0 + 1.00273790935 * h + 0.000026 * (t * t), 24.0) * 15.0;

    return gmst;
}

// Function to compute Greenwich Apparent Sidereal Time (GAST) in degrees
double jd2Gast(double jd)
{
    // Calculate mean sidereal time
    double thetaM = jd2Gmst(jd);

    // Compute the number of centuries since J2000
    double t = (jd - 2451545.0) / 36525;

    // Compute mean obliquity of the ecliptic (epsilonM) in degrees
    double epsilonM = 23.439291 - 0.0130111 * t - 1.64e-07 * (t * t) + 5.04e-07 * (t * t * t);

    // Compute nutations in obliquity and longitude (degrees)
    double L = 280.4665 + 36000.7698 * t;
    double dL = 218.3165 + 481267.8813 * t;
    double omega = 125.04452 - 1934.136261 * t;

    // Calculate nutations
    double dPsi = -17.20 * sin(omega) - 1.32 * sin(2 * L) - 0.23 * sin(2 * dL) + 0.21 * sin(2 * omega);
    double dEpsilon = 9.20 * cos(omega) + 0.57 * cos(2 * L) + 0.10 * cos(2 * dL) - 0.09 * cos(2 * omega);

    // Convert units from arc-seconds to degrees
    dPsi /= 3600.0;
    dEpsilon /= 3600.0;

    // Calculate Greenwich Apparent Sidereal Time (GAST) in degrees
    double gast = fmod(thetaM + dPsi * cos(epsilonM + dEpsilon), 360.0);

    return gast;
}

Vector2d measurementModel(double jd, const Vector3d &rho)
{
    Vector2d z;

    // Calculate GAST
    double gast = jd2Gast(jd);
    cout << "Julian date: " << setprecision(16) << jd << endl;
    cout << "Greenwich Mean Sidereal Time (GMSA): " << jd2Gmst(jd) << " degrees" << endl;
    cout << "Greenwich Apparent Sidereal Time (GAST): " << gast << " degrees" << endl;
    cout << "Greenwich Apparent Sidereal Time (GAST): " << gast / 180 * M_PI << " radians" << endl;

    // right ascension angle
    z(0) = atan2(rho(1), rho(0)) + gast / 180 * M_PI;
    if (z(0) < 0)
    {
        z(0) += 2 * M_PI;
    };

    // declination angle
    z(1) = asin(rho(2) / rho.norm());

    // cout << "calculated measurements: " << z(0) << "\t" << z(1) << endl;

    return z;
}

MatrixXd h2xJacobian(Vector3d rho)
{
    MatrixXd hMatrix(2, 3);
    double rho01 = sqrt(rho(0) * rho(0) + rho(1) * rho(1));
    double rho2 = rho.norm() * rho.norm();
    hMatrix << -rho(1) / rho01, rho(0) / rho01, 0,
        -rho(0) * rho(2) / (rho2 * rho01), -rho(1) * rho(2) / (rho2 * rho01), rho01 / rho2;

    return hMatrix;
}

int main(int argc, char *argv[])
{
    string measFile = "ccdata/meas_data_id_46984.csv";
    cout << "Reading measurement file: " << measFile << endl;

    // Read angular measurements from the exiting file
    int headerLinesToSkip = 1; // Number of lines to skip as header
    MatrixXd matMeas = readCSV(measFile, headerLinesToSkip);

    int dimMeas = 2;
    MatrixXd angMeas(matMeas.rows(), dimMeas + 1);

    // Save the seventh column of MJD
    angMeas.col(0) = matMeas.col(6);
    // Save the last two columns of angles into radians
    angMeas.rightCols(dimMeas) = matMeas.rightCols(dimMeas) / 180 * M_PI;
    Matrix2d wMatrix;
    wMatrix << 1 / (5 / 3600 / 180 * M_PI) / (5 / 3600 / 180 * M_PI), 0,
        0, 1 / (5 / 3600 / 180 * M_PI) / (5 / 3600 / 180 * M_PI);

    // cout << "angMeas:\t" << endl
    //      << angMeas << endl;

    string orbFile = "refdata/od_ecef_id_46984.csv";
    cout << "Reading Sentinel 6A ephemeris file (ECEF): " << orbFile << endl;

    int dimState = 3;
    // Read reference satellite orbit ephemeris from the exiting file
    headerLinesToSkip = 1; // Number of lines to skip as header
    MatrixXd matOrb = readCSV(orbFile, headerLinesToSkip);
    // Save the MJD and ECEF coordinates
    MatrixXd orbECEF = matOrb.block(0, 1, matOrb.rows(), dimState + 1);
    MatrixXd orbECEFInterp(angMeas.rows(), dimState + 1);

    int nOrd = 10;
    for (int i = 0; i < angMeas.rows(); i++)
    {
        int ind = findClosestIndex(orbECEF.col(0), angMeas(i, 0));
        MatrixXd xData(nOrd + 1, 2);
        xData.col(0) = orbECEF.block(ind - 5, 0, nOrd + 1, 1);
        xData.col(1) = orbECEF.block(ind - 5, 1, nOrd + 1, 1);
        double x = interpolate(xData, angMeas(i, 0), nOrd);

        MatrixXd yData(nOrd + 1, 2);
        yData.col(0) = orbECEF.block(ind - 5, 0, nOrd + 1, 1);
        yData.col(1) = orbECEF.block(ind - 5, 2, nOrd + 1, 1);
        double y = interpolate(yData, angMeas(i, 0), nOrd);

        MatrixXd zData(nOrd + 1, 2);
        zData.col(0) = orbECEF.block(ind - 5, 0, nOrd + 1, 1);
        zData.col(1) = orbECEF.block(ind - 5, 3, nOrd + 1, 1);
        double z = interpolate(zData, angMeas(i, 0), nOrd);

        orbECEFInterp(i, 0) = angMeas(i, 0);
        orbECEFInterp(i, 1) = x;
        orbECEFInterp(i, 2) = y;
        orbECEFInterp(i, 3) = z;
    }

    // Station ECEF coordinate
    Vector3d stnECEF0;
    stnECEF0 << -2730164.38085497, 3714370.03545681, 4393632.92952871;
    // Initial information matrix for station ECEF coordinate
    Matrix3d Lambda;
    Lambda << 1 / 1000, 0, 0,
        0, 1 / 1000, 0,
        0, 0, 1 / 500;

    int numMatrices = angMeas.rows();
    int numRowsPerMatrix = dimMeas;
    // Determine the total number of rows in the concatenated matrix
    int totalRows = numMatrices * numRowsPerMatrix;

    // Create a matrix to hold the concatenated H and matrices
    MatrixXd concatenatedHMatrix(totalRows, dimState);
    VectorXd concatenateddZVector(totalRows);
    MatrixXd concatenatedWMatrix(totalRows, totalRows);

    Vector3d deltaRStn0 = Vector3d::Zero(3);
    int iter = 0;
    // while (deltaRStn0.norm() > 1e-8 && iter < 1)
    while (iter < 1)
    {
        // for (int i = 0; i < angMeas.rows(); i++)
        for (int i = 0; i < 1; i++)
        {
            Vector3d rho = orbECEFInterp.row(i).tail(3).transpose() - stnECEF0;
            cout << "rho:\t" << rho.transpose() << endl;
            cout << "orbit:\t" << orbECEFInterp.row(i) << endl;
            Vector2d zVector = measurementModel(orbECEFInterp(i, 0) + 2400000.5, rho);
            cout << "calculated z:\t" << zVector.transpose() << endl;

            MatrixXd hMatrix = h2xJacobian(rho);

            // Calculate the starting row index for concatenation
            int startRow = i * numRowsPerMatrix;

            // Copy the hMatrix to the corresponding block in the concatenatedHMatrix
            concatenatedHMatrix.block(startRow, 0, numRowsPerMatrix, hMatrix.cols()) = hMatrix;
            cout << "observed z:\t" << angMeas.row(i).tail(2) << endl;
            concatenateddZVector.block(startRow, 0, numRowsPerMatrix, 1) = angMeas.row(i).tail(2).transpose() - zVector;

            // Concatenate the W matrix into a diagonal matrix
            concatenatedWMatrix.setZero();
            concatenatedWMatrix.block(startRow, startRow, 2, 2) = wMatrix;
        }

        MatrixXd tmpMat = Lambda + concatenatedHMatrix.transpose() * concatenatedWMatrix * concatenatedHMatrix;
        // deltaRStn0 = tmpMat.inverse();
        deltaRStn0 = tmpMat.inverse() * (Lambda * deltaRStn0 + concatenatedHMatrix.transpose() * concatenatedWMatrix * concatenateddZVector);

        iter += 1;
        stnECEF0 += deltaRStn0;
        cout << "The " << iter << "th iteration." << endl;
        cout << "The updated station ECEF coordinate:" << endl
             << stnECEF0.transpose() << endl;
    }
}

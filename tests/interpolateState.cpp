#include "testOrbDet.hpp"
#include <iostream>
using namespace std;

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

void writeCSV(const string &filename, const MatrixXd &data, const vector<string> &header)
{
    ofstream file(filename);
    if (file.is_open())
    {
        // Write header
        for (const auto &col : header)
        {
            file << col << ",";
        }
        file << "\n";

        // Set precision
        file << fixed << setprecision(DBL_DIG);

        // Write data
        for (int i = 0; i < data.rows(); ++i)
        {
            for (int j = 0; j < data.cols(); ++j)
            {
                file << data(i, j) << ",";
            }
            file << "\n";
        }

        file.close();
        cout << "Data written to " << filename << endl;
    }
    else
    {
        cerr << "Failed to open file: " << filename << endl;
    }
}

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

int main(int argc, char *argv[])
{
    string orignialFile = "refdata/od_eci_id_46984.csv";
    int headerLinesToSkip = 1; // Number of lines to skip as header
    MatrixXd originalState = readCSV(orignialFile, headerLinesToSkip);

    // read angular measurements from the exiting file
    string measFile = "ccdata/meas_data_id_46984.csv";
    headerLinesToSkip = 1; // Number of lines to skip as header
    MatrixXd matMeas = readCSV(measFile, headerLinesToSkip);

    MatrixXd interpolatedState(matMeas.rows(), 7);
    interpolatedState.col(0) = matMeas.col(6).array();

    int nOrd = 10;
    for (int i = 0; i < matMeas.rows(); i++)
    {
        int ind = findClosestIndex(originalState.col(0), matMeas(i, 0));

        MatrixXd xData(nOrd + 1, 2);
        xData.col(0) = orignalState.block(ind - 5, 0, nOrd + 1, 1);
        xData.col(1) = orignalState.block(ind - 5, 1, nOrd + 1, 1);
        interpolatedState(i, 1) = interpolate(xData, orignalState(i, 0), nOrd);

        MatrixXd yData(nOrd + 1, 2);
        yData.col(0) = orignalState.block(ind - 5, 0, nOrd + 1, 1);
        yData.col(1) = orignalState.block(ind - 5, 2, nOrd + 1, 1);
        interpolatedState(i, 2) = interpolate(yData, orignalState(i, 0), nOrd);

        MatrixXd zData(nOrd + 1, 2);
        zData.col(0) = orignalState.block(ind - 5, 0, nOrd + 1, 1);
        zData.col(1) = orignalState.block(ind - 5, 3, nOrd + 1, 1);
        interpolatedState(i, 3) = interpolate(zData, orignalState(i, 0), nOrd);

        MatrixXd vxData(nOrd + 1, 2);
        vxData.col(0) = orignalState.block(ind - 5, 0, nOrd + 1, 1);
        vxData.col(1) = orignalState.block(ind - 5, 4, nOrd + 1, 1);
        interpolatedState(i, 4) = interpolate(vxData, orignalState(i, 0), nOrd);

        MatrixXd vyData(nOrd + 1, 2);
        vyData.col(0) = orignalState.block(ind - 5, 0, nOrd + 1, 1);
        vyData.col(1) = orignalState.block(ind - 5, 5, nOrd + 1, 1);
        interpolatedState(i, 5) = interpolate(vyData, orignalState(i, 0), nOrd);

        MatrixXd zData(nOrd + 1, 2);
        vzData.col(0) = orignalState.block(ind - 5, 0, nOrd + 1, 1);
        vzData.col(1) = orignalState.block(ind - 5, 6, nOrd + 1, 1);
        interpolatedState(i, 6) = interpolate(vzData, orignalState(i, 0), nOrd);
    }

    vector<string> header = {"MJD", "Interpolated_X", "Interpolated_Y", "Interpolated_Z", "Interpolated_VX", "Interpolated_VY", "Interpolated_VZ"};
    writeCSV("refdata/od_ref_id_46984.csv", satECIData, header);

    return 0;
}
#ifndef EIGEN_CSV_H
#define EIGEN_CSV_H

#include <Eigen/Dense>

#include <cfloat>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
using namespace Eigen;
using namespace std;

namespace EigenCSV
{

    template <typename T>
    void write(const DenseBase<T> &A,
               const string &filename)
    {

        using namespace std;

        ofstream file;

        file.open(filename);

        file << setprecision(DBL_DIG);

        for (int i = 0; i < A.rows(); i++)
        {
            for (int j = 0; j < A.cols(); j++)
                file << A(i, j) << ',';
            file << endl;
        }

        file.close();
    }

    template <typename T>
    void write(const DenseBase<T> &A,
               const vector<string> &header,
               const string &filename,
               bool append = false)
    {

        using namespace std;

        ofstream file;
        if (append == true)
        {
            file.open(filename, ios::app);
        }
        else
        {
            file.open(filename);
        }

        for (unsigned int i = 0; i < header.size(); i++)
            file << header[i] << ',';
        file << endl;

        file << setprecision(DBL_DIG);

        for (int i = 0; i < A.rows(); i++)
        {
            for (int j = 0; j < A.cols(); j++)
                file << A(i, j) << ',';
            file << endl;
        }

        file.close();
    }

    template <typename T>
    void read(const string &filename,
              bool header, bool resize, DenseBase<T> &A)
    {

        using namespace std;

        ifstream file;

        file.open(filename);

        string line, cell;

        if (header)
        {
            getline(file, line);
        }

        vector<vector<double>> vals;

        while (getline(file, line))
        {

            stringstream str(line);

            vector<double> valr;

            while (getline(str, cell, ','))
                valr.push_back(stod(cell));

            vals.push_back(valr);
        }

        int rows, cols;

        if (resize)
        {
            rows = vals.size();
            cols = vals[0].size();
            T &B = (T &)A;
            B.resize(rows, cols);
        }
        else
        {
            rows = A.rows();
            cols = A.cols();
        }

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                A(i, j) = vals[i][j];
    }

}

#endif

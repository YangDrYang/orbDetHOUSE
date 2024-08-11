#include "testOrbDet.hpp"
#include <iostream>
using namespace std;

erp_t erpt;
IERS iersInstance;

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

string getFileExtension(const string &filename)
{
    size_t dotPos = filename.rfind('.');
    if (dotPos != string::npos)
    {
        return filename.substr(dotPos + 1);
    }
    return ""; // Return an empty string if no file extension is found
}

MatrixXd readFile(const string &filename, int numHeaderLinesToSkip)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        return Eigen::MatrixXd();
    }

    string line;
    for (int i = 0; i < numHeaderLinesToSkip; ++i)
    {
        getline(file, line); // Skip header lines
    }

    // string line;
    // getline(file, line); // Skip header line 1
    // getline(file, line); // Skip header line 2

    vector<vector<double>> satelliteData;
    while (getline(file, line))
    {
        if (line.find("99") == 0)
        {
            continue; // Skip line starting with "99"
        }
        istringstream iss(line);
        vector<double> data;
        double value;
        while (iss >> value)
        {
            data.push_back(value);
        }
        satelliteData.push_back(data);
    }

    file.close();

    // Convert satelliteData to Eigen MatrixXd
    int numDataPoints = satelliteData.size();
    int numColumns = satelliteData[0].size();
    MatrixXd satelliteMatrix(numDataPoints, numColumns);
    for (int i = 0; i < numDataPoints; ++i)
    {
        for (int j = 0; j < numColumns; ++j)
        {
            satelliteMatrix(i, j) = satelliteData[i][j];
        }
    }

    return satelliteMatrix;
}

int main(int argc, char *argv[])
{
    // default  ECEF file
    string orbECEFFile;
    int flagGPS = 0, flagCPF = 0;
    if (argc == 2)
    {
        orbECEFFile = argv[1];
        // string orbECEFFile = "refdata/od_ecef_id_46984.csv";
        // string orbECEFFile = "refdata/sentinel6a_cpf_210515_13501.eum";
        string fileExtension = getFileExtension(orbECEFFile);
        if (fileExtension == "csv")
            flagGPS = 1;
        if (fileExtension == "eum")
            flagCPF = 1;
    };

    cout << "flag\t" << flagGPS << "\t" << flagCPF << endl;
    erpt = {.n = 0};
    string erpFile = "auxdata/cod21587.erp";
    readerp(erpFile, &erpt);

    // read angular measurements from the exiting file
    string measFile = "ccdata/meas_data_id_46984.csv";
    int headerLinesToSkip = 1; // Number of lines to skip as header
    MatrixXd matMeas = readCSV(measFile, headerLinesToSkip);

    // extract the corresponding MJD for each set measurement
    VectorXd tMJD = matMeas.col(6).array();
    // cout << setprecision(18) << tMJD << endl;

    // Jinlin, Jinlin, ECEF coordinate, unit: m, m/s
    VectorXd stnECEF(6);
    stnECEF << -2730000, 3714000, 4393000, 0, 0, 0;
    // cout << stnECEF << endl;

    // Create a matrix to store MJD and sntECI values
    MatrixXd stnECIData(tMJD.size(), 7);

    for (int i = 0; i < tMJD.size(); i++)
    {
        VectorXd stnECI = VectorXd::Zero(6);
        // cout << tMJD[i] << endl;
        getIERS(tMJD[i], erpt, iersInstance);
        ecef2eciVec_sofa(tMJD[i], iersInstance, stnECEF, stnECI);
        // cout << stnECI << endl;

        // Format the output row using std::ostringstream
        ostringstream oss;
        oss << fixed << setprecision(10) << tMJD[i];
        for (int j = 0; j < stnECI.size(); j++)
        {
            oss << " " << fixed << setprecision(DBL_DIG) << stnECI[j];
        }

        // Assign the formatted row to stnECIData
        istringstream issStn(oss.str());
        VectorXd tempRow(stnECIData.cols());
        for (int j = 0; j < tempRow.size(); j++)
        {
            issStn >> tempRow[j];
        }
        stnECIData.row(i) = tempRow;
    }

    // // Print stnECIData rows
    // for (int i = 0; i < stnECIData.rows(); i++)
    // {
    //     cout << stnECIData.row(i) << endl;
    // }
    // Prepare the header
    vector<string> header = {"MJD", "X_ECI", "Y_ECI", "Z_ECI", "VX_ECI", "VY_ECI", "VZ_ECI"};

    // Write MJD and stnECI to CSV file
    writeCSV("ccdata/stn_eci_coordinates.csv", stnECIData, header);

    // perform coordinate transformation for satellite via a CPF file
    if (flagGPS)
    {
        // read reference orbit in ECEF from the exiting file
        headerLinesToSkip = 1; // Number of lines to skip as header
        MatrixXd matODECEF = readCSV(orbECEFFile, headerLinesToSkip);
        // cout << "matODECEF" << matODECEF << endl;

        // extract the corresponding MJD for each state vector
        tMJD = matODECEF.col(1).array();
        // Create a matrix to store MJD and satECI values
        MatrixXd satECIData(tMJD.size(), 7);

        for (int i = 0; i < tMJD.size(); i++)
        {
            VectorXd satECI = VectorXd::Zero(6);
            VectorXd satECEF = matODECEF.row(i).tail(6).transpose();
            // cout << "satECEF" << satECEF << endl;
            ecef2eciVec_sofa(tMJD[i], iersInstance, satECEF, satECI);
            // cout << "satECI" << satECI << endl;

            // Format the output row using std::ostringstream
            ostringstream oss;
            oss << fixed << setprecision(10) << tMJD[i];
            for (int j = 0; j < satECI.size(); j++)
            {
                oss << " " << fixed << setprecision(DBL_DIG) << satECI[j];
            }

            // Assign the formatted row to satECIData
            istringstream issSat(oss.str());
            VectorXd tempRow(satECIData.cols());
            for (int j = 0; j < tempRow.size(); j++)
            {
                issSat >> tempRow[j];
            }
            satECIData.row(i) = tempRow;
        }
        // Write MJD and satECI to CSV file
        writeCSV("refdata/od_eci_id_46984.csv", satECIData, header);
    }

    // perform coordinate transformation for satellite via a CPF file
    if (flagCPF)
    {
        int headerLinesToSkip = 3;
        // read the cpf file
        MatrixXd orbCPF = readFile(orbECEFFile, headerLinesToSkip);
        MatrixXd satECIData(orbCPF.rows(), 4);
        // cout << "orbCPF" << orbCPF << endl;

        for (int i = 0; i < orbCPF.rows(); i++)
        // for (int i = 0; i < 1; i++)
        {
            VectorXd satECEF = orbCPF.row(i).segment<3>(5);
            cout << "satECEF:   " << satECEF.transpose() << endl;
            VectorXd satECI = VectorXd::Zero(3);
            double epochMJD = orbCPF(i, 2) + orbCPF(i, 3) / 86400;
            cout << "epochMJD:  " << epochMJD << endl;
            getIERS(epochMJD, erpt, iersInstance);
            ecef2eciVec_sofa(epochMJD, iersInstance, satECEF, satECI);
            cout << "satECI" << satECI.transpose() << endl;

            // Format the output row using std::ostringstream
            ostringstream oss;
            oss << fixed << setprecision(10) << epochMJD;
            for (int j = 0; j < satECI.size(); j++)
            {
                oss << " " << fixed << setprecision(DBL_DIG) << satECI[j];
            }

            // Assign the formatted row to stnECIData
            istringstream issSat(oss.str());
            VectorXd tempRow(satECIData.cols());
            for (int j = 0; j < tempRow.size(); j++)
            {
                issSat >> tempRow[j];
            }
            satECIData.row(i) = tempRow;
        }
        vector<string> header = {"MJD", "X_ECI", "Y_ECI", "Z_ECI"};
        // Write MJD and satECI to CSV file
        writeCSV("refdata/od_eci_id_46984_from_cpf.csv", satECIData, header);
    }

    return 0;
}
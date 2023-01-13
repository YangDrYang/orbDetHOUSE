#include <yaml-cpp/yaml.h>
#include <iostream>
#include <Eigen/Dense>

Eigen::VectorXd stdVec2EigenVec(const std::vector<double>& stdVec);

int main(){
    using namespace std;
    using namespace Eigen;
    YAML::Node config = YAML::LoadFile("config.yaml");


    std::vector<double> tempVec = config["orbtial_parameters"]["initial_state"].as<std::vector<double>>();
    VectorXd initialState = stdVec2EigenVec(tempVec);


    cout << initialState << endl;
}

Eigen::VectorXd stdVec2EigenVec(const std::vector<double>& stdVec){
    Eigen::VectorXd eigenVec(stdVec.size());
    for (int i = 0; i < 6; i++)
    {
        eigenVec(i) = stdVec[i];
    }
    return eigenVec;
}
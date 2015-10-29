

#ifndef ___DPanelCode__Geom__
#define ___DPanelCode__Geom__

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <string>


using namespace Eigen;

class Geom{
public:
    void getCoords(MatrixXd &coords, std::string &file_name);
    
    void findControlPts(MatrixXd &coords, MatrixXd &controlPts, VectorXd &dl);
    
    void calcPanelAngles(MatrixXd &coords, VectorXd &theta);
};


#endif /* defined(___DPanelCode__Geom__) */
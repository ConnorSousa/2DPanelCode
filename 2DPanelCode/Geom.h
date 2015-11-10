

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
    
    void findwakepans(VectorXd &wakePan1, VectorXd &wakePan2, double &c_w, const double &Qinf, MatrixXd &coords, double &dt, double &wakeTheta);

};


#endif /* defined(___DPanelCode__Geom__) */
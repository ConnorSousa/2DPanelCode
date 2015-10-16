

#ifndef ___DPanelCode__Geom__
#define ___DPanelCode__Geom__

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

class Geom{
public:
    void findControlPts(VectorXd &x, VectorXd &y, VectorXd &xControlPt, VectorXd &yControlPt, VectorXd &dl, int &n);
    
    void calcPanelAngles(VectorXd &x, VectorXd &y, VectorXd &theta, int &n);
};


#endif /* defined(___DPanelCode__Geom__) */
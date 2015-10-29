//
//  MatlabOutput.h
//  2DPanelCode
//
//  Created by Connor Sousa on 10/19/15.
//  Copyright (c) 2015 Connor Sousa. All rights reserved.
//

#ifndef ___DPanelCode__MatlabOutput__
#define ___DPanelCode__MatlabOutput__

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class MatlabOutput{
public:
    
    void print(MatrixXd &coords, const double &alpha, MatrixXd &controlPts, VectorXd &theta, int &n, MatrixXd &A, MatrixXd &b, MatrixXd &c, VectorXd &dl, VectorXd &Qtan, VectorXd &QinfTan, VectorXd &delCl,  VectorXd &mu, VectorXd &RHS,VectorXd &Cp, VectorXd &sigma);

    void GridPrint(MatrixXd &gridCoords, VectorXd &uvec, VectorXd &wvec);
};



#endif /* defined(___DPanelCode__MatlabOutput__) */

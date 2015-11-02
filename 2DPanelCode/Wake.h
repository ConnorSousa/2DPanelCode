//
//  Wake.h
//  2DPanelCode
//
//  Created by Connor Sousa on 10/19/15.
//  Copyright (c) 2015 Connor Sousa. All rights reserved.
//

#ifndef ___DPanelCode__Wake__
#define ___DPanelCode__Wake__

#include <iostream>
#include "SourceDoublet.h"

using namespace std;
using namespace Eigen;

class Wake{
public:
    
    void wake2panInfl(Vector3d &wakePan1, Vector3d &wakePan2, VectorXd &Cvec, MatrixXd &controlPts);
    
    void PartToPanInfluence(MatrixXd &particles, MatrixXd &controlPts, VectorXd &D);

    void PartInflonPart(MatrixXd &part2partVel, MatrixXd &particles);
    
    void wake2partVel(Vector3d &wakePan1, Vector3d &wakePan2, MatrixXd &particles, MatrixXd &wake2part);

    void panInflonPart(MatrixXd &coords, VectorXd &theta, VectorXd &sigma, VectorXd &mu, MatrixXd &particles);

};

#endif /* defined(___DPanelCode__Wake__) */



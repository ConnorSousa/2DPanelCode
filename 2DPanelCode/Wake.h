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
    
    void wake2influence(Vector3d bufferWake2, VectorXd Cvec, MatrixXd controlPts);
    
    void PartToPanInfluence(MatrixXd particles, MatrixXd controlPts, VectorXd D);

    
};

#endif /* defined(___DPanelCode__Wake__) */



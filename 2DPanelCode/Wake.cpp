//
//  Wake.cpp
//  2DPanelCode
//
//  Created by Connor Sousa on 10/19/15.
//  Copyright (c) 2015 Connor Sousa. All rights reserved.
//

#include <Eigen/Dense>

#include "Wake.h"
#include "SourceDoublet.h"

using namespace std;
using namespace Eigen;

void Wake::wake2influence(Vector3d bufferWake2, VectorXd Cvec, MatrixXd controlPts){
    
    // bufferWake2  = [panel start, panel end, doublet strength]
    
    int n = Cvec.rows();
    
    for(int i = 0; i < n; i++){ //n is panels including wake panel
        // Doublet influence coeff
        double xt = controlPts(i,0) - bufferWake2(0);
        double yt = controlPts(i,1) - 0;
        
        double px = xt;   //xt*cos(theta(j)) + yt*sin(theta(j));  the wake panel is horizontal, so theta = 0
        double py = yt;   //-xt*sin(theta(j))+ yt*cos(theta(j));
        double x2 = bufferWake2(1)-bufferWake2(0);
        
        double y2 = 0;
        double x1 = 0;
        double y1 = 0;
        
        double phiD;
        double mu = bufferWake2(2);
        
        
        SourceDoublet wakeobj;
        wakeobj.findPhiDoublet(mu, px, py, x1, y1, x2, y2, phiD);

        
        Cvec(i) = phiD;
    }
    
}


void Wake::PartToPanInfluence(MatrixXd particles, MatrixXd controlPts, VectorXd D){
    
    MatrixXd phiV(controlPts.rows(),particles.rows());
    
    for(int i = 0; i < controlPts.rows();i++){
        for(int j = 0; j < particles.rows(); j++){
            
            // phiV = -(gamm/2pi)atan(z-z0/x-x0). ref frame on pdf p.247
            phiV(i,j) = -(particles(j,2)/(2*3.14159265))*atan2(controlPts(i,2)-particles(j,1),controlPts(i,0)-particles(j,0)); // Eq 10.8
            
        }
    }
    
    D = phiV.rowwise().sum(); // Make sure this is summing the right thing. Treat each col/row as a different vortex and sum.
    
    
};





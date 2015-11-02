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

const double pi = 3.141592653589;

void Wake::wake2panInfl(Vector3d &wakePan1, Vector3d &wakePan2, VectorXd &Cvec, MatrixXd &controlPts){
        
    for(int i = 0; i < Cvec.rows(); i++){ //n is panels including wake panel
        
        double xt = controlPts(i,0) - wakePan1(0);
        double yt = controlPts(i,1) - 0; // Assuming that the wake panel is always straight out from the trailin edge and that any alpha will only affect the flow, not actually tilt the geometry. Might need to change this if you want to follow the trailling edge for a cambered/cusped airfoil not pointing relatively straight
        
        double px = xt;
        double py = yt;
        double x2 = wakePan1(1)-wakePan1(0);
        
        double y2 = 0;
        double x1 = 0;
        double y1 = 0;
        
        double phiD1;
        double mu = wakePan1(2);
        
        
        SourceDoublet wakeobj;
        wakeobj.findPhiDoublet(mu, px, py, x1, y1, x2, y2, phiD1);
        
        
        xt = controlPts(i,0) - wakePan2(0);
        yt = controlPts(i,1) - 0; //assuming that the wake panel is always straight out from the trailin edge and that any alpha will only affect the flow, not actually tilt the geometry. Might need to change this if you want to follow the trailling edge for a cambered/cusped airfoil not pointing relatively straight
        
        px = xt;
        py = yt;
        x2 = wakePan2(1)-wakePan2(0);
        
        y2 = 0;
        x1 = 0;
        y1 = 0;
        
        double phiD2;
        mu = wakePan2(2);
        
        wakeobj.findPhiDoublet(mu, px, py, x1, y1, x2, y2, phiD2);
        
        Cvec(i) = phiD1 + phiD2;
    }
    
}


void Wake::PartToPanInfluence(MatrixXd &particles, MatrixXd &controlPts, VectorXd &D){
    
    MatrixXd phiV(controlPts.rows(),particles.rows());
    
    for(int i = 0; i < controlPts.rows();i++){
        for(int j = 0; j < particles.rows(); j++){
            
            // phiV = -(gamm/2pi)atan(z-z0/x-x0). ref frame on pdf p.247
            phiV(i,j) = -(particles(j,2)/(2*3.14159265))*atan2(controlPts(i,2)-particles(j,1),controlPts(i,0)-particles(j,0));
            // Eq 10.8
            
        }
    }
    
    D = phiV.rowwise().sum(); // Make sure this is summing the right thing. Treat each col/row as a different vortex and sum.
    
    
};


void Wake::PartInflonPart(MatrixXd &part2partVel, MatrixXd &particles){
    
    MatrixXd upart(particles.rows(),particles.rows());
    MatrixXd wpart(particles.rows(),particles.rows());
    //x0,z0 will be the current particle(i) and x,y will be the one giving the influence (j) and will have the strength in the equation.
    
    for(int i = 0; i <particles.rows(); i++){
        for(int j=0; j<particles.rows(); j++){
            upart(i,j) = (particles(j,4)/(2*pi))*((particles(j,1)-particles(i,1))/(pow(particles(j,0)-particles(i,0),2)+pow(particles(j,1)-particles(i,1),2)));
            
            wpart(i,j) = (-particles(j,4)/(2*pi))*((particles(j,0)-particles(i,0))/(pow(particles(j,0)-particles(i,0),2)+pow(particles(j,1)-particles(i,1),2)));
        }
    }
    // My attempt at assigning the below directly. Make sure this produces the same results...
    part2partVel.col(0) = upart.rowwise().sum();
    part2partVel.col(1) = wpart.rowwise().sum();
    
    /*
     // There probably is a way to assign this directly...
     VectorXd u(particles.rows()) = upart.rowwise().sum();
     VectorXd w(particles.rows()) = wpart.rowwise().sum();
     
     for(int i=0; i < particles.rows(); i++){
     part2partVel(i,0) = u(i);
     part2partVel(i,1) = w(i);
     }
     
     */
};

void Wake::wake2partVel(Vector3d &wakePan1, Vector3d &wakePan2, MatrixXd &particles, MatrixXd &wake2part){
    
    for(int i = 0; i < particles.rows(); i++){
        
        // Converting grid and control points to panel coords.
        double xt = particles(i,0) - wakePan1(0);
        double yt = particles(i,1) - 0; // Wakepan y is zero.
        double dx = wakePan1(1) - wakePan1(0);
        double dy = 0;
        
        double x = xt;
        double z = yt;
        
        double x2 = sqrt(dx*dx + dy*dy);
        double x1 = 0;
        
        double dif = atan2(z,(x-x2)) - atan2(z,(x-x1)); //theta2-theta1
        double val;
        
        if(dif < -pi){
            val = dif + 2*pi;
        }else if(dif > pi){
            val = dif-2*pi;
        }else{
            val = dif;
        }
        
        // In book on page 236 and 281, the signs of u and w are switched.
        double uWake1 = (-wakePan1(2)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z));
        double wWake1 = (wakePan1(2)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z));
        
        //---------------------------------------Second Wake Panel-------------------------------------------------------//
        
        xt = particles(i,0) - wakePan2(0);
        yt = particles(i,1) - 0; // Wakepan y is zero.
        dx = wakePan2(1) - wakePan2(0);
        dy = 0;
        
        x = xt;
        z = yt;
        
        x2 = sqrt(dx*dx + dy*dy);
        x1 = 0;
        
        dif = atan2(z,(x-x2)) - atan2(z,(x-x1)); //theta2-theta1
        
        if(dif < -pi){
            val = dif + 2*pi;
        }else if(dif > pi){
            val = dif-2*pi;
        }else{
            val = dif;
        }
        
        // In book on page 236 and 281, the signs of u and w are switched.
        double uWake2 = (-wakePan2(2)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z));
        double wWake2 = (wakePan2(2)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z));
        
        wake2part(i,0) = uWake1 + uWake2;
        wake2part(i,1) = wWake1 + wWake2;
        
    }
};


void Wake::panInflonPart(MatrixXd &coords, VectorXd &theta, VectorXd &sigma, VectorXd &mu, MatrixXd &particles){
    
    double u_local;
    double w_local;
    
    MatrixXd u_global(particles.rows(),coords.rows());
    MatrixXd w_global(particles.rows(),coords.rows());
    
    for(int i = 0; i < particles.rows(); i++){
        for(int j = 0; j < coords.rows(); j++){ // through body panels excluding wake
            
            // Converting grid and control points to panel coords.
            double xt = particles(i,0) - coords(j,0);
            double yt = particles(i,1) - coords(j,1);
            double dx = coords(j+1,0) - coords(j,0);
            double dy = coords(j+1,1) - coords(j,1);
            
            double x = xt*cos(theta(j)) + yt*sin(theta(j));
            double z = xt*-sin(theta(j))+ yt*cos(theta(j));
            
            double x2 = sqrt(dx*dx + dy*dy);
            double x1 = 0;
            
            double dif = atan2(z,(x-x2)) - atan2(z,(x-x1)); //theta2-theta1
            double val;
            
            if(dif < -pi){
                val = dif + 2*pi;
            }else if(dif > pi){
                val = dif-2*pi;
            }else{
                val = dif;
            }
            
            // In book on page 236 and 281, the signs of u and w are switched.
            u_local = (-mu(j)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z)) + (sigma(j)/(4*pi))*log((pow(x-x1,2) + z*z)/(pow(x-x2,2) + z*z));
            w_local = (mu(j)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z)) + (sigma(j)/(2*pi))*(val);
            
            
            // Transform every velocity into global coords. Transformation 11.23 & 11.23a on pp.277. These equations look to be labeled backwards from what I derive. I use my derived version in earlier code. Example code in the appendix pp. 555 uses negative thetas which, using trig identities, is equivalent to my version.
            u_global(i,j) = u_local*cos(theta(j)) - w_local*sin(theta(j));
            w_global(i,j) = u_local*sin(theta(j)) + w_local*cos(theta(j));
        }
    }
    

    particles.col(0) = u_global.rowwise().sum();
    particles.col(1) = w_global.rowwise().sum();
    cout << "check this too. it's in panInflonPart in wake class\n";
};



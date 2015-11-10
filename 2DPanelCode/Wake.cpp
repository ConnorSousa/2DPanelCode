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

void Wake::wake2bodyInfl(VectorXd &wakePan1, VectorXd &wakePan2, VectorXd &C, MatrixXd &controlPts){
        
    for(int i = 0; i < C.rows(); i++){ //n is panels including wake panel
        
        
        double wakeTheta = atan2(wakePan1(3) - wakePan1(1),wakePan1(2) - wakePan1(0));
        

        // First panel
        // Converting grid and control points to panel coords
        double xt = controlPts(i,0) - wakePan1(0);
        double yt = controlPts(i,1) - wakePan1(1);
        double dx = wakePan1(2) - wakePan1(0);
        double dy = wakePan1(3) - wakePan1(1);
        
        double px = xt*cos(wakeTheta) + yt*sin(wakeTheta);
        double py = -xt*sin(wakeTheta)+ yt*cos(wakeTheta);
        double x2 = sqrt(dx*dx + dy*dy);
        
        double y2 = 0;
        double x1 = 0;
        double y1 = 0;

        double phiD1;
        double mu = wakePan1(4);
        
        
        SourceDoublet wakeobj;
        wakeobj.findPhiDoublet(mu, px, py, x1, y1, x2, y2, phiD1);
        
        // Second panel
         xt = controlPts(i,0) - wakePan2(0);
         yt = controlPts(i,1) - wakePan2(1);
         dx = wakePan2(2) - wakePan2(0);
         dy = wakePan2(3) - wakePan2(1);
        
         px = xt*cos(wakeTheta) + yt*sin(wakeTheta);
         py = -xt*sin(wakeTheta)+ yt*cos(wakeTheta);
         x2 = sqrt(dx*dx + dy*dy);
        
         y2 = 0;
         x1 = 0;
         y1 = 0;
        
        double phiD2;
        mu = wakePan2(4);
        
        wakeobj.findPhiDoublet(mu, px, py, x1, y1, x2, y2, phiD2);
        
        C(i) = phiD1 + phiD2;
    }
    
}


void Wake::PartToPanInfluence(MatrixXd &particles, MatrixXd &controlPts, VectorXd &D){
    
    MatrixXd phiV(controlPts.rows(),particles.rows());
    
    for(int i = 0; i < controlPts.rows();i++){
        for(int j = 0; j < particles.rows(); j++){
            
            // phiV = -(gamm/2pi)atan(z-z0/x-x0). ref frame on pdf p.247
            phiV(i,j) = -(particles(j,2)/(2*pi))*atan2(controlPts(i,1)-particles(j,1),controlPts(i,0)-particles(j,0)); // Eq 10.8
            
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
            
            if(i == j){ //particles don't influence themselves
                upart(i,j) = 0;
                wpart(i,j) = 0;
            }else{
            upart(i,j) = (particles(j,2)/(2*pi))*((particles(j,1)-particles(i,1))/(pow(particles(j,0)-particles(i,0),2)+pow(particles(j,1)-particles(i,1),2)));
            
            wpart(i,j) = (-particles(j,2)/(2*pi))*((particles(j,0)-particles(i,0))/(pow(particles(j,0)-particles(i,0),2)+pow(particles(j,1)-particles(i,1),2)));
            }
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


void Wake::body2partVel(MatrixXd &coords, VectorXd &theta, VectorXd &sigma, VectorXd &mu, MatrixXd &particles, MatrixXd &pan2partVel, const double Qinf){
    
    MatrixXd u_global(particles.rows(),coords.rows()-1);
    MatrixXd w_global(particles.rows(),coords.rows()-1);
    double u_local;
    double w_local;
    
    
    for(int i = 0; i < particles.rows(); i++){
        for(int j = 0; j < coords.rows()-1; j++){ // number of panels
            
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
            
            // In book on page 233-36
            if(abs(z) < 1.0e-9){ // if point is on panel surface, then velocity is different
                
                u_local = 0 + (sigma(j)/(2*pi))*log(abs((x-x1)/(x-x2))); // eq. 10.32 & 10.23
                w_local = (mu(j)/(2*pi))*((1/(x-x1))-(1/(x-x2))) + sigma(j)/2; //eq 10.33 & 10.24
                
            }else{
                
                u_local = (-mu(j)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z)) + (sigma(j)/(4*pi))*log((pow(x-x1,2) + z*z)/(pow(x-x2,2) + z*z)); // eq 10.29 & 10.20
                w_local = (mu(j)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z)) + (sigma(j)/(2*pi))*(val); // eq. 10.30 & 10.21
            };
            

            // Transform every velocity into global coords. Transformation 11.23 & 11.23a on pp.277. These equations look to be labeled backwards from what I derive. I use my derived version in earlier code. Example code in the appendix pp. 555 uses negative thetas which, using trig identities, is equivalent to my version.
            u_global(i,j) = u_local*cos(theta(j)) - w_local*sin(theta(j));
            w_global(i,j) = u_local*sin(theta(j)) + w_local*cos(theta(j));
        }
    }
    
    VectorXd uvec(particles.rows());
    VectorXd wvec(particles.rows());
    pan2partVel.col(0) = u_global.rowwise().sum()*Qinf;
    pan2partVel.col(1) = w_global.rowwise().sum()*Qinf;

};


void Wake::wake2partVel(VectorXd &wakePan1, VectorXd &wakePan2, MatrixXd &particles, MatrixXd &wake2part, double &wakeTheta){
    

    double u_local;
    double w_local;
        
    for(int i = 0; i < particles.rows(); i++){
        
        // First panel
        // Converting grid and control points to panel coords
        double xt = particles(i,0) - wakePan1(0);
        double yt = particles(i,1) - wakePan1(1);
        double dx = wakePan1(2) - wakePan1(0);
        double dy = wakePan1(3) - wakePan1(1);
        
        double x = xt*cos(wakeTheta) + yt*sin(wakeTheta);
        double z = xt*-sin(wakeTheta)+ yt*cos(wakeTheta);
        
        double x2 = sqrt(dx*dx + dy*dy);
        double x1 = 0;
        
        // In book on page 233-36
        if(abs(z) < 1.0e-9){ // if point is on panel surface, then velocity is different
            u_local = 0; // eq. 10.32
            w_local = (wakePan1(4)/(2*pi))*((1/(x-x1))-(1/(x-x2))); // eq 10.33
            
        }else{
            u_local = (-wakePan1(4)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z)); // eq. 10.29
            w_local = (wakePan1(4)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z)); // eq. 10.30
            
        };
        // Second Panel
        xt = particles(i,0) - wakePan2(0);
        yt = particles(i,1) - wakePan2(1);
        dx = wakePan2(2) - wakePan2(0);
        dy = wakePan2(3) - wakePan2(1);
        
        x = xt*cos(wakeTheta) + yt*sin(wakeTheta);
        z = xt*-sin(wakeTheta)+ yt*cos(wakeTheta);
        
        x2 = sqrt(dx*dx + dy*dy);
        
        if(abs(z) < 1.0e-9){ // if point is on panel surface, then velocity is different
            u_local += 0; // += adds to the current ulocal velocity
            w_local += (wakePan2(4)/(2*pi))*((1/(x-x1))-(1/(x-x2))); //eq. 10.33
            
        }else{
            u_local += (-wakePan2(4)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z)); // eq. 10.29
            w_local += (wakePan2(4)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z)); //eq. 10.30
            
        };
        
        // Transform every velocity into global coords. Transformation 11.23 & 11.23a on pp.277. These equations look to be labeled backwards from what I derive. I use my derived version in earlier code. Example code in the appendix pp. 555 uses negative thetas which, using trig identities, is equivalent to my version.
        wake2part(0) = u_local*cos(wakeTheta) - w_local*sin(wakeTheta);
        wake2part(1) = u_local*sin(wakeTheta) + w_local*cos(wakeTheta);
    }
    
};





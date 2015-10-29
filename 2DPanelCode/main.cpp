//
//  main.cpp
//  2DPanelCode
//
//  Created by Connor Sousa on 4/20/15.
//  Copyright (c) 2015 Connor Sousa. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "Geom.h"
#include "SourceDoublet.h"
#include "MatlabOutput.h"
#include "Wake.h"

using namespace std;
using namespace Eigen;
const double pi = 3.141592653589;


int main()
{
    const double Qinf = 1;
    const double alpha = 0/57.2958;
    MatrixXd coords(1,2);
    
    
    Geom object;
    string file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/4412.txt";
    object.getCoords(coords, file_name);
    int n = coords.rows()-1;

    // Find control points, panel lengths, and panel angles
    MatrixXd controlPts(n,2);
    VectorXd dl(n);
    VectorXd theta(n);

    object.findControlPts(coords, controlPts, dl);
    object.calcPanelAngles(coords, theta);
    
    // Find source strength
    VectorXd sigma(n);
    SourceDoublet object2;
    object2.findSigma(alpha, theta, sigma, Qinf);
    
    
    // Find Doublet coefficient matrix
    MatrixXd A(n,n);
    MatrixXd c(n+1,n+1);
    MatrixXd b(n,n);

    object2.InfluenceCoeffs(coords, controlPts, A, c, theta, b, dl);
    
    // Find right had side vector and solve
    VectorXd RHS(n);
    RHS = -b*sigma;
    VectorXd mu = A.colPivHouseholderQr().solve(RHS);


// ---------------------------------------------Wake---------------------------------------------//
    //First time through. No loop yet.

    // RHS2 = -b*sigma - Cvec*wakedoublet - d*vortex parts;
    // RHS2 =   RHS    -        C         -       D
/*
    
    
    // Finding strength and size of wake doublet panel.
    double c_w = 0.3;
    double dt = .1;
    int timesteps = 2;
    int n_part = timesteps-1;
    
    Vector3d bufferWake2; // Vector containing [x1,x2,mu].
    bufferWake2(0) = coords(0,0) + c_w*Qinf*dt; // trailing edge plus the length of the first panel.
    bufferWake2(1) = bufferWake2(0) + Qinf*dt; // x2
    bufferWake2(2) = mu(n-1)-mu(0); // Eq. 11.36. mu1 - mu_n + mu_w = 0
    
    
    
    // Influence of wake panel onto body
    VectorXd C(n);
    Wake time2;
    time2.wake2influence(bufferWake2, C, controlPts);

    
    // Particle strength
    MatrixXd particles(n_part,3); // (x,y,strength)
    particles(0,0) = bufferWake2(0);
    particles(0,1) = bufferWake2(1);
    particles(0,2) = bufferWake2(2);

    
    
    
    // Velocity influence on particles from particles
    //      U_inf + part_to_part + pan_to_part       
    
    
    
    // Body velocity influence on particles ===> From the velocity infl. written before
    MatrixXd bodyInflonParts(n_part,2);
    
    
    
    
    
    
    
    // Potential influence from particles onto body
    VectorXd D(n); // The vector with the sum of particle strengths on each panel.
    time2.PartToPanInfluence(particles, controlPts, D);
    
    
    
    
    
    
    
    // time2.part_strength //influence from body as well as vortex damping.
    // time2.trackpart
    
    
    //time2.part influence
    

    VectorXd RHS2(n);
   // RHS2 = -b*sigma - c*wakedoublet - d*vortex parts;
*/
    
 


   //---------------------------------------------------------------------------------------------------------
    
    // Calc tangential Q_inf
    VectorXd QinfTan(n);
    for(int i =0 ; i < n; i++){
        QinfTan(i) = Qinf*(cos(alpha)*cos(theta(i))+sin(alpha)*sin(theta(i)));
    }
    // Calc tangential velocity
    VectorXd Qtan(n);
    for(int i =0 ; i < n-1; i++){
        Qtan(i) = ((mu(i)-mu(i+1))/dl(i)*Qinf + QinfTan(i));
    }
    Qtan(n-1) = ((mu(n-2)-mu(n-1))/dl(n-1)*Qinf + QinfTan(n-1)); // Backwards differencing the nth panel is same as forward differencing the nth-1 panel.
    

    
    // Calc Cp
    VectorXd Cp(n);
    for (int i = 0; i < n; i++){
        Cp(i) = 1- Qtan(i)*Qtan(i)/(Qinf*Qinf);
    }
    
    // Calc Coefficients
    VectorXd delCl(n);
    VectorXd delCd(n);
    VectorXd delCm(n);
    
    for (int i = 0; i < n; i++){
        delCl(i) = -Cp(i)*dl(i)*cos(theta(i))/1;
        delCd(i) = -Cp(i)*dl(i)*sin(theta(i))/1;
        delCm(i) = -delCl(i)*(controlPts(i,0)-.25) + delCd(i)*controlPts(i,1);
        
    }
    
    cout << "Cl = " << delCl.sum() << endl;
    cout << "Cd = " << delCd.sum() << endl;
    // cout <<"Cm25 = "<< delCm.sum() << endl;


//--------------------------------------------Grid calcs------------------------------------------//

//file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/SixteenGridSmall.txt";
file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/NearField.txt";

    
    MatrixXd gridCoords(1,2);
    Geom gobject;
    gobject.getCoords(gridCoords, file_name);
    
    MatrixXd u_global(gridCoords.rows(),n+1);
    MatrixXd w_global(gridCoords.rows(),n+1);
    double u_local;
    double w_local;
    
    n = n+1; //now including wake panel.
    theta.conservativeResize(n);
    theta(n-1) = 0; //Wake panel is flat and fixed for now
    mu.conservativeResize(n);
    mu(n-1) = mu(n-2)-mu(0); // Kutta condition eq. 11.36
    
    
    for(int i = 0; i < gridCoords.rows(); i++){
        for(int j = 0; j < n; j++){ //n is panels including wake panel
            
            if (j < n-1){
                // Converting grid and control points to panel coords.
                double xt = gridCoords(i,0) - coords(j,0);
                double yt = gridCoords(i,1) - coords(j,1);
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

            }else{

                // Add wake panel
                double xt = gridCoords(i,0) - coords(j,0);
                double yt = gridCoords(i,1) - coords(j,1);
                
                double x = xt*cos(theta(j)) + yt*sin(theta(j));
                double z = xt*sin(theta(j))+ yt*cos(theta(j));
                
                double x1 = 0;
                
                u_local = (-mu(j)/(2*pi))*(z/(pow(x-x1,2) + z*z));
                w_local = (mu(j)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z));

            }
            
// Transform every velocity into global coords. Transformation 11.23 & 11.23a on pp.277. These equations look to be labeled backwards from what I derive. I use my derived version in earlier code. Example code in the appendix pp. 555 uses negative thetas which, using trig identities, is equivalent to my version.
            u_global(i,j) = u_local*cos(theta(j)) - w_local*sin(theta(j));
            w_global(i,j) = u_local*sin(theta(j)) + w_local*cos(theta(j));
        }
    }
    
    VectorXd uvec(gridCoords.rows());
    VectorXd wvec(gridCoords.rows());
    uvec = u_global.rowwise().sum() + MatrixXd::Ones(gridCoords.rows(),1)*Qinf ;
    wvec = w_global.rowwise().sum();
    


 
    
    
    
    
    
MatlabOutput object3;
object3.print(coords,alpha, controlPts, theta, n, A, b, c, dl, Qtan, QinfTan, delCl, mu, RHS, Cp, sigma);
object3.GridPrint(gridCoords, uvec, wvec);
    
    

    
    
cout << "%";
}










/*
 
 //--------------------------------------------panel test------------------------------------------//
 
 file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/NineHundredGrid.txt";
 //file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/NearField.txt";
 
 
 
 MatrixXd gridCoords(1,2);
 Geom gobject;
 gobject.getCoords(gridCoords, file_name);
 
 MatrixXd u_global(gridCoords.rows(),n+1);
 MatrixXd w_global(gridCoords.rows(),n+1);
 double u_local;
 double w_local;
 
 n = 1; //now including the last panel.
 
 double thet = 0;
 double mu2 = 1;
 double sig2 = 0;
 
 
 
 for(int i = 0; i < gridCoords.rows(); i++){
 for(int j = 0; j < n; j++){ //n is panels including wake panel
 
 
 // Converting grid and control points to panel coords.
 double xt = gridCoords(i,0) - 0;
 double yt = gridCoords(i,1) - 0;
 double dx = 1 - 0;
 double dy = 0;
 
 double x = xt*cos(thet) + yt*sin(thet);
 double z = xt*-sin(thet)+ yt*cos(thet);
 
 double x2 = sqrt(dx*dx + dy*dy);
 double x1 = 0;
 
 double dif = atan2(z,(x-x2)) - atan2(z,(x-x1));
 double val;
 
 if(dif < -pi){
 val = dif + 2*pi;
 }else if(dif > pi){
 val = dif-2*pi;
 }else{
 val = dif;
 }
 
 
 // Book is confusing. on page 236 and 281, the signs of u and w are switched?? eq. 10.29
 u_local = (-mu2/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z)) + (sig2/(4*pi))*log((pow(x-x1,2) + z*z)/(pow(x-x2,2) + z*z));
 
 w_local = (mu2/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z)) + (sig2/(2*pi))*(val);
 
 
 
 // Transform every velocity into global coords. Transformation 11.23 & 11.23a on pp.277. These equations look to be labeled backwards from what I derive. I use my derived version in earlier code. Example code in the appendix pp. 555 uses negative thetas which, using trig identities, is equivalent to my version.
 u_global(i,j) = u_local*cos(thet) - w_local*sin(thet);
 w_global(i,j) = u_local*sin(thet) + w_local*cos(thet);
 }
 }
 
 VectorXd uvec(gridCoords.rows());
 VectorXd wvec(gridCoords.rows());
 uvec = u_global.rowwise().sum();
 wvec = w_global.rowwise().sum();

 
 //-----------------------------------------------------end---------------------------------------//
 
 
 
 
 cout << "u_global = [";
 for(int i = 0; i < u_global.rows(); i++){
 for (int j = 0; j < u_global.cols(); j++){
 cout << u_global(i,j) << " ";
 }
 if(i < u_global.rows()){
 cout << ";";
 }
 }
 cout << "];" << endl;
 
 
 cout << "w_global = [";
 for(int i = 0; i < w_global.rows(); i++){
 for (int j = 0; j < w_global.cols(); j++){
 cout << w_global(i,j) << " ";
 }
 if(i < w_global.rows()){
 cout << ";";
 }
 }
 cout << "];" << endl;
 
 





 VectorXd sigInfUtot(gridCoords.rows());
 VectorXd sigInfWtot(gridCoords.rows());

 sigInfUtot = sigInfU.rowwise().sum();
 sigInfWtot = sigInfW.rowwise().sum();
 
 
 
 cout << "SourceInf = [";
 for(int i = 0; i< gridCoords.rows(); i++){
 cout << sigInfUtot(i) << "," << sigInfUtot(i);
 if(i < gridCoords.rows()){
 cout << ";";
 }
 }
 cout << "];" <<endl;
 

 */


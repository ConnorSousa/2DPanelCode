//
//  main.cpp
//  2DPanelCode
//
//  Created by Connor Sousa on 4/20/15.
//  Copyright (c) 2015 Connor Sousa. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "Geom.h"
#include "SourceDoublet.h"
#include "MatlabOutput.h"
#include "Wake.h"
#include <iomanip> // to set output precision

using namespace std;
using namespace Eigen;
const double pi = 3.14159265;

int main()
{
    const double Qinf = 1;
    const double alpha = 5/57.2958;
    MatrixXd coords(1,2);
    
    Geom geomObj;
    string file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/4412.txt";
    geomObj.getCoords(coords, file_name);
    int n = coords.rows()-1;

    // Find control points, panel lengths, and panel angles
    MatrixXd controlPts(n,2);
    VectorXd dl(n);
    VectorXd theta(n);

    geomObj.findControlPts(coords, controlPts, dl);
    geomObj.calcPanelAngles(coords, theta);
    
    // Find source strength
    VectorXd sigma(n);
    SourceDoublet sdobj;
    sdobj.findSigma(alpha, theta, sigma, Qinf);
    
    
    // Find Doublet coefficient matrix
    MatrixXd A(n,n);
    MatrixXd c(n+1,n+1);
    MatrixXd b(n,n);

    sdobj.InfluenceCoeffs(coords, controlPts, A, c, theta, b, dl);
    
    // Find right had side vector and solve
    VectorXd RHS(n);
    RHS = -b*sigma;
    VectorXd mu = A.colPivHouseholderQr().solve(RHS);

//===================================================== Creating Wake Panels ===========================================================

    //   RHS2 = -b*sigma - Cvec*BufferWake;
    //   RHS2 =  (1)RHS  -      (2)C
 
   // Find the wake panel geom
    double c_w = 0.3;
    double dt = .1;
    double wakeTheta;
    VectorXd wakePan1(5); // Vectors containing [x1,y1,x2,y2,mu]
    VectorXd wakePan2(5); // [x1,y1,x2,y2,mu]
    wakePan2(4) = mu(n-1)-mu(0); // Wake panel 2 is 1 from previous timestep
    geomObj.findwakepans(wakePan1, wakePan2, c_w,  Qinf,  coords, dt, wakeTheta);
    
    
   // Find new A matrix because panel length has changed
    sdobj.doubletMatA(controlPts, A, c, wakePan1);
    
    // Influence of wake panel onto body
    VectorXd C(n);
    Wake wakeobj;
    wakeobj.wake2bodyInfl(wakePan1, wakePan2, C, controlPts);

    VectorXd RHS2(n);
    RHS2 = -b*sigma - C;

    VectorXd mu2(RHS2.rows());
    mu2 = A.colPivHouseholderQr().solve(RHS2);
//this includes the influence from the second wake panel and then enforces it by the first.
    wakePan1(4) = mu2(n-1)-mu2(0);
    
    //mu = mu2;
 
//===================================================== Plus First Particle ===========================================================
    
    
   // Calculate strength of first particle and shed it into flow using the velocity at the TE. This steps everthing by dt
    int n_part = 1;
    MatrixXd particles(n_part,3); // (x,y,strength)
    particles(n_part-1,0) = wakePan2(2); // new particle goes at the end of the matrix so no data shifting is required.
    particles(n_part-1,1) = wakePan2(3);
    particles(n_part-1,2) = wakePan2(4);

    // Convect first particle downstream using wake.
    particles(0,0) += Qinf*cos(wakeTheta)*dt;
    particles(0,1) += Qinf*sin(wakeTheta)*dt;

    VectorXd D(n);
    // Potential influence from particles onto body
    wakeobj.PartToPanInfluence(particles, controlPts, D);
    
    // Influence of wake panel onto body
    wakeobj.wake2bodyInfl(wakePan1, wakePan2, C, controlPts);
    
    // Saving wake panel strength for next time step
    wakePan2(4) = wakePan1(4);
//========================================================= Define stuff =============================================================

    VectorXd RHS3(n); VectorXd mu3(n); MatrixXd UinfVel(n_part,2); MatrixXd pan2partVel(n_part,2);
    MatrixXd wake2part(n_part,2); MatrixXd part2partVel(n_part,2); MatrixXd partVel(n_part,2);

//======================================================== Start the loop =============================================================
    
    // RHS3 =  -b*sigma  -   C*BufferWake  -   d*vortex parts; // source strengths depend only on geom which doesn't change.
    // RHS3 =  (1)RHS3   -      (2)C       -       (3)D

    cout << "close all; figure;";
    for(int i = 0; i < 35; i++){

        // Find body panel strengths and first wake panel
        RHS3 = -b*sigma - C - D;
       // mu3 = A.colPivHouseholderQr().solve(RHS3);
        wakePan1(4) = mu(n-1)-mu(0);
        
        
       // Velocity influence on particles
        // U_inf + pan_to_part + wake_to_part + part_to_part
        
        // Freestream influence including alpha. First particle is influenced from last panel velocity.
        UinfVel.resize(n_part,2);
        UinfVel.col(0) = Qinf*VectorXd::Ones(n_part)*cos(alpha);
        UinfVel.col(1) = Qinf*VectorXd::Ones(n_part)*sin(alpha);
        
        
        // Body velocity influence on particles
        pan2partVel.resize(n_part,2);
        wakeobj.body2partVel(controlPts, theta, sigma, mu, particles,pan2partVel, Qinf);
        
        
        // Wake to particle velocity
        wake2part.resize(n_part,2);
        wakeobj.wake2partVel(wakePan1, wakePan2, particles, wake2part, wakeTheta);
        
        
        // Particle to Particle
        part2partVel.resize(n_part,2); // Matrix containing [du0,dw0;du1,dw1]
        wakeobj.PartInflonPart(part2partVel, particles);
        
        
        // Add these influences
        partVel.resize(n_part,2);
        partVel = UinfVel + pan2partVel + wake2part + part2partVel;
    

        // Convect particles downstream
        particles.col(0) += partVel.col(0)*dt;
        particles.col(1) += partVel.col(1)*dt;
        
        
        // New particle
        n_part++;
        particles.conservativeResize(n_part,3);
        particles(n_part-1,0) = wakePan2(2); // New particle goes at the end of the matrix to eliminate need for data shifting
        particles(n_part-1,1) = wakePan2(3);
        particles(n_part-1,2) = wakePan2(4);
        
        
        // Convect new particle downstream using wake vel.
        particles(n_part-1,0) += Qinf*cos(wakeTheta)*dt;
        particles(n_part-1,1) += Qinf*sin(wakeTheta)*dt;
        
        
        // Matlab visual
        cout << "clf; hold on; plot(coords(:,1),coords(:,2));" << endl;
        for(int z=0; z<n_part; z++){
            cout << "plot(" << particles(z,0) << "," << particles(z,1) << ",'r*');";
        }
        cout << "hold off;\n title([num2str(" << n_part << ") ' Particles']); axis([0 3.5 -1.5 1.5]); drawnow; pause(1/30);" << endl; //;
        
        
        // Potential influence from particles onto body
        wakeobj.PartToPanInfluence(particles, controlPts, D);
        
        
        // Influence of wake panel onto body
        wakeobj.wake2bodyInfl(wakePan1, wakePan2, C, controlPts);
        
        
        // Saving wake panel strength for next time step
        wakePan2(4) = wakePan1(4);
    }
   // mu = mu3;
 
   //------------------------------------------------POST PROCESS--------------------------------------------------------------
    
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
    Cp = 1 -Qtan.array().square()/(Qinf*Qinf);
    
    // Calc Coefficients
    VectorXd delCl(n);
    VectorXd delCd(n);
    VectorXd delCm(n);
    
    for (int i = 0; i < n; i++){
        delCl(i) = -Cp(i)*dl(i)*cos(theta(i));
        delCd(i) = -Cp(i)*dl(i)*sin(theta(i));
        delCm(i) = -delCl(i)*(controlPts(i,0)-.25) + delCd(i)*controlPts(i,1);
    }
    
    cout << "Cl = " << delCl.sum() << endl;
    cout << "Cd = " << delCd.sum() << endl;
    // cout <<"Cm25 = "<< delCm.sum() << endl;


//--------------------------------------------Grid calcs------------------------------------------//

   //file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/SixteenGridSmall.txt";
    file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/NearField.txt";
   // file_name = "/Users/C_Man/Desktop/Thesis/2DPanelCode/4412longnear.txt";

 
    MatrixXd gridCoords(1,2); Geom gobject; gobject.getCoords(gridCoords, file_name);
    
   // MatrixXd gridCoords(controlPts.rows(),2); gridCoords = controlPts;
    
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
                double zt = gridCoords(i,1) - coords(j,1);
                double dx = coords(j+1,0) - coords(j,0);
                double dz = coords(j+1,1) - coords(j,1);
                
                double x = xt*cos(theta(j)) + zt*sin(theta(j));
                double z = xt*-sin(theta(j))+ zt*cos(theta(j));
                
                double x2 = sqrt(dx*dx + dz*dz);
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
                  
                    u_local = 0 + (sigma(j)/(4*pi))*log(abs((x-x1)/(x-x2)));
                    w_local = (mu(j)/(2*pi))*((1/(x-x1))-(1/(x-x2))) + sigma(j)/2;
                    
                    cout << "z = " << z << ";\n";
                    
                }else{
                    u_local = (-mu(j)/(2*pi))*(z/(pow(x-x1,2) + z*z) - z/(pow(x-x2,2) + z*z)) + (sigma(j)/(4*pi))*log((pow(x-x1,2) + z*z)/(pow(x-x2,2) + z*z));
                    w_local = (mu(j)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z) - (x-x2)/(pow(x-x2,2) + z*z)) + (sigma(j)/(2*pi))*(val);
                    
                };
                
                
            }else{

                // Add wake panel
                double xt = gridCoords(i,0) - coords(j,0);
                double zt = gridCoords(i,1) - coords(j,1);
                
                double x = xt;
                double z = zt;

                double x1 = 0;
                double x2 = 67998;
                
                if( abs(z) < 1.0e-9){ // if point is on panel surface, then velocity is different
                    u_local = 0; // Eq. 10.32
                    w_local = (mu(j)/(2*pi))*(1/(x-x1) - 1/(x-x2)); // Eq. 10.33
                }else{
                    u_local = (-mu(j)/(2*pi))*(z/(pow(x-x1,2) + z*z)); // Eq. 10.29
                    w_local = (mu(j)/(2*pi))*((x-x1)/(pow(x-x1,2) + z*z)); // Eq. 10.30
                    
                };
            }
            
// Transform every velocity into global coords. Transformation 11.23 & 11.23a on pp.277. These equations look to be labeled backwards from what I derive. I use my derived version in earlier code. Example code in the appendix pp. 555 uses negative thetas which, using trig identities, is equivalent to my version.
            u_global(i,j) = u_local*cos(theta(j)) - w_local*sin(theta(j));
            w_global(i,j) = u_local*sin(theta(j)) + w_local*cos(theta(j));
        }
    }
    
    VectorXd uvec(gridCoords.rows());
    VectorXd wvec(gridCoords.rows());
    uvec = u_global.rowwise().sum()*Qinf + VectorXd::Ones(gridCoords.rows())*Qinf*cos(alpha);
    wvec = w_global.rowwise().sum()*Qinf + VectorXd::Ones(gridCoords.rows())*Qinf*sin(alpha);
    


    
  
    MatlabOutput object3;
    object3.print(coords,alpha, controlPts, theta, n, A, b, c, dl, Qtan, QinfTan, delCl, mu, RHS, Cp, sigma);
    object3.WakePrint(wakePan1, wakePan2, wakeTheta, mu2, C, RHS2, particles, UinfVel, part2partVel, pan2partVel, partVel, D, RHS3, mu3);
    object3.GridPrint(gridCoords, uvec, wvec);
    
cout << "%";
}
/*
cout << "Uinfvel" << i << "(1) = " << UinfVel(i,0) << ";";
cout << "Uinfvel" << i << "(2) = " << UinfVel(i,1) << ";";
cout << "part2part" << i << "(1) = " << part2partVel(i,0) << ";";
cout << "part2part" << i << "(2) = " << part2partVel(i,1) << ";";
cout << "partvel" << i << "(1) = " << partVel(i,0) << ";";
cout << "partvel" << i << "(2) = " << partVel(i,1) << ";";
cout << "wake2part" << i << "(1) = " << wake2part(i,0) << ";";
cout << "wake2part" << i << "(2) = " << wake2part(i,1) << ";";
cout << "pan2part" << i << "(1) = " << pan2partVel(i,0) << ";";
cout << "pan2part" << i << "(2) = " << pan2partVel(i,1) << ";";
*/
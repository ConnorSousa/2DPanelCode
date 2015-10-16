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

using namespace std;
using namespace Eigen;
const double pi = 3.14159265;


void MATLABoutput(VectorXd &x, VectorXd &y,const double &alpha, VectorXd &xControlPt, VectorXd &yControlPt, VectorXd &theta, int &numpts, int &n, MatrixXd &A, MatrixXd &b, MatrixXd &c, VectorXd &dl,VectorXd &Qtan, VectorXd &QinfTan, VectorXd &delCl,  VectorXd &doublets, VectorXd &RHS,VectorXd &Cp, VectorXd &sigma);

//================================================================================

void getCoords(VectorXd &x, VectorXd &y, int &numpts){
    ifstream thefile;
    thefile.open("/Users/C_Man/Desktop/Thesis/2DPanelCode/4412.txt");
    //thefile.open("/Users/C_Man/Desktop/Thesis/2DPanelCode/11pts.txt");
    
    int fline;
    thefile >> fline;
    numpts = fline;
    
    double xpt;
    double ypt;
    
    int count = 0;
    while(thefile >> xpt >> ypt){
        x[count] = xpt;
        y[count] = ypt;
        count++;
    }
    thefile.close();
}

void findSigma(const double &AoA, VectorXd &theta, VectorXd &sigma, int &n, double Qinf){
    for(int i = 0; i < n; i++){
        sigma[i] = Qinf*(cos(AoA)*sin(theta[i]) - sin(AoA)*cos(theta[i]));
    }
}

void findPhiSource(double &sigma, double &px, double &py, double &x1, double &y1, double &x2, double &y2, double &phiS ){
    
    phiS = (sigma/(4*pi))*((px-x1)*log(pow(px-x1,2)+pow(py,2)) - (px-x2)*log(pow(px-x2,2) + pow(py,2)) + 2*py*(atan2(py,(px-x2)) - atan2(py,(px-x1))));
}

void findPhiDoublet(double &mu, double &px, double &py, double &x1, double &y1, double &x2, double &y2, double &phiD){
    
    phiD = (-mu/(2*pi))*(atan2(py,(px-x2)) - atan2(py,(px-x1)));
}

void InfluenceCoeffs(VectorXd &x, VectorXd &y, VectorXd &xControlPt, VectorXd &yControlPt, MatrixXd &A, int &n, MatrixXd &c, VectorXd &theta, MatrixXd &b,VectorXd &dl){
    
    
    // Populate matrix
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            
            // Doublet influence coeff
            double phiD;
            double phiS;
            double mu = 1;
            double sig = 1;
            
            // Converting everything to local coords
            double xt = xControlPt(i) - x(j);
            double yt = yControlPt(i) - y(j);
            double dx= x(j+1) - x(j);
            double dy= y(j+1) - y(j);
            
            double px = xt*cos(theta(j)) + yt*sin(theta(j));
            double py = -xt*sin(theta(j))+ yt*cos(theta(j));
            double x2 = sqrt(dx*dx + dy*dy);
            // dl(j) = x2;
            
            double y2 = 0;
            double x1 = 0;
            double y1 = 0;
            
            
            
            if(i==j){
                c(i,i) = 0.5;
            }else{
                findPhiDoublet(mu,px,py,x1,y1,x2,y2,phiD);
                c(i,j) = phiD;
            }
            
            findPhiSource(sig,px,py,x1,y1,x2,y2,phiS);
            b(i,j) = phiS;
            
        }
        
        // Add wake influence coeff.
        double mu = 1;
        double phiD;
        
        // Converting everything to local coords
        double xt = xControlPt(i) - x(n);
        double yt = yControlPt(i) - y(n);
        
        double px = xt*cos(0) + yt*sin(0);
        double py = -xt*sin(0)+ yt*cos(0);
        double x2 = numeric_limits<double>::infinity();
        double y2 = 0;
        double x1 = 0;
        double y1 = 0;
        
        findPhiDoublet(mu,px,py,x1,y1,x2,y2,phiD);
        c(i,n) = phiD; //This is returning same value as book func
        
    }
    // new wake panel takes out all circulation? or lift?
    // Explicit Kutta condition: Since the condition requires circulation at TE = 0, a wake panel is added which extends to x = inf.  (mu_1 - mu_N) + mu_w = 0. This is included in the numpans+1 column because the wake doublet is added onto the numpans+1 row.
    c(n,0) = 1;
    c(n,n-1) = -1;
    c(n, n) = 1;
    
    // Reducing the c matrix from an n+1 by n+1 matrix to an n by n matrix called A.
    //     a(i,1) = c(i,1) - c(i,n+1)
    //     a(i,n) = c(i,n) + c(i,n+1)
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(j == 0){
                A(i,0) = c(i,0) - c(i,n);
            }else if(j == n-1){
                A(i,n-1) = c(i,n-1) + c(i,n);
            }else{
                A(i,j) = c(i,j);
            }
        }
    }
    
}


int main()
{
    const double Qinf = 10;
    const double alpha = 5/57.2958;
    int numpts = 161;
    VectorXd x(numpts);
    VectorXd y(numpts);
    getCoords(x, y, numpts);
    int n;
    n = numpts-1;
    
    
    // Find control points, panel lengths, and panel angles
    VectorXd xControlPt(n);
    VectorXd yControlPt(n);
    VectorXd dl(n);
    VectorXd theta(n);

    Geom object;
    object.findControlPts(x, y, xControlPt, yControlPt, dl, n);
    object.calcPanelAngles(x, y, theta, n);
    
    // Find source strength
    VectorXd sigma(n);
    findSigma(alpha, theta, sigma, n, Qinf);
    
    
    // Find Doublet coefficient matrix
    MatrixXd A(n,n);
    MatrixXd c(n+1,n+1);
    MatrixXd b(n,n);
    
    InfluenceCoeffs(x, y, xControlPt, yControlPt, A, n, c, theta, b, dl);
    
    // Find right had side vector and solve
    VectorXd RHS(n);
    RHS = -b*sigma;
    VectorXd mu = A.colPivHouseholderQr().solve(RHS);
    
    
    // Calc tangential Q_inf
    VectorXd QinfTan(n);
    for(int i =0 ; i < n; i++){
        QinfTan(i) = Qinf*(cos(alpha)*cos(theta(i))+sin(alpha)*sin(theta(i)));
    }
    
    // Calc tangential velocity
    VectorXd Qtan(n);
    for(int i =0 ; i < n-1; i++){
        Qtan(i) = (mu(i)-mu(i+1))/dl(i) + QinfTan(i);
    }
    Qtan(n-1) = (mu(n-2)-mu(n-1))/dl(n-1) + QinfTan(n-1); // Backwards differencing the nth panel is same as forward differencing the nth-1 panel.
    
    
    // Calc Cp
    VectorXd Cp(n);
    for (int i = 0; i < n; i++){
        Cp(i) = 1- Qtan(i)*Qtan(i);
    }
    
    
    // Calc Coefficients
    VectorXd delCl(n);
    VectorXd delCd(n);
    VectorXd delCm(n);
    
    for (int i = 0; i < n; i++){
        delCl(i) = -Cp(i)*dl(i)*cos(theta(i))/1;
        delCd(i) = -Cp(i)*dl(i)*sin(theta(i))/1;
        delCm(i) = -delCl(i)*(xControlPt(i)-.25) + delCd(i)*yControlPt(i);
        
    }
    
    cout << "Cl = " << delCl.sum() << endl;
    cout << "Cd = " << delCd.sum() << endl;
    // cout <<"Cm25 = "<< delCm.sum() << endl;
    
    
    
    
    /////////////////////////////////////////////////// Grid calcs
    
    ifstream myfile;
    //myfile.open("/Users/C_Man/Desktop/Thesis/2DPanelCode/SixteenGridSmall.txt");
    //myfile.open("/Users/C_Man/Desktop/Thesis/2DPanelCode/NineHundredGrid.txt");
    myfile.open("/Users/C_Man/Desktop/Thesis/2DPanelCode/NearField.txt");
    //myfile.open("/Users/C_Man/Desktop/Thesis/2DPanelCode/4412.txt");
    
    
    int f1line;
    int num;
    myfile >> f1line;
    num = f1line;
    
    double xptg;
    double yptg;
    VectorXd xg(num);
    VectorXd yg(num);
    
    
    int count = 0;
    while(myfile >> xptg >> yptg){
        xg[count] = xptg;
        yg[count] = yptg;
        count++;
    }
    
    myfile.close();
    
    
    MatrixXd u_p(num,n+1);
    MatrixXd w_p(num,n+1);
    MatrixXd uvec(num,n+1);
    MatrixXd wvec(num,n+1);
    MatrixXd uvectemp(num,n+1);
    MatrixXd wvectemp(num,n+1);

    
    // Wake Panel Addition
    mu.conservativeResize(n+1);
    mu(n) = mu(n-1)-mu(0);
    
    // calc strengths for a random point
    for(int i = 0; i < num; i++){
        for(int j = 0; j < n+1; j++){
            double xt; double yt; double dx; double dy; double px; double py; double x2; double y2; double x1; double y1;

            if (j < n){
                // Converting everything to local coords
                xt = xg(i) - x(j);
                yt = yg(i) - y(j);
                dx = x(j+1) - x(j);
                dy = y(j+1) - y(j);
                
                px = xt*cos(theta(j)) + yt*sin(theta(j));
                py = -xt*sin(theta(j))+ yt*cos(theta(j));
                
                x2 = sqrt(dx*dx + dy*dy);
                y2 = 0;
                x1 = 0;
                y1 = 0;
                
                u_p(i,j) = (-mu(j)/(2*pi))*(py/(pow(px-x1,2)+py*py) - py/(pow(px-x2,2)+py*py)) + (sigma(j)/(4*pi))*log((pow(px-x1,2)+py*py)/(pow(px-x2,2)+py*py));
                
                w_p(i,j) = (mu(j)/(2*pi))*((px-x1)/(pow(px-x1,2)+py*py) - (px-x2)/(pow(px-x2,2) + py*py)) + sigma(j)/(2*pi)*(atan2(py,(px-x2))-atan2(py,px-x1));
            }else{
                xt = xg(i) - x(j);
                yt = yg(i) - y(j);
                dx=  numeric_limits<double>::infinity() - x(j);
                dy= 0 - y(j);
                count = count+1;
                px = xt*cos(0) + yt*sin(0);
                py = -xt*sin(0)+ yt*cos(0);
                
                x2 = sqrt(dx*dx + dy*dy);
                y2 = 0;
                x1 = 0;
                y1 = 0;
                //cout << mu(j) << "  " << py << "  " << px << "  " << x1 << "  " << x2 << endl;
                //Wake panel is only doublet. I think
                u_p(i,j) = (-mu(j)/(2*pi))*(py/(pow(px-x1,2)+py*py));
                w_p(i,j) = (mu(j)/(2*pi))*((px-x1)/(pow(px-x1,2)+py*py));
            }

            uvectemp(i,j) = u_p(i,j)*(cos(alpha) - sin(alpha));
            wvectemp(i,j) = w_p(i,j)*(sin(alpha) + cos(alpha));

        }
    }
    
    for(int i = 0; i< num;i++){         //sum the influences on point from each panel
        uvec(i) = uvectemp.row(i).sum();
        wvec(i) = wvectemp.row(i).sum();
    }

    cout << "VelGrid = [";
    for(int i = 0; i<num; i++){
        cout << xg(i) << "," << yg(i) << "," << uvec(i) << "," << wvec(i);
        if(i < num-1){
            cout << ";";
        }
    }
    cout << "];%";
  
    
    
    
    
    
    
    //----------------------------------------WAKE PANEL----------------------------------------------//
    
    
    double WakePanelc = 0.3;
    double TimeStep = 1;
    double WakePanelLength = WakePanelc*Qinf*TimeStep;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //MATLAB output
    MATLABoutput(x,y,alpha,xControlPt,yControlPt,theta,numpts,n,A,b,c,dl,Qtan,QinfTan,delCl,mu,RHS,Cp,sigma);
    
}




void MATLABoutput(VectorXd &x, VectorXd &y, const double &alpha, VectorXd &xControlPt, VectorXd &yControlPt, VectorXd &theta, int &numpts, int &n, MatrixXd &A, MatrixXd &b, MatrixXd &c, VectorXd &dl, VectorXd &Qtan, VectorXd &QinfTan, VectorXd &delCl,  VectorXd &mu, VectorXd &RHS,VectorXd &Cp, VectorXd &sigma){
    
    cout << /*"clear all; close all;*/ "format compact;" << endl;
    cout << "alpha = " << alpha*57.2958 << ";"<< endl;
    cout << "x = [";
    
    for(int i = 0; i < n+1; i++){
        cout << x[i] << ",";
    }
    cout << "];" << endl;
    
    cout << "y = [";
    for(int i = 0; i < n+1; i++){
        cout << y[i] << ",";
    }
    cout << "];" << endl;
    
    cout << "xControlPt = [";
    for(int i = 0; i < n; i++){
        cout << xControlPt[i] << ",";
    }
    cout << "];" << endl;
    
    cout << "yControlPt = [";
    for(int i = 0; i < n; i++){
        cout << yControlPt[i] << ",";
    }
    cout << "];" << endl;
    
    cout << "theta = [";
    for(int i = 0; i < n; i++){
        cout << theta[i] << ",";
    }
    cout << "];" << endl;
    
    /*    cout << "beta = [";
     for(int i = 0; i < n; i++){
     for (int j = 0; j < n; j++){
     cout << beta(i,j) << " ";
     }
     if(i < n-1){
     cout << ";";
     
     }
     
     }
     cout << "];"<< endl;*/
    
    cout << "A = [";
    for(int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A(i,j) << " ";
        }
        if(i < n){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    cout << "b = [";
    for(int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << b(i,j) << " ";
        }
        if(i < n){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    
    cout << "c = [";
    for(int i = 0; i < n+1; i++){
        for (int j = 0; j < n+1; j++){
            cout << c(i,j) << " ";
        }
        if(i < n){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    cout << "dl = [";
    for(int i = 0; i < n; i++){
        cout << dl[i];
        if(i < n){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    cout << "Qtan = [";
    for(int i = 0; i < n; i++){
        cout << Qtan[i];
        if(i < n-1){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    cout << "QinfTan = [";
    for(int i = 0; i < n; i++){
        cout << QinfTan[i];
        if(i < n-1){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    cout << "delCl = [";
    for(int i = 0; i < n; i++){
        cout << delCl[i];
        if(i < n-1){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    /*  cout << "radius = [";
     for(int i = 0; i < n; i++){
     for (int j = 0; j < n+1; j++){
     cout << radius(i,j) << " ";
     }
     if(i < n){
     cout << ";";
     }
     }
     cout << "];" << endl;*/
    
    cout << "mu = [";
    for(int i = 0; i < n+1; i++){
        cout << mu[i];
        if(i < n){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    cout << "RHS = [";
    for(int i = 0; i < n; i++){
        cout << RHS[i];
        if(i < n-1){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    cout << "Cp = [";
    for(int i = 0; i < n; i++){
        cout << Cp[i];
        if(i < n-1){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    cout << "sigma = [";
    for(int i = 0; i < n; i++){
        cout << sigma[i];
        if(i < n-1){
            cout << ",";
        }
    }
    cout << "];" << endl;
    
    /*  cout << "figure" << endl;
     cout << "plot(x,y)" << endl;
     cout << "axis([-.2 1 -1 1.5]);" << endl;
     cout << "figure" << endl;*/
    cout << "plot(xControlPt,Cp);" << endl;
    cout << "set(gca,'ydir','reverse')" << endl;
    cout << "xlabel('X');" << endl;
    cout << "ylabel('Cp')" << endl;
    cout << "axis([-.2 1.2 -2 1])" << endl;

    cout << "% ";
}














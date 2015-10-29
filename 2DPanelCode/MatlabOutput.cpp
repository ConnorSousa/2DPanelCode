//
//  MatlabOutput.cpp
//  2DPanelCode
//
//  Created by Connor Sousa on 10/19/15.
//  Copyright (c) 2015 Connor Sousa. All rights reserved.
//

#include "MatlabOutput.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void MatlabOutput::print(MatrixXd &coords, const double &alpha, MatrixXd &controlPts, VectorXd &theta, int &n, MatrixXd &A, MatrixXd &b, MatrixXd &c, VectorXd &dl, VectorXd &Qtan, VectorXd &QinfTan, VectorXd &delCl,  VectorXd &mu, VectorXd &RHS,VectorXd &Cp, VectorXd &sigma){
    
    
    cout << "clear all; close all; format compact;" << endl;
    cout << "alpha = " << alpha*57.2958 << ";"<< endl;
    
    
    cout << "coords = [";
    for(int i = 0; i < coords.rows(); i++){
        for (int j = 0; j < coords.cols(); j++){
            cout << coords(i,j) << " ";
        }
        if(i < coords.rows()){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    cout << "controlPts = [";
    for(int i = 0; i < controlPts.rows(); i++){
        for (int j = 0; j < controlPts.cols(); j++){
            cout << controlPts(i,j) << " ";
        }
        if(i < controlPts.rows()){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    
    cout << "theta = [";
    for(int i = 0; i < theta.size(); i++){
        cout << theta[i] << ",";
    }
    cout << "];" << endl;
    
    cout << "A = [";
    for(int i = 0; i < A.rows(); i++){
        for (int j = 0; j < A.cols(); j++){
            cout << A(i,j) << " ";
        }
        if(i < A.rows()){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    cout << "b = [";
    for(int i = 0; i < b.rows(); i++){
        for (int j = 0; j < b.cols(); j++){
            cout << b(i,j) << " ";
        }
        if(i < b.rows()){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    
    cout << "c = [";
    for(int i = 0; i < c.rows(); i++){
        for (int j = 0; j < c.cols(); j++){
            cout << c(i,j) << " ";
        }
        if(i < c.rows()){
            cout << ";";
        }
    }
    cout << "];" << endl;
    
    cout << "dl = [";
    for(int i = 0; i < dl.size(); i++){
        cout << dl[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "Qtan = [";
    for(int i = 0; i < Qtan.size(); i++){
        cout << Qtan[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "QinfTan = [";
    for(int i = 0; i < QinfTan.size(); i++){
        cout << QinfTan[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "delCl = [";
    for(int i = 0; i < delCl.size(); i++){
        cout << delCl[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "mu = [";
    for(int i = 0; i < mu.size(); i++){
        cout << mu[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "RHS = [";
    for(int i = 0; i < RHS.size(); i++){
        cout << RHS[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "Cp = [";
    for(int i = 0; i < Cp.size(); i++){
        cout << Cp[i] << ",";
    }
    cout << "];" << endl;
    
    
    cout << "sigma = [";
    for(int i = 0; i < sigma.size(); i++){
        cout << sigma[i] << ",";
    }
    cout << "];" << endl;
    /*
     cout << "figure" << endl;
     cout << "plot(controlPts(:,1),Cp);" << endl;
     cout << "set(gca,'ydir','reverse')" << endl;
     cout << "xlabel('X');" << endl;
     cout << "ylabel('Cp')" << endl;
     cout << "axis([-.2 1.2 -2 1])" << endl;
     */
    
}


void MatlabOutput::GridPrint(MatrixXd &gridCoords, VectorXd &uvec, VectorXd &wvec){
    
    cout << "VelGrid = [";
    for(int i = 0; i< gridCoords.rows(); i++){
        cout << gridCoords(i,0) << "," << gridCoords(i,1) << "," << uvec(i) << "," << wvec(i);
        if(i < gridCoords.rows()){
            cout << ";";
        }
    }
    cout << "];" <<endl;
    cout << "figure; hold on;" << " plot(coords(:,1),coords(:,2),'*'); quiver(VelGrid(:,1),VelGrid(:,2),VelGrid(:,3),VelGrid(:,4)); axis equal;" << endl;

    
    
}










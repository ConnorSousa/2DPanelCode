
#include "Geom.h"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;


void Geom::getCoords(MatrixXd &coords, string &file_name){
    
    
    ifstream thefile;
    thefile.open(file_name);
    
    int number_of_lines = 0;
    string line;
    ifstream myfile(file_name);
    
    
    while (getline(myfile, line))
        ++number_of_lines;
    
    coords.conservativeResize(number_of_lines,2);
    
    double x;
    double y;
    int count = 0;
    while(thefile >> x >> y){
        coords(count,0) = x;
        coords(count,1) = y;
        count++;
    }
    thefile.close();
    
    //  n = number_of_lines-1;
}



void Geom::findControlPts(MatrixXd &coords, MatrixXd &controlPts, VectorXd &dl){
    

    for(int i = 0; i < controlPts.rows(); i++){
        
        controlPts(i,0) = (coords(i+1,0) + coords(i,0))/2;
        controlPts(i,1) = (coords(i+1,1) + coords(i,1))/2;
    }
    
    for(int i = 0;i<controlPts.rows()-1;i++){
        dl(i) = sqrt(pow(controlPts(i+1,0)-controlPts(i,0),2) + pow(controlPts(i+1,1)-controlPts(i,1),2));
    }
    dl(controlPts.rows()-1) = dl(controlPts.rows()-2);
}

void Geom::calcPanelAngles(MatrixXd &coords, VectorXd &theta){
    
    for(int i = 0; i < coords.rows()-1; i++){
        double dx = coords(i+1,0) - coords(i,0);
        double dy = coords(i+1,1) - coords(i,1);
        theta[i] = atan2(dy,dx);
    }
}

void Geom::findwakepans(VectorXd &wakePan1, VectorXd &wakePan2, double &c_w, const double &Qinf, MatrixXd &coords, double &dt, double &wakeTheta){

    // Find the angle of last and first panel. Find the average and subtract from first panel to find the angle of bisection. Propagate wake to find panel length and coords.
    
    //"It is usually sufficient to assume that the wake leaves the trailing edge at a median angle delta/2 as shown in fig 13.2" Katz- pg. 376 Buffer Wake

    double delta1 = atan((coords(coords.rows()-2,1)-coords(0,1))/(coords(coords.rows()-2,0)-coords(0,0)));
    double delta2 = atan((coords(0,1)-coords(1,1))/(coords(1,0)-coords(0,0)));
    wakeTheta = delta1-(delta1+delta2)/2;
    
    wakePan1(0) = coords(0,0); // Probably going to need some tricky geom logic when dealing with a harder case.
    wakePan1(1) = coords(0,1);
    wakePan1(2) = wakePan1(0) + c_w*dt*Qinf*cos(wakeTheta);//wakePan1 x2
    wakePan1(3) = wakePan1(1) + c_w*dt*Qinf*sin(wakeTheta);//wakePan1 y2
    
    wakePan2(0) = wakePan1(2);
    wakePan2(1) = wakePan1(3);
    wakePan2(2) = wakePan2(0) + Qinf*dt*cos(wakeTheta);
    wakePan2(3) = wakePan2(1) + Qinf*dt*sin(wakeTheta);
    
}


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

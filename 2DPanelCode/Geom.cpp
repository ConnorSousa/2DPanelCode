
#include "Geom.h"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void Geom::findControlPts(VectorXd &x, VectorXd &y, VectorXd &xControlPt, VectorXd &yControlPt, VectorXd &dl, int &n){
    
    
    for(int i = 0; i < n; i++){
        
        xControlPt[i] = (x[i+1] + x[i])/2;
        yControlPt[i] = (y[i+1] + y[i])/2;
    }
    
    for(int i = 0;i<n-1;i++){
        dl(i) = sqrt(pow(xControlPt(i+1)-xControlPt(i),2) + pow(yControlPt(i+1)-yControlPt(i),2));
    }
    dl(n-1) = dl(n-2); //last panel same as
}

void Geom::calcPanelAngles(VectorXd &x, VectorXd &y, VectorXd &theta, int &n){
    
    for(int i = 0; i < n; i++){
        double dx = x[i+1] - x[i];
        double dy = y[i+1] - y[i];
        theta[i] = atan2(dy,dx);
    }
}

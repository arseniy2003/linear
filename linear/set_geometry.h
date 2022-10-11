// ТУТ ФУНКЦИИ ДЛЯ СОЗДАНИЯ РАЗНЫХ ГЕОМЕТРИЙ, ИХ МОЖЕШЬ НЕ СМОТРЕТЬ ПОКА, НО ТУТ ЕСТЬ ФУНКЦИЯ ДЛЯ СТАРТА ИЗ РАНДОМНОЙ ТОЧКИ
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "assistive_func.h"
using namespace std;
#ifndef set_geometry_h
#define set_geometry_h


vector<vector<double>>  start_init(int n,int dim){
    vector<vector<double>> start;
    start.resize(n);
    for (int ii = 0; ii < n; ii++) {
        vector<double> dot;
        dot.resize(dim);
        for (int i = 0; i < dim; i++)
        {
            dot.at(i) = rand()*((rand()%2)*2-1)/10000000;
        }
        start.at(ii) = dot;
    }
    return start;
}

vector<vector<double>> rand_points_init(int n,int dim){
    vector <vector<double>> points;
    points.resize(n);
    for (int i = 0; i < n; i++)
    {
        vector<double> point;
        point.resize(dim);
        for (int j = 0; j < dim; j++)
        {
            double x;
            x = rand()*((rand()%2)*2-1)/100000;
            point.at(j) = x;
        }
        points.at(i) = point;
    }
    return points;
}

vector<vector<double>> circle_area_points_init(int n,int dim,double r, vector<double> center){
    vector <vector<double>> points;
    points.resize(n);
    int k=0;
    while(k<n)
    {
        
        vector<double> point;
        point.resize(dim);
        double rr;
        rr = ((rand()%((int)r)));
        if(dim==2){
            double angle = rand();
            point.at(0) = rr*sin(angle)+center.at(0);
            point.at(1) = rr*cos(angle)+center.at(1);
            
        }
        if(dim==3){
            double alpha=rand();
            double beta=rand();
            point.at(0)=rr*cos(alpha)+center.at(0);
            point.at(1)=rr*cos(beta)+center.at(1);
            point.at(2)=sqrt(abs(rr*rr-pow(point.at(0),2)-pow(point.at(0),2)))+center.at(2);
            
            
        }
        points.at(k) = point;
        bool a =false;
        for (int j=0; j<k; j++) {
            if(ro(points.at(j),points.at(k))<r/2){
                a=true;
            }
        }
        if(a==true){
            k--;
        }
        k++;
    }
    return points;
}


vector<vector<double>> arc_points_init(int n,int dim, double r, vector<double> center, double a,vector<double> vv){
    vector <vector<double>> points;
    points.resize(n);
    double da=a/(n-1);
    double a0=atan(vv.at(0)/vv.at(1))-a/2;
    for (int i = 0; i < n; i++)
    {
        points.at(i)=center;
    }
    for (int i = 0; i < n; i++)
    {
        
        points.at(i).at(0)+=r*cos(a0);
        points.at(i).at(1)+=r*sin(a0);
        a0+=da;
    }
    return points;
}


#endif /* set_geometry_h */

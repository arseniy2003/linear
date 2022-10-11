// ТУТ ФУНКЦИИ, СЧИТАЮЩИЕ ГРАДИЕНТ
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "assistive_func.h"
#include "error_functions.h"
#include "recalculate_from_broken.h"
using namespace std;

#ifndef gradients_h
#define gradients_h

vector<double> space_gradcalc( vector<double> point, vector<vector<double> > points,  vector<double> deltas_g, vector<double>delays) {
    double eps=0.01;
    vector<double> grad;
    int dim = point.size();
    grad.resize(dim);
    double err0= ro(full_syst_improve({get_space_deltas(dim, points, point)}, delays).at(0), deltas_g);
    for (int i = 0; i < dim; i++)
    {
        vector<double> pointx = point;
        pointx.at(i) = point.at(i) + eps;
        double errcur= ro( full_syst_improve({get_space_deltas(dim, points, pointx)}, delays).at(0), deltas_g);
        double ddx = (errcur - err0) / eps;
        grad.at(i) = ddx;
    }
    return grad;
}


vector<vector<double> > usual_gradcalc1(int dim, vector<vector<double> > points,vector<vector<double> > spots,vector<vector<double> > deltas_g){
    double eps=0.001;
    vector<vector<double> > grad;
    grad.resize(points.size());
    double err0 = usual_errcalc(set_full_deltas(dim, points, spots), deltas_g);
    for (int i=0; i<grad.size(); i++) {
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j =0; j<dim; j++) {
            vector<vector<double> > temppoints = points;
            temppoints.at(i).at(j)+=eps;
            double errcur = usual_errcalc(set_full_deltas(dim, temppoints, spots), deltas_g);
            tempgrad.at(j)=(errcur-err0)/eps;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > usual_gradcalc2(int dim, vector<vector<double> > points,vector<vector<double> > spots,vector<vector<double> > deltas_g){
    double eps=0.01;
    vector<vector<double> > grad;
    grad.resize(points.size());
    double err0 = usual_errcalc(set_full_deltas(dim, spots, points), deltas_g);
    for (int i=0; i<grad.size(); i++) {
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j =0; j<dim; j++) {
            vector<vector<double> > temppoints = points;
            temppoints.at(i).at(j)+=eps;
            double errcur = usual_errcalc(set_full_deltas(dim, spots, temppoints), deltas_g);
            tempgrad.at(j)=(errcur-err0)/eps;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > usual_gradcalc3(int dim, vector<vector<double> > points,vector<vector<double> > spots,vector<vector<double> > deltas_g1,vector<vector<double> > deltas_g2){
    double eps=0.01;
    vector<vector<double> > grad;
    grad.resize(points.size());
    double err0 = usual_errcalc(set_full_deltas(dim, spots, points), deltas_g2)+usual_errcalc(set_full_deltas(dim, points, spots), deltas_g1);
    for (int i=0; i<grad.size(); i++) {
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j =0; j<dim; j++) {
            vector<vector<double> > temppoints = points;
            temppoints.at(i).at(j)+=eps;
            double errcur =usual_errcalc(set_full_deltas(dim, temppoints, spots), deltas_g1) + usual_errcalc(set_full_deltas(dim, spots, temppoints), deltas_g2);
            tempgrad.at(j)=(errcur-err0)/eps;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > halftime_gradcalc1(int dim, vector<vector<double> > points,vector<vector<double> > spots, vector<vector<double> > ranges, vector<vector<double> > deltas_g1,vector<vector<double> > deltas_g2){
    double eps=0.01;
    vector<vector<double> > grad;
    grad.resize(points.size());
    double err0 = halftime_errcalc1(dim, points, spots, ranges, deltas_g1, deltas_g2);
    for (int i=0; i<grad.size(); i++) {
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j =0; j<dim; j++) {
            vector<vector<double> > temppoints = points;
            temppoints.at(i).at(j)+=eps;
            double errcur = halftime_errcalc1(dim, temppoints, spots, ranges, deltas_g1, deltas_g2);
            tempgrad.at(j)=(errcur-err0)/eps;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > halftime_gradcalc2(int dim, vector<vector<double> > points,vector<vector<double> > spots, vector<vector<double> > ranges, vector<vector<double> > deltas_g1,vector<vector<double> > deltas_g2, vector<int> order){
    double eps=0.01;
    vector<vector<double> > grad;
    grad.resize(points.size());
    double err0 = halftime_errcalc2(dim, points, spots, ranges, deltas_g1, deltas_g2,order);
    for (int i=0; i<grad.size(); i++) {
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j =0; j<dim; j++) {
            vector<vector<double> > temppoints = points;
            temppoints.at(i).at(j)+=eps;
            double errcur = halftime_errcalc2(dim, temppoints, spots, ranges, deltas_g1, deltas_g2,order);
            tempgrad.at(j)=(errcur-err0)/eps;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > abstract_gradcalc(int dim, vector<vector<double> > points, vector<vector<double> > spots, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2, vector<vector<int>> orders){
    double eps=0.01;
    int sz=points.size();
    vector<vector<double> > grad;
    vector<vector<double> > temp_points;
    
    grad.resize(sz);
    double err0=0;
    for(int kk=0;kk<orders.size();kk++){
        err0 += abstract_errcalc(dim,points,spots,broken_deltas1,broken_deltas2,orders.at(kk));
    }
    for (int i=0; i<sz; i++) {
        temp_points=points;
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j = 0; j < dim; j++)
        {
            temp_points.at(i).at(j) += eps;
            double errcur=0;
            for(int kk=0;kk<orders.size();kk++){
                errcur += abstract_errcalc(dim,temp_points,spots,broken_deltas1,broken_deltas2,orders.at(kk));
            }
            double ddx = (errcur - err0) / eps;
            tempgrad.at(j) = ddx;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > abstract_gradcalc_two_halfs(int dim, vector<vector<double> > points, vector<vector<double> > spots, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2, vector<int> order1,vector<int> order2){
    double eps=0.01;
    int sz=points.size();
    vector<vector<double> > grad;
    vector<vector<double> > temp_points;
    
    grad.resize(sz);
    double err0=0;
    
    err0 += abstract_errcalc_half(dim, points, spots, broken_deltas1, order1);
    err0 += abstract_errcalc_half(dim, spots, points, broken_deltas2, order2);
    for (int i=0; i<sz; i++) {
        temp_points=points;
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j = 0; j < dim; j++)
        {
            temp_points.at(i).at(j) += eps;
            double errcur=0;
            errcur +=  abstract_errcalc_half(dim, temp_points, spots, broken_deltas1, order1);
            errcur +=  abstract_errcalc_half(dim, spots, temp_points, broken_deltas2, order2);
            double ddx = (errcur - err0) / eps;
            tempgrad.at(j) = ddx;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

vector<vector<double> > abstract_gradcalc_one_half1(int dim, vector<vector<double> > points, vector<vector<double> > spots, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2, vector<int> order1,vector<int> order2){
    double eps=0.01;
    int sz=points.size();
    vector<vector<double> > grad;
    vector<vector<double> > temp_points;
    
    grad.resize(sz);
    double err0=0;
    
    err0 += 10*usual_errcalc(set_full_deltas(dim, points, spots), broken_deltas1);
    err0 += abstract_errcalc_half(dim, spots, points, broken_deltas2, order2);
    for (int i=0; i<sz; i++) {
        temp_points=points;
        vector<double> tempgrad;
        tempgrad.resize(dim);
        for (int j = 0; j < dim; j++)
        {
            temp_points.at(i).at(j) += eps;
            double errcur=0;
            errcur +=  usual_errcalc(set_full_deltas(dim, temp_points, spots), broken_deltas1);
            errcur +=  abstract_errcalc_half(dim, spots, temp_points, broken_deltas2, order2);
            double ddx = (errcur - err0) / eps;
            tempgrad.at(j) = ddx;
        }
        grad.at(i)=tempgrad;
    }
    return grad;
}

#endif /* gradients_h */

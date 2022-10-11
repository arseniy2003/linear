// ТУТ ФУНКЦИИ ОШИБОК
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "assistive_func.h"
#include "recalculate_from_broken.h"
using namespace std;

#ifndef error_functions_h
#define error_functions_h

double usual_errcalc(vector<vector<double> > found_deltas,vector<vector<double> > deltas_g){
    double err=0;
    for (int i = 0; i<deltas_g.size(); i++) {
        err+=pow(ro(found_deltas.at(i),deltas_g.at(i)),2);
    }
    return sqrt(err);
}

double abstract_errcalc(int dim,vector<vector<double> >points_c, vector<vector<double> > spots,vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2,vector<int> order){
    double err=0;
    vector<vector<double> > syst1=abstract_syst_improve(broken_deltas1, order);
    vector<vector<double> > syst2=abstract_syst_improve(broken_deltas2, order);
    vector<vector<double> > deltas1=set_full_deltas(dim, points_c, spots);
    vector<vector<double> > deltas2=set_full_deltas(dim, spots, points_c);
    vector<vector<double> > improved_syst1=abstract_syst_improve(deltas1,order);
    vector<vector<double> > improved_syst2=abstract_syst_improve(deltas2,order);
    for (int i=0; i<syst1.size(); i++) {
        err+=pow(ro(improved_syst1.at(i),syst1.at(i)),2)+pow(ro(improved_syst2.at(i),syst2.at(i)),2);
    }
    
    return sqrt(err);
}

double abstract_errcalc_half(int dim,vector<vector<double> >points_c, vector<vector<double> > spots,vector<vector<double> > broken_deltas1,vector<int> order){
    double err=0;
    vector<vector<double> > syst1=abstract_syst_improve(broken_deltas1, order);
    vector<vector<double> > deltas1=set_deltas(dim, points_c, spots);
    vector<vector<double> > improved_syst1=abstract_syst_improve(deltas1,order);
    for (int i=0; i<syst1.size(); i++) {
        err+=pow(ro(improved_syst1.at(i),syst1.at(i)),2);
    }
    
    return sqrt(err);
}

double halftime_errcalc1(int dim,vector<vector<double> >points_c, vector<vector<double> > spots,vector<vector<double> > ranges,vector<vector<double> > deltas1,vector<vector<double> > deltas2){
    double err=0;
    vector<vector<double> > found_deltas1=set_full_deltas(dim, points_c, spots);
    vector<vector<double> > found_deltas2=set_full_deltas(dim, spots, points_c);
    for (int i = 0; i < points_c.size(); i++) {
        for(int j = 0; j < points_c.size(); j++){
            err+=10*pow(ro(points_c.at(i),points_c.at(j))-ranges.at(i).at(j),2);
        }
    }
    for (int i = 0; i<deltas1.size(); i++) {
        err+=pow(ro(found_deltas1.at(i),deltas1.at(i)),2);
    }
    for (int i = 0; i<deltas2.size(); i++) {
        err+=pow(ro(found_deltas2.at(i),deltas2.at(i)),2);
    }
    return sqrt(err);
}

double halftime_errcalc2(int dim,vector<vector<double> >points_c, vector<vector<double> > spots,vector<vector<double> > ranges,vector<vector<double> > deltas1,vector<vector<double> > deltas2, vector<int> order){
    double err=0;
    vector<vector<double> > found_deltas1=set_full_deltas(dim, points_c, spots);
    vector<vector<double> > found_deltas2=set_full_deltas(dim, spots, points_c);
    for (int i = 0; i < points_c.size(); i++) {
        for(int j = 0; j < points_c.size(); j++){
            err+=10*pow(ro(points_c.at(i),points_c.at(j))-ranges.at(i).at(j),2);
        }
    }
    for (int i = 0; i<deltas2.size(); i++) {
        err+=pow(ro(found_deltas2.at(i),deltas2.at(i)),2);
    }
    err+=abstract_errcalc_half(dim, points_c, spots, found_deltas1, order);
    return sqrt(err);
}


#endif /* error_functions_h */

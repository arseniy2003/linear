// ТУТ ФУНКЦИИ, КОТОРЫЕ МОДЕЛИРУЮТ ИЗМЕРЕНИЯ И ПОРТЯТ ИХ
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "assistive_func.h"
using namespace std;

#ifndef simulating_dim_h
#define simulating_dim_h


vector<double> get_space_deltas( int dim, vector<vector<double> > points,vector<double> point) {
    vector<double> deltas_c;
    int m=points.size();
    deltas_c.resize(m-1);
    for (int i = 0; i < m-1; i++)
    {
        deltas_c.at(i) = ro(point, points.at(0)) - ro(point, points.at(i+1));
    }
    return deltas_c;
}

vector<double> set_delays(int n){
    vector<double> delays;
    delays.resize(n);
    for (int i=0; i<n; i++) {
        delays.at(i)=((double)(rand()%100000));
    }
    return delays;
}

vector<double> set_zero_delays(int n){
    vector<double> delays;
    delays.resize(n);
    for (int i=0; i<n; i++) {
        delays.at(i)=0;
    }
    return delays;
}

vector<double> set_small_delays(int n, double SKO_delays){
    vector<double> delays;
    delays.resize(n);
    for (int i=0; i<n; i++) {
        delays.at(i)=gaussrand(0, SKO_delays);
    }
    return delays;
}

vector<double> make_diff_delays(vector<double> delays){
    int n=delays.size();
    vector<double> diff_delays;
    diff_delays.resize(n-1);
    for (int i=0; i<n-1; i++) {
        diff_delays.at(i)=delays.at(0)-delays.at(i+1);
    }
    return diff_delays;
}

vector<double> make_full_diff_delays(vector<double> delays){
    int n=delays.size();
    vector<double> diff_delays;
    diff_delays.resize(0);
    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            diff_delays.push_back(delays.at(i)-delays.at(j));
        }
    }
    return diff_delays;
}

vector<vector<double>> set_length (vector<vector<double> > points, vector<vector<double> > spots, double SKO, vector<double> delays1,vector<double> delays2){
    vector<vector<double>> length;
    length.resize(points.size());
    for(int i=0;i<points.size();i++){
        vector<double> sublength;
        sublength.resize(spots.size());
        for(int j=0; j<spots.size(); j++){
            sublength.at(j)= gaussrand(ro(points.at(i),spots.at(j)), SKO) + delays1.at(i) + delays2.at(j);
        }
        length.at(i)=sublength;
    }
    return length;
}

vector<vector<double> > set_deltas_correct(vector<vector<double> > points, vector<vector<double> > spots, double SKO, vector<double> delays1,vector<double> delays2){
    vector<vector<double> > length=set_length(points,spots,SKO,delays1,delays2);
    vector<vector<double> > deltas;
    deltas.resize(length.size());
    for(int i = 0; i <length.size();i++){
        vector<double> delta;
        delta.resize(length.at(i).size()-1);
        for(int j = 0; j <length.at(i).size()-1;j++){
            delta.at(j)=length.at(i).at(0)-length.at(i).at(j+1);
        }
        deltas.at(i)=delta;
    }
    return deltas;
}

vector<vector<double> > set_full_deltas_correct(vector<vector<double> > points, vector<vector<double> > spots, double SKO, vector<double> delays1,vector<double> delays2){
    vector<vector<double> > length=set_length(points,spots,SKO,delays1,delays2);
    vector<vector<double> > deltas;
    deltas.resize(length.size());
    for(int i = 0; i <length.size();i++){
        vector<double> delta;
        delta.resize(0);
        for(int j = 0; j <length.at(i).size()-1;j++){
            for (int k = j+1; k<length.at(i).size(); k++) {
                delta.push_back(length.at(i).at(j)-length.at(i).at(k));
            }
            
        }
        deltas.at(i)=delta;
    }
    return deltas;
}

vector<vector<double> > set_deltas( double dim, vector<vector<double> > points, vector<vector<double> > spots) {
    int n1=points.size();
    int n2=spots.size();
    vector<vector<double> > deltas_c;
    deltas_c.resize((n1));
    for (int i = 0; i < n1; i++)
    {
        vector<double> syst;
        syst.resize(n2-1);
        for (int j=1;j<n2;j++){
            syst.at(j-1)=ro(points.at(i),spots.at(0))-ro(points.at(i),spots.at(j));
        }
        deltas_c.at(i)=syst;
    }
    return deltas_c;
}

vector<vector<double> > set_full_deltas( double dim, vector<vector<double> > points, vector<vector<double> > spots) {
    int n1=points.size();
    int n2=spots.size();
    vector<vector<double> > deltas_c;
    deltas_c.resize((n1));
    for (int i = 0; i < n1; i++)
    {
        vector<double> syst;
        syst.resize(0);
        for (int k = 0; k < n2-1; k++) {
            for (int j = k+1; j < n2; j++){
                syst.push_back(ro(points.at(i),spots.at(k))-ro(points.at(i),spots.at(j)));
            }
        }
        deltas_c.at(i)=syst;
    }
    return deltas_c;
}

vector<vector<double>> set_ranges(int dim,vector<vector<double>> points){
    vector<vector<double>> ranges;
    ranges.resize(points.size());
    for (int i = 0; i < points.size(); i++) {
        vector<double> temp_ranges;
        temp_ranges.resize(points.size());
        for(int j = 0; j < points.size(); j++){
            temp_ranges.at(j)=ro(points.at(i),points.at(j));
        }
        ranges.at(i)=temp_ranges;
    }
    return ranges;
}


vector<vector<double> > break_deltas(vector<vector<double> > deltas,vector<double> delays,double SKO){
    vector<vector<double> > broken_deltas=deltas;
    vector<double> diff_delays=make_full_diff_delays(delays);
    for(int i = 0;i<deltas.size();i++){
        for(int j =0;j<deltas.at(i).size();j++){
            broken_deltas.at(i).at(j)=gaussrand(broken_deltas.at(i).at(j), SKO)+diff_delays.at(j);
        }
    }
    return broken_deltas;
}

vector<vector<double> > break_ranges(vector<vector<double> > ranges,double SKO){
    vector<vector<double> > broken_ranges = ranges;
    for(int i = 0; i < ranges.size(); i++){
        for (int j = 0; j < ranges.at(i).size(); j++) {
            broken_ranges.at(i).at(j)=gaussrand(broken_ranges.at(i).at(j), SKO);
        }
    }
    return(broken_ranges);
}

vector<vector<double>> break_ranges_with_delays(vector<vector<double>> ranges1, double SKO, vector<double> delays){
    vector<vector<double>> ranges = ranges1;
    for (int i = 0; i < ranges.size(); i++) {
        for (int j = 0; j < ranges.at(i).size(); j++) {
            ranges.at(i).at(j)=gaussrand(ranges.at(i).at(j), SKO)+delays.at(i)-delays.at(j);
        }
    }
    return ranges;
}


#endif /* simulating_dim_h */

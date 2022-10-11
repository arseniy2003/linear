//ТУТ ФУНКЦИИ, КОТОРЫЕ СОВЕРШАЮТ ОПЕРАЦИИ С УРАВНЕНИЯМИ: ВЫЧИТАЮТ  ИХ ДЖРУГ ИЗ ДРУГА, УБИРАЮТ ЗАДЕРЖКИ И ПРОЧЕЕ
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "assistive_func.h"
using namespace std;
#ifndef recalculate_from_broken_h
#define recalculate_from_broken_h

vector<vector<double> > syst_improve(vector<vector<double> > deltas){
    vector<vector<double> > syst;
    int n1=deltas.size();
    syst.resize(n1-1);
    for (int i=1; i<n1; i++) {
        syst.at(i-1)=deltas.at(i);
        for (int j =0; j<deltas.at(i).size(); j++) {
            syst.at(i-1).at(j)=syst.at(i-1).at(j)-deltas.at(0).at(j);
        }
    }
    return syst;
}

vector<vector<double> > better_syst_improve(vector<vector<double> > deltas){
    vector<vector<double> > syst;
    syst=deltas;
    double sz1=deltas.size();
    double sz2=deltas.at(0).size();
    for (int j =0; j<sz2; j++) {
        for (int i=0; i<sz1; i++) {
            syst.at(i).at(j)=syst.at(i).at(j)-deltas.at(j).at(j);
        }
    }
    return syst;
}

vector<vector<double> > full_syst_improve(vector<vector<double> > deltas,vector<double> delays){
    vector<vector<double> > syst;
    syst=deltas;
    double sz1=deltas.size();
    double sz2=delays.size();
    for (int j =0; j<sz2; j++) {
        for (int i=0; i<sz1; i++) {
            syst.at(i).at(j)=syst.at(i).at(j)+delays.at(j);
        }
    }
    return syst;
}

vector<vector<double> > full_syst_improve_minus(vector<vector<double> > deltas,vector<double> delays){
    vector<vector<double> > syst;
    syst=deltas;
    double sz1=deltas.size();
    double sz2=delays.size();
    for (int j =0; j<sz2; j++) {
        for (int i=0; i<sz1; i++) {
            syst.at(i).at(j)=syst.at(i).at(j)-delays.at(j);
        }
    }
    return syst;
}

vector<vector<double> > abstract_syst_improve(vector<vector<double> > deltas,vector<int> order){
    vector<vector<double> > syst;
    syst=deltas;
    for (int i=0; i<deltas.size(); i++) {
        for (int j =0; j<deltas.at(i).size(); j++) {
            syst.at(i).at(j)=syst.at(i).at(j)-deltas.at(order.at(j)).at(j);
        }
    }
    return syst;
}

vector<vector<double>> calculate_ranges(vector<vector<double>> broken_ranges){
    vector<vector<double>> ranges = broken_ranges;
    for (int i = 0; i < ranges.size(); i++) {
        for (int j = 0; j < ranges.at(i).size(); j++) {
            ranges.at(i).at(j)=(broken_ranges.at(i).at(j)+broken_ranges.at(j).at(i))/2;
        }
    }
    return ranges;
}

vector<double> calculate_diff_delays(vector<vector<double>> broken_ranges){
    vector<double> diff_delays;
    diff_delays.resize(0);
    for (int i = 0; i < broken_ranges.size()-1; i++) {
        for (int j = i+1; j < broken_ranges.at(i).size(); j++) {
            diff_delays.push_back((broken_ranges.at(i).at(j)-broken_ranges.at(j).at(i))/2);
        }
    }
    return diff_delays;
}

vector<double> recalculate_delays(vector<vector<double>> diff_delays){
    vector<double> delays;
    delays.resize(diff_delays.size()+1);
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return delays;
}

vector<vector<double> > recalculate_ranges(vector<vector<double>> broken_ranges,vector<double> delays){
    vector<vector<double>> ranges = broken_ranges;
    for (int i = 0; i < ranges.size(); i++) {
        for (int j = 0; j < ranges.at(i).size(); j++) {
            ranges.at(i).at(j)=(broken_ranges.at(i).at(j)-delays.at(i)+delays.at(j));
        }
    }
    return ranges;
}

vector<vector<double> > recalculate_deltas(vector<vector<double>> broken_deltas, vector<double>  diff_delays){
    vector<vector<double> > calculated_deltas=broken_deltas;
    for (int i = 0; i < calculated_deltas.size(); i++) {
        for (int  j = 0; j < calculated_deltas.at(i).size(); j++) {
                calculated_deltas.at(i).at(j)=calculated_deltas.at(i).at(j)-diff_delays.at(j);
        }
    }
    return calculated_deltas;
}


#endif /* recalculate_from_broken_h */

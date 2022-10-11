//ТУТ ОПИСАНЫ САМИ ГРАДИЕНТНЫЕ СПУСКИ
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "assistive_func.h"
#include "gradients.h"
#include "error_functions.h"
#include "recalculate_from_broken.h"
using namespace std;

#ifndef gradient_descents_h
#define gradient_descents_h

vector<double> space_search_points(int dim, vector<double>  starts,  vector<vector<double> > points,  vector<double> broken_deltas,vector<double> delays) {
    int times=1;
    double faslamda=0.005;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    vector<double> syst=full_syst_improve_minus({broken_deltas},delays).at(0);
    vector<double>  grad;
    vector<vector<double>> v=zero_vect(1,dim);
    vector<double>  point_c;
    vector<vector<vector<vector<double>>>> founded;
    founded.resize(times);
    for (int k=0; k<times; k++) {
        point_c=starts;
        bool a = true;
        while (a) {
            grad=space_gradcalc(point_c, points, syst,delays);
            // printvec(grad);
            v=v_change(gamma, v, {grad});
            //cout<<endl<<abs_multivect(v);
            cach=cache_change(beta, cach, {grad});
            //cout<<endl<<cach;
                for (int j = 0; j < dim; j++)
                {
                    point_c.at(j) = point_c.at(j) - faslamda * v.at(0).at(j)/sqrt(cach+0.0001);
                }
            
            //printvec(point_c);
            step++;
            if(step>500){if((abs_multivect(v)<0.0005)){a=false;}}
        }
    }
    return point_c;
}


vector<vector<double>> usual_searchpoint(int dim,int n, int l,vector<vector<double> > starts,  vector<vector<double> > spots,  vector<vector<double>> deltas_g1, vector<vector<double>> deltas_g2){
    int times=1;
    double faslamda=0.003;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    bool a = true;
    vector<vector<double>>  grad;
    vector<vector<double>> v=zero_vect(n,dim);
    vector<vector<double>>  points_c = starts;
    while(a){
        if(l==1){
        grad=usual_gradcalc1(dim, points_c, spots, deltas_g1);
        }
        if(l==2){
            grad=usual_gradcalc2(dim, points_c, spots, deltas_g2);
        }
        if(l==3){
            grad=usual_gradcalc3(dim, points_c, spots, deltas_g1, deltas_g2);
        }
        v=v_change(gamma, v, grad);
        //cout<<endl<<abs_multivect(v);
        cach=cache_change(beta, cach, grad);
        for (int jj=0; jj<n; jj++) {
            for (int j = 0; j < dim; j++)
            {
                points_c.at(jj).at(j) = points_c.at(jj).at(j) - faslamda * v.at(jj).at(j)/sqrt(cach+0.0000001);
            }
        }
        step++;
        if(step>500){if((abs_multivect(v)<0.0008)){a=false;}}
    }
    
    return(points_c);
}

vector<vector<double>> halftime_searchpoints1(int dim,int n, vector<vector<double> > starts,  vector<vector<double> > spots, vector<vector<double>> ranges,  vector<vector<double>> deltas_g1, vector<vector<double>> deltas_g2){
    int times=1;
    double faslamda=0.003;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    bool a = true;
    vector<vector<double>>  grad;
    vector<vector<double>> v=zero_vect(n,dim);
    vector<vector<double>>  points_c = starts;
    while(a){
        grad=halftime_gradcalc1(dim, points_c, spots,ranges, deltas_g1, deltas_g2);
        v=v_change(gamma, v, grad);
        //cout<<endl<<abs_multivect(v);
        cach=cache_change(beta, cach, grad);
        for (int jj=0; jj<n; jj++) {
            for (int j = 0; j < dim; j++)
            {
                points_c.at(jj).at(j) = points_c.at(jj).at(j) - faslamda * v.at(jj).at(j)/sqrt(cach+0.0000001);
            }
        }
        step++;
        if(step>500){if((abs_multivect(v)<0.0008)){a=false;}}
    }

return(points_c);
}

vector<vector<double>> halftime_searchpoints2(int dim,int n, vector<vector<double> > starts,  vector<vector<double> > spots, vector<vector<double>> ranges,  vector<vector<double>> deltas_g1, vector<vector<double>> deltas_g2, vector<int> order){
    int times=1;
    double faslamda=0.003;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    bool a = true;
    vector<vector<double>>  grad;
    vector<vector<double>> v=zero_vect(n,dim);
    vector<vector<double>>  points_c = starts;
    while(a){
        grad=halftime_gradcalc2(dim, points_c, spots,ranges, deltas_g1, deltas_g2,order);
        v=v_change(gamma, v, grad);
        //cout<<endl<<abs_multivect(v);
        cach=cache_change(beta, cach, grad);
        for (int jj=0; jj<n; jj++) {
            for (int j = 0; j < dim; j++)
            {
                points_c.at(jj).at(j) = points_c.at(jj).at(j) - faslamda * v.at(jj).at(j)/sqrt(cach+0.0000001);
            }
        }
        step++;
        if(step>500){if((abs_multivect(v)<0.0008)){a=false;}}
    }

return(points_c);
}

vector<vector<double>> abstract_search_points(int dim,int n, vector<vector<double> > starts,  vector<vector<double> > spots,  vector<vector<double>> broken_deltas1,vector<vector<double>> broken_deltas2, vector<vector<int > > orders) {
    int times=1;
    double faslamda=0.005;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    vector<vector<double>> grad;
    vector<vector<double>> v=zero_vect(n,dim);
    vector<vector<double>>  points_c;
    for (int k=0; k<times; k++) {
        points_c=starts;
        bool a = true;
        while (a) {
            grad=abstract_gradcalc(dim, points_c, spots, broken_deltas1, broken_deltas2,orders);
            v=v_change(gamma, v, grad);
            cach=cache_change(beta, cach, grad);
            for (int jj=0; jj<n; jj++) {
                for (int j = 0; j < dim; j++)
                {
                    points_c.at(jj).at(j) = points_c.at(jj).at(j) - faslamda * v.at(jj).at(j)/sqrt(cach+0.0000001);
                }
            }
            step++;
            
            if(step>500){if((abs_multivect(v)<0.001)){a=false;}}
        }
    }
    return points_c;
}

vector<vector<double>> abstract_search_points_two_halfs(int dim,int n, vector<vector<double> > starts,  vector<vector<double> > spots,  vector<vector<double>> broken_deltas1,vector<vector<double>> broken_deltas2, vector<int >  order1, vector<int >  order2) {
    int times=1;
    double faslamda=0.005;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    vector<vector<double>> grad;
    vector<vector<double>> v=zero_vect(n,dim);
    vector<vector<double>>  points_c;
    for (int k=0; k<times; k++) {
        points_c=starts;
        bool a = true;
        while (a) {
            grad=abstract_gradcalc_two_halfs(dim, points_c, spots, broken_deltas1, broken_deltas2, order1, order2);
            v=v_change(gamma, v, grad);
            cach=cache_change(beta, cach, grad);
            for (int jj=0; jj<n; jj++) {
                for (int j = 0; j < dim; j++)
                {
                    points_c.at(jj).at(j) = points_c.at(jj).at(j) - faslamda * v.at(jj).at(j)/sqrt(cach+0.0000001);
                }
            }
            step++;
            
            if(step>500){if((abs_multivect(v)<0.001)){a=false;}}
        }
    }
    return points_c;
}

vector<vector<double>> abstract_search_points_one_half(int dim,int n, vector<vector<double> > starts,  vector<vector<double> > spots,  vector<vector<double>> broken_deltas1,vector<vector<double>> broken_deltas2, vector<int >  order1, vector<int >  order2) {
    int times=1;
    double faslamda=0.005;
    double gamma=0.9;
    double beta=0.999;
    double cach=0;
    double err=0;
    int step = 0;
    vector<vector<double>> grad;
    vector<vector<double>> v=zero_vect(n,dim);
    vector<vector<double>>  points_c;
    for (int k=0; k<times; k++) {
        points_c=starts;
        bool a = true;
        while (a) {
            grad=abstract_gradcalc_one_half1(dim, points_c, spots, broken_deltas1, broken_deltas2, order1, order2);
            v=v_change(gamma, v, grad);
            cout<<endl<<abs_multivect(v)<<endl;
            cach=cache_change(beta, cach, grad);
            for (int jj=0; jj<n; jj++) {
                for (int j = 0; j < dim; j++)
                {
                    points_c.at(jj).at(j) = points_c.at(jj).at(j) - faslamda * v.at(jj).at(j)/sqrt(cach+0.0000001);
                }
            }
            step++;
            
            if(step>500){if((abs_multivect(v)<0.001)){a=false;}}
        }
    }
    return points_c;
}


#endif /* gradient_descents_h */

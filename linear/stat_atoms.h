//ЭТО СОЗДАНО ЧИСТО РАДИ МНОГОПОТОЧКИ, ПОКА ПО СУТИ Я ЗДЕСЬ ПРОСТО ВЫЗЫВАЮ НЕКОТОРЫЕ ГРАДИЕНТНЫЕ СПУСКИ
#include <math.h>
#include <fstream>
#include <thread>
#include <iostream>
#include <vector>
#include "assistive_func.h"
#include "simulating_dim.h"
#include "set_geometry.h"
#include "recalculate_from_broken.h"
#include "error_functions.h"
#include "gradients.h"
#include "gradient_descents.h"
using namespace std;

#ifndef stat_atoms_h
#define stat_atoms_h


vector<vector<double>> stat_atom1( int dim, vector<vector<double>> points, vector<vector<double>> spots,vector<vector<double> > broken_ranges, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2){
    int n1=points.size();
    vector<vector<double> > calculated_ranges = calculate_ranges(broken_ranges);
    vector<double> calculated_diff_delays = calculate_diff_delays(broken_ranges); // тут их гораздо больше и есть нули
    vector<vector<double> > calculated_deltas2=recalculate_deltas(broken_deltas2, calculated_diff_delays);
    
    return halftime_searchpoints1(dim, n1, points, spots, calculated_ranges, broken_deltas1, calculated_deltas2);
}

vector<vector<double>> stat_atom2( int dim, vector<vector<double>> points, vector<vector<double>> spots,vector<vector<double> > broken_ranges, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2){
    int n1=points.size();
    int n2=spots.size();
    vector<int> order1 = set_order1(n2, n1);
    vector<vector<double> > calculated_ranges = calculate_ranges(broken_ranges);
    vector<double> calculated_diff_delays = calculate_diff_delays(broken_ranges); // тут их гораздо больше и есть нули
    vector<vector<double> > calculated_deltas2=recalculate_deltas(broken_deltas2, calculated_diff_delays);
    
    return halftime_searchpoints2(dim, n1, points, spots, calculated_ranges, broken_deltas1, calculated_deltas2,order1);
}

vector<vector<double>> stat_atom3( int dim, vector<vector<double>> points, vector<vector<double>> spots,vector<vector<double> > broken_ranges, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2){
    int n1=points.size();
    int n2=spots.size();
    vector<int> order1 = set_order1(n2, n1);
    vector<int> order2 = set_order1(n1, n2);
    return abstract_search_points_two_halfs(dim, n1, points, spots, broken_deltas1, broken_deltas2,order1, order2);
}

vector<vector<double>> stat_atom4( int dim, vector<vector<double>> points, vector<vector<double>> spots,vector<vector<double> > broken_ranges, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2){
    int n1=points.size();
    int n2=spots.size();
    vector<int> order1 = set_order1(n2, n1);
    vector<int> order2 = set_order1(n1, n2);
    return abstract_search_points_one_half(dim, n1, points, spots, broken_deltas1, broken_deltas2,order1, order2);
}

vector<vector<double>> stat_atom121(double SKO, int dim, vector<vector<double>> points, vector<vector<double>> spots){
    int n1=points.size();
    int n2=spots.size();
    vector <vector<int>> orders1;
    orders1.resize(0);
    vector <vector<int>> orders2;
    orders2.resize(0);
    for(int ii=0;ii<5;ii++){
        for(int iii=0;iii<5;iii++){
            for(int iiii=0;iiii<5;iiii++){
                for(int iiiii=0;iiiii<5;iiiii++){
                    orders1.push_back({ii,iii,iiii,iiiii});
                    if((ii!=iii)&(ii!=iiii)&(ii!=iiiii)&(iii!=iiii)&(iii!=iiiii)&(iiii!=iiiii)){
                        orders2.push_back({ii,iii,iiii,iiiii});
                    }
                }
            }
        }
    }
    
    
    cout<<endl<<orders2.size();
    vector<double> delays1=set_delays(n1);
    vector<double> delays2=set_delays(n2);
    vector<double> diff_delays1=make_diff_delays(delays1);
    vector<double> diff_delays2=make_diff_delays(delays2);
    vector<vector<vector<double>>> stat;
    stat.resize(100);
    
    int threads_num=4;
    for (int j = 0; j<1; j++) {
        vector<vector<double> > broken_deltas1=set_deltas_correct(points, spots, SKO, delays1, delays2);
        vector<vector<double> > broken_deltas2=set_deltas_correct(spots, points, SKO, delays2, delays1);
        vector<vector<vector<double>>> temp_stat;
        temp_stat.resize(0);

        vector <vector<int>> orders={{orders2.at(0)},{orders2.at(1)},{orders2.at(2)},{orders2.at(3)}};
        

        for(int i = 0; i < orders2.size(); i +=8){
            thread t1([&]()
                {
                vector<vector<double>> a1=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i)});
                temp_stat.at(i)=a1;
            
            });
            
            thread t2([&]()
                {
                temp_stat.at(i+1)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+1)});
            });
            
            thread t3([&]()
                {
                temp_stat.at(i+2)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+2)});
            });
            
            thread t4([&]()
                {
                temp_stat.at(i+3)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+3)});
            });
            
            thread t5([&]()
                {
                temp_stat.at(i+4)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+4)});
            });
            
            thread t6([&]()
                {
                temp_stat.at(i+5)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+5)});
            });
            
            thread t7([&]()
                {
                temp_stat.at(i+6)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+6)});
            });
            
            thread t8([&]()
                {
                temp_stat.at(i+7)=abstract_search_points(dim, n1, points, spots, broken_deltas1, broken_deltas2, {orders2.at(i+7)});
            });
            
            t1.join();
            t2.join();
            t3.join();
            t4.join();
            t5.join();
            t6.join();
            t7.join();
            t8.join();
        }
        double k11=0.0;
        double k12=0.0;
        double k21=0.0;
        double k22=0.0;
        double k31=0.0;
        double k32=0.0;
        double k41=0.0;
        double k42=0.0;
        double k51=0.0;
        double k52=0.0;
        double l = temp_stat.size();
        for (int u=0; u<temp_stat.size();u++){
            k11=k11+temp_stat[u][0][0];
            k12=k12+temp_stat[u][0][1];
            k21=k21+temp_stat[u][1][0];
            k22=k22+temp_stat[u][1][1];
            k31=k31+temp_stat[u][2][0];
            k32=k32+temp_stat[u][2][1];
            k41=k41+temp_stat[u][3][0];
            k42=k42+temp_stat[u][3][1];
            k51=k51+temp_stat[u][4][0];
            k52=k52+temp_stat[u][4][1];
        }
        vector<vector<double>> centers={{k11/l,k12/l},{k21/l,k22/l},{k31/l,k32/l},{k41/l,k42/l},{k51/l,k52/l}};
        printvect(centers);
        stat.at(j)=centers;
    }
    return {{0}};
}

vector<vector<double>> choose_method(int ind,int dim, vector<vector<double>> points, vector<vector<double>> spots,vector<vector<double> > broken_ranges, vector<vector<double> > broken_deltas1,vector<vector<double> > broken_deltas2){
    
    vector<vector<double>>answ;
    int n=points.size();
    switch (ind) {
        case 1:
            answ = usual_searchpoint(dim, n, 1, points, spots, broken_deltas1, broken_deltas2);
            break;
        case 2:
            answ = usual_searchpoint(dim, n, 2, points, spots, broken_deltas1, broken_deltas2);
            break;
        case 3:
            answ = usual_searchpoint(dim, n, 3, points, spots, broken_deltas1, broken_deltas2);
            break;
        case 4:
            answ = stat_atom1(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2);
            break;
        case 5:
            answ = stat_atom2(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2);
            break;
        case 6:
            answ = stat_atom3(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2);
            break;
        case 7:
            answ = stat_atom4(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2);
        default:
            break;
    }
    return answ;
}

#endif /* stat_atoms_h */

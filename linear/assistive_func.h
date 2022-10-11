// ТУТ ВСЯКАЯ ВСПОМОГАТЕЛЬНАЯ ФИГНЯ
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#ifndef assistive_func_h
#define assistive_func_h

void printvec( vector<double> vect) {
    int n = vect.size();
    cout <<endl<< "[";
    for (int i = 0; i < n; i++)
    {
        cout << vect.at(i);
        if (i != n-1) {
            cout << ",";
        }
    }
    cout << "]"<<endl;
}

void printvect(vector<vector<double> > vect) {
    int n = vect.size();
    cout<<endl;
    for (int i = 0; i < n; i++)
    {
        cout << "[";
        int m = vect.at(i).size();
        for (int j = 0; j < m; j++)
        {
            cout << vect[i][j];
            if (j != m - 1)
            {
                cout << ",";
            }
        }
        cout << "]";
        if (i != n - 1)
        {
            cout << ",";
        }

    }
    cout<<endl;
}

void printvecto(vector<vector<int> > vect) {
    int n = vect.size();
    cout<<endl;
    for (int i = 0; i < n; i++)
    {
        cout << "(";
        int m = vect.at(i).size();
        for (int j = 0; j < m; j++)
        {
            cout << vect[i][j];
            if (j != m - 1)
            {
                cout << ",";
            }
        }
        cout << ")";
        if (i != n - 1)
        {
            cout << ",";
        }

    }
    cout<<endl;
}

void print3vect(vector<vector<vector<double> >> vect){
    int n =vect.size();
    cout<<endl<<"[";
    for (int i=0; i<n; i++) {
        printvect(vect.at(i));
        if(i<n-1){
            cout<<",";
        }
        
    }
    cout<<"]"<<endl;
}

double gaussrand(double MO, double sko)
{
    double sum = 0, x;

    for (int i = 0; i < 28; i++)
        sum += 1.0 * rand() / RAND_MAX;
    x = (sqrt(2.0) * (sko) * (sum - 14.)) / 2.11233 + MO;

    return x;
}

double ro(vector<double> dot1, vector<double> dot2) {
    int k=dot1.size();
    double r=0;
    for (int i = 0; i < k; i++)
    {
        r += pow((dot2[i] - dot1[i]), 2);
    }
    return sqrt(r);
}

vector<vector<double>>  zero_vect(int n,int dim){
    vector<vector<double>> v;
    v.resize(n);
    for (int i =0; i<n; i++) {
        vector<double> temp;
        temp.resize(dim);
        for (int ii=0; ii<dim; ii++) {
            temp.at(ii)=0;
        }
        v.at(i)=temp;
    }
    return v;
}


double abs_multivect(vector<vector<double>> multivect){
    double abs=0;
    for (int i=0; i<multivect.size(); i++) {
        for (int j=0; j<multivect.at(i).size(); j++) {
            abs+=pow(multivect.at(i).at(j),2);
        }
    }
    return sqrt(abs);
}

vector<vector<double>> v_change (double gamma,vector<vector<double>>v,vector<vector<double>>grad){
    for (int i =0; i<v.size(); i++) {
        for (int ii=0; ii<v.at(i).size(); ii++) {
            v.at(i).at(ii)=gamma*v.at(i).at(ii)+(1-gamma)*grad.at(i).at(ii);
        }
    }
    return v;
}

double cache_change (double beta,double cache,vector<vector<double>>grad){
    cache=beta*cache+(1-beta)*pow(abs_multivect(grad),2);
    return cache;
}

vector<int> set_order1(int n1, int n2){
    vector<int> order;
    order.resize(0);
    for(int i = 0; i < n1*(n1-1)/2; ++i){
        order.push_back(i % n2);
    }
    return order;
}

#endif /* assistive_func_h */

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
#include "stat_atoms.h"
using namespace std;
char up, down;
int methods;

vector<vector<double>> stat_analysis(vector<vector <vector<double>>> stat,vector <vector<double>> points){
    vector<vector<double>> results;
    vector<double> buckets;
    results.resize(2);
    buckets=zero_vect(1, 7).at(0);
    double sko=0;
    for (int i=0; i<stat.size(); i++) {
        for (int j=0; j<stat.at(i).size(); j++) {
            if(ro(stat.at(i).at(j),points.at(j))<5){
                buckets.at(0)+=1;
            }
            if(ro(stat.at(i).at(j),points.at(j))<10){
                buckets.at(1)+=1;
            }
            if(ro(stat.at(i).at(j),points.at(j))<20){
                buckets.at(2)+=1;
            }
            if(ro(stat.at(i).at(j),points.at(j))<30){
                buckets.at(3)+=1;
            }
            if(ro(stat.at(i).at(j),points.at(j))<50){
                buckets.at(4)+=1;
            }
            if(ro(stat.at(i).at(j),points.at(j))<100){
                buckets.at(5)+=1;
            }
            if(ro(stat.at(i).at(j),points.at(j))<200){
                buckets.at(6)+=1;
            }
            for(int k=0;k<stat.at(i).at(j).size();k++){
                sko+=pow(stat.at(i).at(j).at(k)-points.at(j).at(k),2);
            }
        }
    }
    results.at(0)={sqrt(sko/(stat.size()*stat.at(0).size()))};
    results.at(1)=buckets;
    return results;
}

void usual_stat(double SKO, int number, int dim, int n, vector<vector<double> > points,  vector<vector<double> > spots,vector<vector<double>> ranges, vector<vector<double>> deltas1,vector<vector<double>> deltas2,vector<double> delays1,vector<double> delays2){
    vector<int> method(0);
    int tmp=methods;
    while(tmp>0){
        method.push_back(tmp%10);
        tmp=tmp/10;
    }
    
    vector<vector<vector <vector<double>>>> multistat(7);
    /*vector<vector <vector<double>>> statistics1;
    statistics1.resize(number);
    vector<vector <vector<double>>> statistics2;
    statistics2.resize(number);
    vector<vector <vector<double>>> statistics3;
    statistics3.resize(number);
    vector<vector <vector<double>>> statistics4;
    statistics4.resize(number);
    vector<vector <vector<double>>> statistics5;
    statistics5.resize(number);
    vector<vector <vector<double>>> statistics6;
    statistics6.resize(number);*/
    for (int i=0; i<number; i++) {

        
        cout<<endl<<i;
        vector<vector<double> > broken_deltas1=set_full_deltas_correct(points, spots, SKO, delays1, delays2);
        vector<vector<double> > broken_deltas2=set_full_deltas_correct(spots, points, SKO, delays2, delays1);
        vector<vector<double> > broken_ranges=break_ranges_with_delays(ranges, SKO, delays1);
        for (auto & l: method) {
            cout<<endl<<l;
            multistat.at(l-1).push_back(choose_method(l, dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2));
        }
        //statistics1.at(i)=usual_searchpoint(dim, n, 1, points, spots, broken_deltas1, broken_deltas2);// прямой метод
        //print3vect(statistics1);
        //statistics2.at(i)=usual_searchpoint(dim, n, 2, points, spots, broken_deltas1, broken_deltas2);// обратный метод
        //print3vect(statistics2);
        //statistics3.at(i)=usual_searchpoint(dim, n, 3, points, spots, broken_deltas1, broken_deltas2);// совмещенный
        //print3vect(statistics3);
        //statistics4.at(i)=stat_atom1(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2); // совмещенный, но известны расстояния сверху, снизу СЕВ
        //print3vect(statistics4);
        //statistics5.at(i)=stat_atom2(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2); // совмещенный, но известны расстояния сверху, снизу нет СЕВ
        //print3vect(statistics5);
        //statistics6.at(i)=stat_atom3(dim, points, spots, broken_ranges, broken_deltas1, broken_deltas2);
        //statistics6.at(i)=abstract_search_points(dim, n, points, spots, broken_deltas1, broken_deltas2, {{0,1,2,3,4,4,3,2,1,0}}); // разность разностей
        

    }
    for (int i = 0; i<multistat.size(); i++) {
        print3vect(multistat.at(i));
    }
    /*print3vect(statistics1);
    cout<<endl<<endl<<endl<<endl<<endl;
    print3vect(statistics2);
    cout<<endl<<endl<<endl<<endl<<endl;
    print3vect(statistics3);
    cout<<endl<<endl<<endl<<endl<<endl;
    print3vect(statistics4);
    cout<<endl<<endl<<endl<<endl<<endl;
    print3vect(statistics5);
    cout<<endl<<endl<<endl<<endl<<endl;
    print3vect(statistics6);
    cout<<endl<<endl<<endl<<endl<<endl;
    printvect(stat_analysis(statistics1, points));
    printvect(stat_analysis(statistics2, points));
    printvect(stat_analysis(statistics3, points));
    printvect(stat_analysis(statistics4, points));
    printvect(stat_analysis(statistics5, points));
    printvect(stat_analysis(statistics6, points));*/
}




int main() {
    srand((unsigned)time(NULL));
    int n1;
    int n2;
    int times;
    double SKO;
    double radius1;
    double radius2;
    double length;
    double angle1;
    double angle2;
    int dim=2;
    
    
    /*cout<<"введите количество пунктов";
    cin>>n2;
    cout<<"введите количество неизвестных точек";
    cin>>n1;
    cout<<"введите радиус круга с пунктами, м";
    cin>>radius2;
    cout<<"введите радиус круга с неизвестными точками, м";
    cin>>radius1;
    cout<<"введите расстояние между центрами кругов, м";
    cin>>length;
    cout<<"введите угол для круга с пунктами, радианы";
    cin>>angle2;
    cout<<"введите угол для круга с неизвестными точками, радины";
    cin>>angle1;
    cout<<"введите СКО измерений, метры";
    cin>>SKO;
    cout<<"введите количество прогонов для статистики";
    cin>>times;*/
    n1=5;
    n2=5;
    radius1=30000;
    radius2=40000;
    length=40000;
    angle1=2;
    angle2=2;
    SKO=3;
    times=1;
    vector<double> spacepoint={40000,0};
    vector <vector<double>> points=arc_points_init(n1, dim, radius1, {0,length}, angle1, {-1,0});
    vector <vector<double>> spots=arc_points_init(n2, dim, radius2, {0,0}, angle2, {-1,0});
        
    vector<vector<double> > deltas1=set_full_deltas(dim, points, spots);
    vector<vector<double> > deltas2=set_full_deltas(dim, spots, points);
    vector<vector<double> > ranges =set_ranges(dim, points);

    vector<double> delays1;
    vector<double> delays2;
    // Далее идут 4 варианта, выбрать один
    std::cout << "Введите 1, если существует СЕВ сверху. Иначе введите 0";
    std::cin >> up;
    std::cout << "Введите 1, если существует СЕВ снизу. Иначе введите 0";
    std::cin >> down;
    if (up == '1') {
        delays1=set_zero_delays(n1);
    } else {
        delays1=set_delays(n1);
    }
    if (down == '1') {
        delays2=set_zero_delays(n2);
    } else {
        delays2=set_delays(n2);
    }
    cout<<endl<<"Выберите методы поиска для их сравнения: "<<endl;
    if( (up == '1') && (down == '1') ){
        cout<<"Чтобы добавить прямой метод введите 1"<<endl;
        cout<<"Чтобы добавить обратный метод введите 2"<<endl;
        cout<<"Чтобы добавить объединение прямого и обратного методов введите 3"<<endl;
        cout<<"Чтобы добавить объединение прямого и обратного методов, использующее расстояния сверху, введите 4"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей и расстояния сверху, введите 5"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей без расстояний сверху, введите 6"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей снизу и просто разности сверху, введите 7"<<endl;
        cout<<"Пример ввода: 1236 (все цифры подряд без пробелов в порядке возрастания)"<<endl;
    }
    if( (up == '0') && (down == '0') ){
        cout<<"Чтобы добавить метод, использующий разности разностей и расстояния сверху, введите 5"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей без расстояний сверху, введите 6"<<endl;
        cout<<"Пример ввода: 56 (все цифры подряд без пробелов в порядке возрастания)"<<endl;
    }
    
    if( (up == '0') && (down == '1') ){
        cout<<"Чтобы добавить объединение прямого и обратного методов, использующее расстояния сверху, введите 4"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей и расстояния сверху, введите 5"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей без расстояний сверху, введите 6"<<endl;
        cout<<"Пример ввода: 456 (все цифры подряд без пробелов в порядке возрастания)"<<endl;
    }
    
    if( (up == '1') && (down == '0') ){
// еще один метод тут
        cout<<"Чтобы добавить метод, использующий разности разностей и расстояния сверху, введите 5"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей без расстояний сверху, введите 6"<<endl;
        cout<<"Пример ввода: 1236 (все цифры подряд без пробелов в порядке возрастания)"<<endl;
        cout<<"Чтобы добавить метод, использующий разности разностей снизу и просто разности сверху, введите 7"<<endl;
    }
    cin>>methods;
    /*//нигде нет СЕВ
    delays1=set_delays(n1);
    delays2=set_delays(n2);
    
    //СЕВ снизу
    delays1=set_delays(n1);
    delays2=set_zero_delays(n2);
    
    //СЕВ сверху
    delays1=set_zero_delays(n1);
    delays2=set_delays(n2);
    
    //везде есть СЕВ
    delays1=set_zero_delays(n1);
    delays2=set_zero_delays(n2);*/
    
    
    usual_stat(SKO, times, dim, n1, points, spots, ranges, deltas1, deltas2, delays1, delays2);
}


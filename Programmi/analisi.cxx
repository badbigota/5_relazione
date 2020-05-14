#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include "functions.h"
using namespace std;

int main()
{
    vector<data_torque> torq;
    vector<x_y> periodi;
    string map_file = "../Dati/mappa.txt";
    load_data(map_file, torq);
    get_zero_time(torq, periodi);
    
    for (auto d : torq)//controllo per i dati
    {
        cout << d.time.size() << "\t" << d.forcing.size() << "\t" << d.a.size() << endl;
    }
    for(auto d: periodi){
        cout<<d.time.size()<<endl;
    }
}
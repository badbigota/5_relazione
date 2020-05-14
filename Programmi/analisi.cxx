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
    vector<x_y> tempi;
    vector<x_y> periodi;
    string map_file = "../Dati/mappa.txt";
    load_data(map_file, torq);
    get_zero_time(torq, tempi);
    get_periods(tempi, periodi);

    //Stampa per prova
    for(int i=0; i<periodi[3].time.size(); i++){
        cout<<periodi[0].time[i]<<"\t"<<periodi[1].time[i]<<"\t"<<periodi[4].time[i]<<endl;
    }

    vector<point_regime> punti_grafico_resonance;
    
}
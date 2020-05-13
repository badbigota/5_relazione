#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include "functions.h"
using namespace std;

int main()
{
    vector<data_torque> tor;
    string map_file = "../Dati/mappa.txt";
    load_data(map_file, tor);

    
    for (auto d : tor)//controllo per i dati
    {
        cout << d.time.size() << "\t" << d.forcing.size() << "\t" << d.a.size() << endl;
    }
}
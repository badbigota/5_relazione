#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
using namespace std;

//Struttura dati grezzi
struct data_torque
{
    int data_file_number;  //numero progressivo
    double data_file_freq; //freq fornita
    vector<double> time;
    vector<double> a;
    vector<double> forcing;
};

//Carica tutti i dati grezzi in strutture contenute in un vettore
void load_data(string map_file, vector<data_torque> &vec_data)
{
    ifstream f_map(map_file);
    int num;
    string file_name;
    double freq;
    while (f_map >> num >> file_name >> freq)
    {
        data_torque temp;
        temp.data_file_number = num;
        temp.data_file_freq = freq;
        ifstream f_data_points("../Dati/" + file_name);
        double time, forcing, a;
        while (f_data_points >> time >> forcing >> a)
        {
            temp.time.push_back(time);
            temp.a.push_back(a);
            temp.forcing.push_back(forcing);
        }
        vec_data.push_back(temp);
    }
}
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include "statistica.h"
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
//Struttura varie tipologie di dato x,y a volte neppure usate entrambe
struct x_y
{
    double freq;
    vector<int> number;
    vector<double> time;
    vector<double> amplitude;
};

struct point_regime{
    double omega_regime;
    double amplitude_regime;
    double err_omega_regime;
    double err_amplitude_regime;
    double omega_th;
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

//Fa interpolazione per prendere la "radice", ovvero quando la sinusoide intercetta asse x
double get_root(int &i, data_torque &d)
{
    vector<double> x = {d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2]};
    vector<double> y = {d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2]};
    double root = -a_intercetta_err_uguali(x, y) / b_angolare_err_uguali(x, y); //trova intersezione retta con asse x
    return root;
}
//Fa interpolazione per prendere la "radice", ovvero quando la sinusoide intercetta asse x ma prende 5 dati piuttosto che 4 poichè voglio come polo lo 0
double get_root_per_0(int &i, data_torque &d)
{
    vector<double> x = {d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
    vector<double> y = {d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.time[i + 3]};
    double root = -a_intercetta_err_uguali(x, y) / b_angolare_err_uguali(x, y); //trova intersezione retta con asse x
    return root;
}
//Genera il vettore di struttura (una per torq) dove ciascuna strutura ha il vettore con il momento di intercetta della sinusoide con lo zero
void get_zero_time(vector<data_torque> &vec_data, vector<x_y> &time_amp)
{
    for (auto d : vec_data)
    {
        x_y zeros_torq;
        zeros_torq.freq = d.data_file_freq;
        for (int i = 0; i < d.a.size(); i++)
        {
            if (d.a[i] > 0)
            {
                if (d.a[i + 1] <= 0 && d.a[i + 2] <= 0 && d.a[i + 3] <= 0) //both of the opposite sign
                {
                    if (d.a[i + 1] == 0)
                    {
                        zeros_torq.time.push_back(get_root_per_0(i, d));
                    }
                    else
                    {
                        zeros_torq.time.push_back(get_root(i, d));
                    }
                }
            }
            else if (d.a[i] < 0)
            {
                if (d.a[i + 1] >= 0 && d.a[i + 2] >= 0 && d.a[i + 3] >= 0) //both of the opposite sign
                {
                    if (d.a[i + 1] == 0)
                    {
                        zeros_torq.time.push_back(get_root_per_0(i, d));
                    }
                    else
                    {
                        zeros_torq.time.push_back(get_root(i, d));
                    }
                }
            }
        }
        time_amp.push_back(zeros_torq);
    }
}

//Genera tutti i periodi prendendo differenze indipendenti
void get_periods(vector<x_y> &times, vector<x_y> &periods)
{
    for (int i = 0; i < times.size(); i++)
    {
        x_y temp;
        temp.freq = times[i].freq;
        for (int j = 0; j < times[i].time.size() - 3; j = j + 3) //meno 3 per evitare out of range o simili
        {
            temp.time.push_back(times[i].time[j + 2] - times[i].time[j]);
        }
        periods.push_back(temp);
    }
}

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include "statistica.h"
using namespace std;

/*
void funzione(lettura, scrittura)
*/

//Struttura dati grezzi
struct data_campione
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
    double avg_time;
    double err_avg_time;
    double avg_amplitude;
    double err_avg_amplitude;
    vector<int> closest_index_zero;
    vector<double> time;
    vector<double> amplitude;
    vector<double> omega_sperim;
    double omega_media1;
    double omega_media2;
    double err_omega_media1;
    double err_omega_media2;
};

struct punti_massimo
{
    vector<double> t_max;
    vector<double> ampl_max;
    vector<double> index; //per una veloce referenza ai dati grezzi
    vector<double> t_ref; //idem sopra
};

struct vettoredoppio
{
    vector<double> vettore2; //per l'indice del massimo
    vector<double> vettore3; //per il valore del massimo
};

//Carica tutti i dati grezzi in strutture contenute in un vettore
void load_data(string map_file, vector<data_campione> &vec_data)
{
    ifstream f_map(map_file);
    int num;
    string file_name;
    double freq;
    while (f_map >> num >> file_name >> freq)
    {
        data_campione temp;
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
double get_root(int &i, data_campione &d)
{

    vector<double> x; // = {d.time[i - 3], d.time[i - 2], d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
    vector<double> y; // = {d.a[i - 3], d.a[i - 2], d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.a[i + 3]};
    if (i > 3)
    {
        x = {d.time[i - 3], d.time[i - 2], d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
        y = {d.a[i - 3], d.a[i - 2], d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.a[i + 3]};
    }
    else
    {
        x = {d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
        y = {d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.a[i + 3]};
    }

    double root = -a_intercetta_err_uguali(x, y) / b_angolare_err_uguali(x, y); //trova intersezione retta con asse x
    return root;
}
//Fa interpolazione per prendere la "radice", ovvero quando la sinusoide intercetta asse x ma prende 5 dati piuttosto che 4 poichè voglio come polo lo 0
double get_root_per_0(int &i, data_campione &d)
{
    vector<double> x = {d.time[i - 2], d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
    vector<double> y = {d.a[i - 2], d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.time[i + 3]};
    double root = -a_intercetta_err_uguali(x, y) / b_angolare_err_uguali(x, y); //trova intersezione retta con asse x
    return root;
}
//Genera il vettore di struttura (una per campione) dove ciascuna strutura ha il vettore con il momento di intercetta della sinusoide con lo zero
void get_zero_time(vector<data_campione> &vec_data, vector<x_y> &time_amp)
{
    for (auto d : vec_data)
    {
        x_y zeros_campione;
        zeros_campione.freq = d.data_file_freq;
        for (int i = 0; i < d.a.size() - 3; i++)
        {
            if (d.a[i] > 0)
            {
                if (d.a[i + 1] <= 0 && d.a[i + 2] <= 0 && d.a[i + 3] <= 0) //both of the opposite sign
                {
                    if (d.a[i + 1] == 0)
                    {
                        zeros_campione.time.push_back(get_root_per_0(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                    else
                    {
                        zeros_campione.time.push_back(get_root(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                }
            }
            else if (d.a[i] < 0)
            {
                if (d.a[i + 1] >= 0 && d.a[i + 2] >= 0 && d.a[i + 3] >= 0) //both of the opposite sign
                {
                    if (d.a[i + 1] == 0)
                    {
                        zeros_campione.time.push_back(get_root_per_0(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                    else
                    {
                        zeros_campione.time.push_back(get_root(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                }
            }
        }
        time_amp.push_back(zeros_campione);
    }
}

//Genera tutti i periodi prendendo differenze indipendenti
void get_periods(vector<x_y> &times, vector<x_y> &periods)
{
    for (int i = 0; i < times.size(); i++)
    {
        x_y temp;
        temp.freq = times[i].freq;                               // si ricorda la frequenza, ovvero quali colonne sto leggendo
        for (int j = 0; j < times[i].time.size() - 3; j = j + 3) //meno 3 per evitare out of range o simili
        {
            temp.time.push_back(times[i].time[j + 2] - times[i].time[j]);
        }
        periods.push_back(temp);
    }
}

void get_periodi_medi(vector<x_y> &periodi, vector<x_y> &media_periodi) //genra la media dei tempi e errore
{
    for (auto d : periodi)
    {
        x_y temp;
        temp.avg_time = media(d.time);
        temp.err_avg_time = dstd_media(d.time);
        media_periodi.push_back(temp);
    }
}

//Trova gli indici di massimo in assoluto, utili anche poi per fare picco picco automaticamente senza approx
void get_index_maxima(vector<data_campione> &dati_grezzi, vector<x_y> &tempi, vector<vettoredoppio> &indici_max_per_campioni)
{
    for (int i = 0; i < dati_grezzi.size(); i++) //per ogni cmapione
    {
        vector<double> &ampiezze_grezze = dati_grezzi[i].a;
        vector<int> &tempi_zero = tempi[i].closest_index_zero;
        vettoredoppio massimi;
        for (int j = 0; j < tempi_zero.size() - 1; j++)
        {
            int index_max = 0;
            double val_max = 0;
            for (int num = tempi_zero[j]; num < tempi_zero[j + 1]; num++)
            {
                if (abs(ampiezze_grezze[num]) > val_max)
                {
                    val_max = ampiezze_grezze[num];
                    index_max = num;
                }
            }
            massimi.vettore2.push_back(index_max);
            massimi.vettore3.push_back(val_max);
        }
        indici_max_per_campioni.push_back(massimi);
    }
}

//Somme con potenze per massima verosimiglianza
double sum_pow(vector<double> x, int power)
{
    double sum_temp = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum_temp += pow(x[i], power);
    }
    return sum_temp;
}

double sum_pow_special(vector<double> x, vector<double> y, int power_x, int power_y)
{
    double sum_temp;
    for (int i = 0; i < x.size(); i++)
    {
        sum_temp += pow(x[i], power_x) * pow(y[i], power_y);
    }
    return sum_temp;
}

//Massima verosimiglianza
void best_fit_mv(vector<double> x, vector<double> y, vector<double> &best_fit)
{
    double a, b, c, d, e, f, g, h;
    a = sum_pow(x, 4);
    b = sum_pow(x, 3);
    c = sum_pow(x, 2);
    d = sum_pow(x, 1);
    e = sum_pow(x, 0);
    f = sum_pow_special(x, y, 2, 1);
    g = sum_pow_special(x, y, 1, 1);
    h = sum_pow_special(x, y, 0, 1);
    double num_a = -pow(d, 2) * f + c * e * f + c * d * g - b * e * g - pow(c, 2) * h + b * d * h;
    double num_b = c * d * f - b * e * f - pow(c, 2) * g + a * e * g + b * c * h - a * d * h;
    double num_c = -pow(c, 2) * f + b * d * f + b * c * g - a * d * g - pow(b, 2) * h + a * c * h;
    double den = pow(c, 3) - 2 * b * c * d + a * pow(d, 2) + pow(b, 2) * e - a * c * e;
    best_fit.push_back(-num_a / den);
    best_fit.push_back(-num_b / den);
    best_fit.push_back(-num_c / den);
}

////Punto di massimo per fit
double max_x_fit(vector<double> parametri_fit)
{
    double a = parametri_fit[0];
    double b = parametri_fit[1];
    double c = parametri_fit[2];
    double max_val = -b / (2 * a);
    return max_val;
}

//////Valore punto di massimo per fit
double max_y_fit(vector<double> parametri_fit)
{
    double a = parametri_fit[0];
    double b = parametri_fit[1];
    double c = parametri_fit[2];
    double max_val = -b / (2 * a);
    double y_max = a * pow(max_val, 2) + b * max_val + c;
    return y_max;
}

//Genera un vettore di punti_massimo struct con tutti i punti di massimo e minimo per ciascun campione
void get_maxima(vector<data_campione> &raw_data, vector<vettoredoppio> &index_maxima, vector<punti_massimo> &maxima)
{
    for (int i = 0; i < raw_data.size(); i++)
    {
        punti_massimo temp_punti_max;
        vector<double> &tempi_grezzi = raw_data[i].time;
        vector<double> &ampiezze_grezze = raw_data[i].a;
        vector<double> &indici_massimi = index_maxima[i].vettore2;
        for (int k = 0; k < indici_massimi.size(); k++) //per ogni indice masimo
        {
            double pivot = indici_massimi[k]; //effettivo punto di massimo
            vector<double> x, y;
            vector<double> best_fitting_par;
            for (int num = -9; num <= 9; num++) //prende i nove prima e i nove dopo
            {
                x.push_back(tempi_grezzi[num + pivot]); //genera i vettori da dare in pasto alla mv
                y.push_back(ampiezze_grezze[num + pivot]);
            }
            best_fit_mv(x, y, best_fitting_par);                            //genera i valori di interpolazione parabola
            temp_punti_max.t_max.push_back(max_x_fit(best_fitting_par));    //trova punto di max
            temp_punti_max.ampl_max.push_back(max_y_fit(best_fitting_par)); //trova valore punto di max
            temp_punti_max.index.push_back(pivot);                          //per un traceback accurato
            temp_punti_max.t_ref.push_back(tempi_grezzi[pivot]);            //per un traceback accurato
        }
        maxima.push_back(temp_punti_max);
    }
}

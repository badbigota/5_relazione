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
    vector<double> picco_picco;
    double omega_media1;
    double omega_media2;
    double err_omega_media1;
    double err_omega_media2;
    double picco_medio;
    double err_picco_medio;
    vector<double> offset;
};

struct punti_massimo
{
    double freq_ref;
    vector<double> t_max;
    vector<double> ampl_max;
    vector<double> index; //per una veloce referenza ai dati grezzi
    vector<double> t_ref; //idem sopra
    vector<double> a_ref;
    vector<double> coeff_a;
    vector<double> coeff_b;
    vector<double> coeff_c;
};

struct vettoredoppio
{
    vector<double> vettore2; //per l'indice del massimo
    vector<double> vettore3; //per il valore del massimo
};

struct punto_regime //perchè useremo questa per avere tutti i puti per il grafico
{
    double omega;
    double err_omega;
    double theta;
    double err_theta;
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
            temp.a.push_back(a * 2. * M_PI);
            temp.forcing.push_back(forcing * 1000. * 2. * M_PI);
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
//Adattato interamente da quello di marc e fab in analisi cxx
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
                    val_max = abs(ampiezze_grezze[num]);
                    index_max = num;
                }
            }
            massimi.vettore2.push_back(index_max);
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
    long double num_a = -pow(d, 2) * f + c * e * f + c * d * g - b * e * g - pow(c, 2) * h + b * d * h;
    long double num_b = c * d * f - b * e * f - pow(c, 2) * g + a * e * g + b * c * h - a * d * h;
    long double num_c = -pow(c, 2) * f + b * d * f + b * c * g - a * d * g - pow(b, 2) * h + a * c * h;
    long double den = pow(c, 3) - 2 * b * c * d + a * pow(d, 2) + pow(b, 2) * e - a * c * e;
    best_fit.push_back(-num_a / den);
    best_fit.push_back(-num_b / den);
    best_fit.push_back(-num_c / den);
}

////Punto di massimo per fit
double max_x_fit(vector<double> parametri_fit)
{
    long double a = parametri_fit[0];
    long double b = parametri_fit[1];
    long double c = parametri_fit[2];
    long double max_val = -b / (2 * a);
    return max_val;
}

//////Valore punto di massimo per fit
double max_y_fit(vector<double> parametri_fit)
{
    long double a = parametri_fit[0];
    long double b = parametri_fit[1];
    long double c = parametri_fit[2];
    long double max_val = -b / (2 * a);
    long double y_max = a * pow(max_val, 2) + b * max_val + c;
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
            temp_punti_max.a_ref.push_back(ampiezze_grezze[pivot]);
            temp_punti_max.freq_ref = raw_data[i].data_file_freq; //per un traceback
            temp_punti_max.coeff_a.push_back(best_fitting_par[0]);
            temp_punti_max.coeff_b.push_back(best_fitting_par[1]);
            temp_punti_max.coeff_c.push_back(best_fitting_par[2]);
        }
        maxima.push_back(temp_punti_max);
    }
}

void offset(vector<punti_massimo> vec_punti_max, vector<x_y> &vec_picco_picco)
{
    for (int i = 0; i < vec_punti_max.size(); i++)
    {                                                          //per ogni campione
        vector<double> &ampiezze_max = vec_punti_max[i].a_ref; //su i punti grezzi
        x_y temp_picchi;
        for (int j = 0; j < ampiezze_max.size() - 1; j = j + 2) //altirmenti core dumped
        {
            double offset = (ampiezze_max[j] + ampiezze_max[j + 1]) / 2.;
            temp_picchi.offset.push_back(offset);
        }
        temp_picchi.freq = vec_punti_max[i].freq_ref;
        vec_picco_picco.push_back(temp_picchi);
    }
}

//Picco picco con valori massimi assoluti
void picco_picco_ass(vector<punti_massimo> vec_punti_max, vector<x_y> &vec_picco_picco)
{
    for (int i = 0; i < vec_punti_max.size(); i++)
    {                                                          //per ogni campione
        vector<double> &ampiezze_max = vec_punti_max[i].a_ref; //su i punti grezzi
        x_y temp_picchi;
        for (int j = 0; j < ampiezze_max.size() - 1; j = j + 2) //altirmenti core dumped
        {
            double picco_picco_ = abs(ampiezze_max[j] - ampiezze_max[j + 1]);
            temp_picchi.picco_picco.push_back(picco_picco_ / 2.);
        }
        temp_picchi.freq = vec_punti_max[i].freq_ref;
        vec_picco_picco.push_back(temp_picchi);
    }
}

//Picco picco a partire da massima verosimiglianza
void picco_picco_mv(vector<punti_massimo> vec_punti_max, vector<x_y> &vec_picco_picco)
{
    for (int i = 0; i < vec_punti_max.size(); i++)
    {                                                             //per ogni campione
        vector<double> &ampiezze_max = vec_punti_max[i].ampl_max; //su i punti grezzi
        x_y temp_picchi;
        for (int j = 0; j < ampiezze_max.size() - 1; j = j + 2) //altirmenti core dumped
        {
            double picco_picco_ = abs(ampiezze_max[j] - ampiezze_max[j + 1]);
            temp_picchi.picco_picco.push_back(picco_picco_ / 2.);
        }
        temp_picchi.freq = vec_punti_max[i].freq_ref;
        vec_picco_picco.push_back(temp_picchi);
    }
}

void picco_picco_chi() {}

//Calcola la media e la dstd dei picco picco mezzi
void add_picco_medio(vector<x_y> picchi, vector<punto_regime> &punti_regime)
{
    punto_regime temp_punto_reg;
    for (auto d : picchi)
    {
        temp_punto_reg.theta = media(d.picco_picco, 0, 20);
        temp_punto_reg.err_theta = dstd_media(d.picco_picco, 0, 20);
        punti_regime.push_back(temp_punto_reg);
    }
}

//Calcola gli omega secondo JCGM. METTERE SOLO DOPO CALCOLO DI THETA
void add_omega_2(vector<x_y> periodi, vector<punto_regime> &punti_regime)
{
    for (int i = 0; i < periodi.size(); i++)
    {
        vector<double> temp_omega;
        for (auto d : periodi[i].time)
        {
            temp_omega.push_back(2. * M_PI / d);
        }

        punti_regime[i].omega = media(temp_omega);
        punti_regime[i].err_omega = dstd_media(temp_omega);
    }
}
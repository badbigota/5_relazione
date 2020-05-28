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
    double err_time;
    double num_misure;
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
    vector<double> err_coeff_a;
    vector<double> err_coeff_b;
    vector<double> err_coeff_c;
    vector<double> err_ampl_max; //solo per la linearizzazione
};

struct vettoredoppio
{
    vector<double> vettore2; //per l'indice del massimo
    vector<double> vettore3; //per il valore del massimo
    double dispersione_amp;
};

struct punto_regime //perchè useremo questa per avere tutti i puti per il grafico
{
    double omega;
    double err_omega;
    double theta;
    double err_theta;
};

struct interpolazione_gamma
{
    double a_intercetta_gamma_max;
    double b_angolare_gamma_max;
    double a_intercetta_gamma_min;
    double b_angolare_gamma_min;
    double r_max;
    double r_min;
    double t_max;
    double t_min;
    double err_a_post_max;
    double err_b_post_max;
    double err_a_post_min;
    double err_b_post_min;
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

//Carica i dati grezzi del decadimento
void load_data_decay(vector<data_campione> &vec_data)
{
    for (int i = 0; i < 3; i++)
    {
        data_campione temp_data;
        ifstream fin_smorz("../Dati/s_" + to_string(i + 1) + ".csv");
        double time, forcing, amp;
        while (fin_smorz >> time >> forcing >> amp)
        {
            temp_data.time.push_back(time);
            temp_data.forcing.push_back(forcing * 1000. * 2. * M_PI);
            temp_data.a.push_back(amp * 2. * M_PI);
        }
        vec_data.push_back(temp_data);
    }
}

//Legge i parametri di root per la fase a regime
void lettura_parametri(vector<punti_massimo> &parametri)
{
    for (int i = 0; i < 20; i++)
    {
        punti_massimo temp_p_max;
        double a, b, c, d, e, f, g, h;
        ifstream fin_param("../Dati/parameters" + to_string(i + 1) + ".txt");
        while (fin_param >> a >> b >> c >> d >> e >> f >> g >> h)
        {
            temp_p_max.coeff_c.push_back(c);
            temp_p_max.err_coeff_c.push_back(d);
            temp_p_max.coeff_b.push_back(e);
            temp_p_max.err_coeff_b.push_back(f);
            temp_p_max.coeff_a.push_back(g);
            temp_p_max.err_coeff_a.push_back(h);
        }
        parametri.push_back(temp_p_max);
    }
}

//Legge i parametri di smorzamento da root
void lettura_parametri_smorz(vector<punti_massimo> &parametri)
{
    for (int i = 0; i < 3; i++)
    {
        punti_massimo temp_p_max;
        double a, b, c, d, e, f, g, h;
        ifstream fin_param("../Dati/parameters_smorz_" + to_string(i) + ".txt");
        while (fin_param >> a >> b >> c >> d >> e >> f >> g >> h)
        {
            temp_p_max.coeff_c.push_back(c);
            temp_p_max.err_coeff_c.push_back(d);
            temp_p_max.coeff_b.push_back(e);
            temp_p_max.err_coeff_b.push_back(f);
            temp_p_max.coeff_a.push_back(g);
            temp_p_max.err_coeff_a.push_back(h);
        }
        parametri.push_back(temp_p_max);
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
    vector<double> y = {d.a[i - 2], d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.a[i + 3]};
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
/*
//Presa dei tempi su intercette di forzante
*/
double get_root_forcing(int &i, data_campione &d)
{

    vector<double> x; // = {d.time[i - 3], d.time[i - 2], d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
    vector<double> y; // = {d.a[i - 3], d.a[i - 2], d.a[i - 1], d.a[i], d.a[i + 1], d.a[i + 2], d.a[i + 3]};
    if (i >= 1)
    {
        x = {d.time[i - 3], d.time[i - 2], d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
        y = {d.forcing[i - 3], d.forcing[i - 2], d.forcing[i - 1], d.forcing[i], d.forcing[i + 1], d.forcing[i + 2], d.forcing[i + 3]};
    }
    else
    {
        x = {d.time[i], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3], d.time[i + 4]};
        y = {d.forcing[i], d.forcing[i], d.forcing[i + 1], d.forcing[i + 2], d.forcing[i + 3], d.forcing[i + 4]};
    }

    double root = -a_intercetta_err_uguali(x, y) / b_angolare_err_uguali(x, y); //trova intersezione retta con asse x
    return root;
}
//Fa interpolazione per prendere la "radice", ovvero quando la sinusoide intercetta asse x ma prende 5 dati piuttosto che 4 poichè voglio come polo lo 0
double get_root_per_0_forcing(int &i, data_campione &d)
{
    vector<double> x = {d.time[i - 2], d.time[i - 1], d.time[i], d.time[i + 1], d.time[i + 2], d.time[i + 3]};
    vector<double> y = {d.forcing[i - 2], d.forcing[i - 1], d.forcing[i], d.forcing[i + 1], d.forcing[i + 2], d.forcing[i + 3]};
    double root = -a_intercetta_err_uguali(x, y) / b_angolare_err_uguali(x, y); //trova intersezione retta con asse x
    return root;
}

//Genera il vettore di struttura (una per campione) dove ciascuna strutura ha il vettore con il momento di intercetta della sinusoide con lo zero
void get_zero_time_forcing(vector<data_campione> &vec_data, vector<x_y> &time_amp)
{
    for (auto d : vec_data)
    {
        x_y zeros_campione;
        zeros_campione.freq = d.data_file_freq;
        for (int i = 0; i < d.forcing.size() - 3; i++)
        {
            if (d.forcing[i] > 0)
            {
                if (d.forcing[i + 1] <= 0 && d.forcing[i + 2] <= 0 && d.forcing[i + 3] <= 0) //both of the opposite sign
                {
                    if (d.forcing[i + 1] == 0)
                    {
                        zeros_campione.time.push_back(get_root_per_0_forcing(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                    else
                    {
                        zeros_campione.time.push_back(get_root_forcing(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                }
            }
            else if (d.forcing[i] < 0)
            {
                if (d.forcing[i + 1] >= 0 && d.forcing[i + 2] >= 0 && d.forcing[i + 3] >= 0) //both of the opposite sign
                {
                    if (d.forcing[i + 1] == 0)
                    {
                        zeros_campione.time.push_back(get_root_per_0_forcing(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                    else
                    {
                        zeros_campione.time.push_back(get_root_forcing(i, d));
                        zeros_campione.closest_index_zero.push_back(i);
                        zeros_campione.amplitude.push_back(d.time[i]); //solo per testare su file
                    }
                }
            }
        }
        time_amp.push_back(zeros_campione);
    }
}

/*FINE USO DI FORCING
*/

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

//genera la media dei tempi e errore
void get_periodi_medi(vector<x_y> &periodi, vector<x_y> &media_periodi)
{
    for (auto d : periodi)
    {
        x_y temp;
        temp.avg_time = media(d.time);
        temp.err_avg_time = dstd_media(d.time);
        temp.err_time = dstd(d.time);
        temp.num_misure = d.time.size();
        media_periodi.push_back(temp);
    }
}

//Trova gli indici di massimo in assoluto, utili anche poi per fare picco picco automaticamente senza approx
void get_index_maxima(vector<data_campione> &dati_grezzi, vector<x_y> &tempi, vector<vettoredoppio> &indici_max_per_campioni)
{
    for (int i = 0; i < dati_grezzi.size(); i++) //per ogni campione
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

//Punto di massimo per fit
double max_x_fit(vector<double> parametri_fit)
{
    long double a = parametri_fit[0];
    long double b = parametri_fit[1];
    long double c = parametri_fit[2];
    long double max_val = -b / (2 * a);
    return max_val;
}

//Valore punto di massimo per fit
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
void get_maxima_ass(vector<data_campione> &raw_data, vector<vettoredoppio> &index_maxima, vector<punti_massimo> &maxima)
{
    for (int i = 0; i < raw_data.size(); i++)
    {
        punti_massimo temp_punti_max;
        vector<double> &tempi_grezzi = raw_data[i].time;
        vector<double> &ampiezze_grezze = raw_data[i].a;
        vector<double> &indici_massimi = index_maxima[i].vettore2;
        for (int k = 0; k < indici_massimi.size(); k++) //per ogni indice masimo
        {
            double pivot = indici_massimi[k];                    //effettivo punto di massimo
            temp_punti_max.index.push_back(pivot);               //per un traceback accurato
            temp_punti_max.t_max.push_back(tempi_grezzi[pivot]); //per un traceback accurato
            temp_punti_max.ampl_max.push_back(ampiezze_grezze[pivot]);
            temp_punti_max.freq_ref = raw_data[i].data_file_freq; //per un traceback
        }
        maxima.push_back(temp_punti_max);
    }
}

//Genera i massimi tramite massima verosimiglianza intorno agli indici di massimo
void get_maxima_mv(vector<data_campione> &raw_data, vector<vettoredoppio> &index_maxima, vector<punti_massimo> &maxima)
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
            temp_punti_max.freq_ref = raw_data[i].data_file_freq;           //per un traceback
            temp_punti_max.coeff_a.push_back(best_fitting_par[0]);
            temp_punti_max.coeff_b.push_back(best_fitting_par[1]);
            temp_punti_max.coeff_c.push_back(best_fitting_par[2]);
        }
        maxima.push_back(temp_punti_max);
    }
}

//Funzione che genera i massimi con i valori di Fabio generati da root
void get_maxima_root(vector<punti_massimo> parametri, vector<punti_massimo> &maxima)
{
    for (int i = 0; i < parametri.size(); i++)
    {
        punti_massimo temp_punti_max;
        for (int j = 0; j < parametri[i].coeff_a.size(); j++)
        {
            double par_a = parametri[i].coeff_a[j];
            double par_b = parametri[i].coeff_b[j];
            double par_c = parametri[i].coeff_c[j];
            temp_punti_max.t_max.push_back(max_x_fit({par_a, par_b, par_c}));
            temp_punti_max.ampl_max.push_back(2. * M_PI * max_y_fit({par_a, par_b, par_c})); //così da uniformare e non dover più mettere il 2 pigreco dopo
        }
        maxima.push_back(temp_punti_max);
    }
}

//Rapido controllo di andamento di offset per andamento a regime, per genrare grafico
void offset(vector<punti_massimo> vec_punti_max, vector<x_y> &vec_picco_picco)
{
    for (int i = 0; i < vec_punti_max.size(); i++)
    {                                                             //per ogni campione
        vector<double> &ampiezze_max = vec_punti_max[i].ampl_max; //su i punti grezzi
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

//Picco picco a partire da massima verosimiglianza
void picco_picco_mv(vector<punti_massimo> vec_punti_max, vector<x_y> &vec_picco_picco)
{
    for (int i = 0; i < vec_punti_max.size(); i++)
    {                                                             //per ogni campione
        vector<double> &ampiezze_max = vec_punti_max[i].ampl_max; //su i punti grezzi
        x_y temp_picchi;
        for (int j = 0; j < ampiezze_max.size() - 1; j = j + 2) //altirmenti core dumped, evita la correlazione con il +2
        {
            double picco_picco_ = abs(ampiezze_max[j] - ampiezze_max[j + 1]);
            temp_picchi.picco_picco.push_back(picco_picco_ / 2.);
        }
        temp_picchi.freq = vec_punti_max[i].freq_ref;
        vec_picco_picco.push_back(temp_picchi);
    }
}

//Stesso nome perchè il procedimento è lo stesso: leggono gli stessi paramteri
void picco_picco_chi(vector<punti_massimo> vec_punti_max, vector<x_y> &vec_picco_picco)
{
    picco_picco_mv(vec_punti_max, vec_picco_picco);
}

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
        for (int j = 0; j < periodi[i].time.size(); j++) //sono già indipendenti, perchè i tempi venivano presi indipendenti
        {
            temp_omega.push_back(2. * M_PI / periodi[i].time[j]);
        }

        punti_regime[i].omega = media(temp_omega);
        punti_regime[i].err_omega = dstd_media(temp_omega);
    }
}

//Genera grafico istogrmma disperisioni misurazioni di theta in risonanza, usa il valore assoluto di tutti i theta
void gauss_hist_root(vector<punti_massimo> &punti_max, vector<vettoredoppio> &hist, vector<double> bins)
{
    for (int i = 0; i < punti_max.size(); i++)
    {
        vettoredoppio temp_hist;
        double intervals = bins[i];
        vector<double> max_ass;
        vector<double> &massimi = punti_max[i].ampl_max;
        vector<double> counts(intervals, 0);
        for (auto d : massimi)
        {
            if (d > 0)
            {
                max_ass.push_back(d); //
            }
        }
        double max_val = *max_element(max_ass.begin(), max_ass.end());
        double min_val = *min_element(max_ass.begin(), max_ass.end());
        double ampiezza_bin = (max_val - min_val) / intervals;
        vector<double> assex;
        for (int j = 0; j < intervals; j++)
        {
            assex.push_back(min_val + ampiezza_bin * j);
        }

        for (auto c : max_ass)
        {
            for (int j = 0; j < intervals; j++)
            {
                if ((c <= min_val + ampiezza_bin * (j + 1)) && (c > min_val + ampiezza_bin * (j)))
                {
                    counts[j] = counts[j] + 1;
                }
            }
        }
        counts[0] += 1; //per tenere conto del minimo che non viene contato nel if sopra
        for (int k = 0; k < counts.size(); k++)
        {
            temp_hist.vettore2.push_back(assex[k]);
            temp_hist.vettore3.push_back(counts[k] / max_ass.size()); //fa ora le freq relative
            temp_hist.dispersione_amp = dstd(max_ass) / media(max_ass) * 100.;
        }
        hist.push_back(temp_hist);
    }
}

//Genera grafico istogramma per tempi
void gauss_hist_tempi(vector<x_y> &tempi, vector<vettoredoppio> &hist, vector<double> bins)
{
    for (int i = 0; i < tempi.size(); i++)
    {
        vettoredoppio temp_hist;
        double intervals = bins[i];
        vector<double> max_ass;
        vector<double> &massimi = tempi[i].time;
        vector<double> counts(intervals, 0);
        for (auto d : massimi)
        {
            if (d > 0)
            {
                max_ass.push_back(d); //
            }
        }
        double max_val = *max_element(max_ass.begin(), max_ass.end());
        double min_val = *min_element(max_ass.begin(), max_ass.end());
        double ampiezza_bin = (max_val - min_val) / intervals;
        vector<double> assex;
        for (int j = 0; j < intervals; j++)
        {
            assex.push_back(min_val + ampiezza_bin * j);
        }

        for (auto c : max_ass)
        {
            for (int j = 0; j < intervals; j++)
            {
                if ((c <= min_val + ampiezza_bin * (j + 1)) && (c > min_val + ampiezza_bin * (j)))
                {
                    counts[j] = counts[j] + 1;
                }
            }
        }
        counts[0] += 1; //per tenere conto del minimo che non viene contato nel if sopra
        for (int k = 0; k < counts.size(); k++)
        {
            temp_hist.vettore2.push_back(assex[k]);
            temp_hist.vettore3.push_back(counts[k] / max_ass.size()); //fa ora le freq relative
            temp_hist.dispersione_amp = dstd(max_ass) / media(max_ass) * 100.;
        }
        hist.push_back(temp_hist);
    }
}

//Calcolo di omega per smorzamento
void omega_s_scamorza(vector<x_y> periodi, vector<x_y> &omega_s_smorzamento)
{
    for (int i = 0; i < periodi.size(); i++)
    {
        x_y temp_omega_str;
        vector<double> temp_omega;
        for (int j = 0; j < periodi[i].time.size(); j++) //sono già indipendenti, perchè i tempi venivano presi indipendenti
        {
            temp_omega.push_back(2. * M_PI / periodi[i].time[j]);
        }

        temp_omega_str.omega_media2 = media(temp_omega);
        temp_omega_str.err_omega_media2 = dstd_media(temp_omega);
        omega_s_smorzamento.push_back(temp_omega_str);
    }
}

//Calcolo di compatibiità fra omega sperimentale e omega regime per fase di risonanza
void compatibilita_omega(vector<punto_regime> &camp_lorent, vector<data_campione> freq_camp, vector<double> &compatib_omega)
{
    for (int i = 0; i < freq_camp.size(); i++)
    {
        double omega_th = freq_camp[i].data_file_freq * 2. * M_PI / 1000.; //convertita in Hz da mHz e in omega da freq
        double omega_sper = camp_lorent[i].omega;
        double err_omega_sper = camp_lorent[i].err_omega;
        double err_omega_th = sigma_dist_uni(0.001, 1); //dist uniforme su più piccola cifra degli Hz
        compatib_omega.push_back(comp_3(omega_th, omega_sper, err_omega_sper, err_omega_th));
    }
}

double y_th(double x, vector<double> parms)
{
    double a, b, c, d;
    a = parms[0];
    b = parms[1];
    c = parms[2];
    d = parms[3];
    return a / sqrt(pow((pow(b, 2) + 2 * pow(c, 2) - pow(x, 2)), 2) + 4 * pow((c * x), 2)) + d;
}

double post_lor(vector<punto_regime> &campana_lor, vector<double> parametri_fit_gnuplot)
{
    double sum_scarto_quad = 0;
    double gdl = campana_lor.size() - parametri_fit_gnuplot.size();
    for (int i = 0; i < campana_lor.size(); i++) //per tutti i punti della lorenziana
    {
        //y_i-y_i,th
        double scarto = campana_lor[i].theta - y_th(campana_lor[i].omega, parametri_fit_gnuplot);
        sum_scarto_quad += pow(scarto, 2);
    }
    return sqrt(sum_scarto_quad / gdl);
}

//Funzione checalcola la linearizzazione sui punti di massimo
void linearize_max(vector<punti_massimo> &theta_generiche, vector<punti_massimo> &ln_theta, double err_post)
{
    for (auto d : theta_generiche)
    {
        punti_massimo temp_ln_maxima;
        for (int j = 0; j < d.ampl_max.size(); j++)
        {
            if (d.ampl_max[j] > 0)
            {
                temp_ln_maxima.t_max.push_back(d.t_max[j]);
                temp_ln_maxima.ampl_max.push_back(log(d.ampl_max[j]));
                temp_ln_maxima.err_ampl_max.push_back(sqrt(pow(1. / d.ampl_max[j] * err_post, 2)));
            }
        }

        ln_theta.push_back(temp_ln_maxima);
    }
}

//Funzione checalcola la linearizzazione sui punti di minimo
void linearize_min(vector<punti_massimo> &theta_generiche, vector<punti_massimo> &ln_theta, double err_post)
{
    for (auto d : theta_generiche)
    {
        punti_massimo temp_ln_maxima;
        for (int j = 0; j < d.ampl_max.size(); j++)
        {
            if (d.ampl_max[j] < 0)
            {
                temp_ln_maxima.t_max.push_back(d.t_max[j]);
                temp_ln_maxima.ampl_max.push_back(-log(-d.ampl_max[j])); //usa il meno com vole la cinzia ;)
                temp_ln_maxima.err_ampl_max.push_back(sqrt(pow(1. / d.ampl_max[j] * err_post, 2)));
            }
        }
        ln_theta.push_back(temp_ln_maxima);
    }
}

//Calcola i parametri del fit per gamma dopo aver fatto il log naturale, usa err a posteriori derivanti da lorenz
void return_angolari(vector<punti_massimo> &ln_theta_max, vector<punti_massimo> &ln_theta_min, vector<interpolazione_gamma> &parametri_intepolazioni)
{
    for (int i = 0; i < ln_theta_max.size(); i++)
    {
        interpolazione_gamma temp_gamma;
        //per i massimi
        temp_gamma.a_intercetta_gamma_max = a_intercetta(ln_theta_max[i].t_max, ln_theta_max[i].ampl_max, ln_theta_max[i].err_ampl_max);
        temp_gamma.b_angolare_gamma_max = b_angolare(ln_theta_max[i].t_max, ln_theta_max[i].ampl_max, ln_theta_max[i].err_ampl_max);
        temp_gamma.r_max = pearson(ln_theta_max[i].t_max, ln_theta_max[i].ampl_max);
        temp_gamma.t_max = student(ln_theta_max[i].t_max, ln_theta_max[i].ampl_max);
        temp_gamma.err_a_post_max = sigma_a_posteriori(ln_theta_max[i].t_max, ln_theta_max[i].ampl_max);
        temp_gamma.err_b_post_max = sigma_b_posteriori(ln_theta_max[i].t_max, ln_theta_max[i].ampl_max);

        //per i minimi
        temp_gamma.a_intercetta_gamma_min = a_intercetta(ln_theta_min[i].t_max, ln_theta_min[i].ampl_max, ln_theta_min[i].err_ampl_max);
        temp_gamma.b_angolare_gamma_min = b_angolare(ln_theta_min[i].t_max, ln_theta_min[i].ampl_max, ln_theta_min[i].err_ampl_max);
        temp_gamma.r_min = pearson(ln_theta_min[i].t_max, ln_theta_min[i].ampl_max);
        temp_gamma.t_min = student(ln_theta_min[i].t_max, ln_theta_min[i].ampl_max);
        temp_gamma.err_a_post_min = sigma_a_posteriori(ln_theta_min[i].t_max, ln_theta_min[i].ampl_max);
        temp_gamma.err_b_post_min = sigma_b_posteriori(ln_theta_min[i].t_max, ln_theta_min[i].ampl_max);

        //salva tutto in vec
        parametri_intepolazioni.push_back(temp_gamma);
    }
}
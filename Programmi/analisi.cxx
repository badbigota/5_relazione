#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <algorithm>
#include "functions.h"
using namespace std;

/*
void funzione(lettura, scrittura)
*/

int main()
{
    vector<data_campione> campione;
    vector<x_y> tempi;
    vector<x_y> periodi;
    vector<x_y> media_periodi;
    vector<x_y> omega1;
    vector<x_y> omega2;

    string map_file = "../Dati/mappa.txt";
    load_data(map_file, campione);
    get_zero_time(campione, tempi);           //calcolo di tutti tempi con interpolazione a 4 punti
    get_periods(tempi, periodi);              //calcolo di tutti i periodi
    get_periodi_medi(periodi, media_periodi); //calcola media dei periodi con errore

    //Stampa file per verifica corrispondenza punti interpolazione con effettiva intercetta con asse x
    //int num_file = 0;
    //for (int j = 0;  j < tempi.size()-1; j++)
    //{
    //	ofstream fout_time("../Dati/zeroes_"+to_string(num_file)+".csv");
    //	for (int i = 0; i < tempi[j].time.size(); i++)
    //	{
    //		fout_time<< tempi[j].time[i]<< "\t"<<tempi[j].amplitude[i]<<"\t0" << endl;
    //	}
    //    num_file++;
    //}

    //Stampa tempi
    //for(int i=0;i<tempi.size();i++)
    //{
    //    cout<<endl<<"Pendolo "<<i+1<<endl;
    //    for(int j=0;j<tempi[i].time.size();j++)
    //    {
    //        cout<<tempi[i].time[j]<<",  ";
    //    }
    //}

    //Stampa per prova dei periodi
    //for (int k = 0; k < periodi.size(); k++)
    //{
    //    ofstream fout("../Test/periodi_" + to_string(k + 1) + ".txt");
    //    for (int i = 0; i < periodi[k].time.size(); i++)
    //    {
    //        fout << periodi[k].time[i] << endl;
    //    }
    //}

    //Stampa media periodi
    //for (auto d : media_periodi)
    //{
    //    cout << d.avg_time << "+/-" << d.err_avg_time << endl;
    //}

    //OMEGA SPERIMENTALE
    //Primo metodo omega, calcola tutti gli omega e poi fa errore su di essi, prediletto da JCGM
    for (int i = 0; i < periodi.size(); i++)
    {
        x_y temp_omega;
        for (auto d : periodi[i].time)
        {
            temp_omega.omega_sperim.push_back(2 * M_PI / d);
        }
        omega1.push_back(temp_omega);
        omega1[i].omega_media1 = media(omega1[i].omega_sperim);
        omega1[i].err_omega_media1 = dstd_media(omega1[i].omega_sperim);
    }

    //Secondo metodo omega, con la media dei periodi
    for (int i = 0; i < media_periodi.size(); i++)
    {
        x_y temp_omega2;
        temp_omega2.omega_media2 = 2 * M_PI / media_periodi[i].avg_time;
        temp_omega2.err_omega_media2 = 2 * M_PI * media_periodi[i].err_avg_time / pow(media_periodi[i].avg_time, 2);
        omega2.push_back(temp_omega2);
    }

    //cout << "OmegaMedia1"
    //     << "\t"
    //     << "ErrOmegaMedia1"
    //     << "\t"
    //     << "OmegaMedia2"
    //     << "\t"
    //     << "ErrOmegaMedia2" << endl;
    //for (int i = 0; i < omega1.size(); i++)
    //{
    //    cout << omega1[i].omega_media1 << "\t" << omega1[i].err_omega_media1 << "\t" << omega2[i].omega_media2 << "\t" << omega2[i].err_omega_media2 << endl;
    //}

    //for (int i = 0; i < campione.size(); i++)
    //{
    //    ofstream lout("../Dati/Range" + to_string(i + 1) + ".txt");
    //    lout << endl
    //         << "Pendolo " << i + 1 << endl;
    //    for (int j = 0; j < pos_max[i].vettore2.size(); j++)
    //    {
    //        lout << "\t" << campione[i].time[pos_max[i].vettore2[j] - 9] << "\t"
    //             << "," << campione[i].time[pos_max[i].vettore2[j] + 9] << "," << endl;
    //    }
    //    cout << endl;
    //}

    //Genrea vettore con le posizioni dei massimi
    vector<vettoredoppio> indici_dei_massimi;
    get_index_maxima(campione, tempi, indici_dei_massimi);

    //Genera vettore di strutture con
    vector<punti_massimo> punti_di_massimo;
    get_maxima(campione, indici_dei_massimi, punti_di_massimo);
    cout << punti_di_massimo[0].t_max.size() << endl;

    return 0;
}
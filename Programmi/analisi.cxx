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
    vector<data_campione> campione_decadimento;
    vector<x_y> tempi;
    vector<x_y> periodi;
    vector<x_y> media_periodi;
    vector<x_y> omega1;
    vector<x_y> omega2;
    vector<punti_massimo> parametri;

    string map_file = "../Dati/mappa.txt";
    load_data(map_file, campione);
    load_data_decay(campione_decadimento);
    get_zero_time(campione, tempi);           //calcolo di tutti tempi con interpolazione a 4 punti
    get_periods(tempi, periodi);              //calcolo di tutti i periodi
    get_periodi_medi(periodi, media_periodi); //calcola media dei periodi con errore

    //Stampa file per verifica corrispondenza punti interpolazione con effettiva intercetta con asse x
    //int num_file = 0;
    //for (int j = 0; j < tempi.size() - 1; j++)
    //{
    //    ofstream fout_time("../Dati/zeroes_" + to_string(num_file) + ".csv");
    //    for (int i = 0; i < tempi[j].time.size(); i++)
    //    {
    //        fout_time << tempi[j].time[i] << "\t" << tempi[j].amplitude[i] << "\t0" << endl;
    //    }
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
    //for (int i = 0; i < periodi.size(); i++)
    //{
    //    x_y temp_omega;
    //    for (auto d : periodi[i].time)
    //    {
    //        temp_omega.omega_sperim.push_back(2 * M_PI / d);
    //    }
    //    omega1.push_back(temp_omega);
    //    omega1[i].omega_media1 = media(omega1[i].omega_sperim);
    //    omega1[i].err_omega_media1 = dstd_media(omega1[i].omega_sperim);
    //}

    //Secondo metodo omega, con la media dei periodi
    //for (int i = 0; i < media_periodi.size(); i++)
    //{
    //    x_y temp_omega2;
    //    temp_omega2.omega_media2 = 2 * M_PI / media_periodi[i].avg_time;
    //    temp_omega2.err_omega_media2 = 2 * M_PI * media_periodi[i].err_avg_time / pow(media_periodi[i].avg_time, 2);
    //    omega2.push_back(temp_omega2);
    //}

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

    //Genrea vettore con le posizioni dei massimi
    vector<vettoredoppio> indici_dei_massimi;
    get_index_maxima(campione, tempi, indici_dei_massimi);
    //for (int j = 0; j < indici_dei_massimi.size(); j++)
    //{
    //    ofstream fout_maxima("../Dati/maxima_" + to_string(j) + ".txt");
    //    for (auto c : indici_dei_massimi[j].vettore2)
    //    {
    //        fout_maxima << campione[j].time[c] << "\t" << campione[j].a[c] << endl;
    //    }
    //}

    //for (int i = 0; i < campione.size(); i++)
    //{
    //    ofstream lout("../Dati/Range_" + to_string(i) + ".txt");
    //    for (int j = 0; j < indici_dei_massimi[i].vettore2.size(); j++)
    //    {
    //        lout << "," << campione[i].time[indici_dei_massimi[i].vettore2[j] - 9] << ","
    //             << " \t " << campione[i].time[indici_dei_massimi[i].vettore2[j] + 9] << endl;
    //    }
    //}
    //Genera vettore di strutture con
    vector<punti_massimo> punti_di_massimo;      //usati sia per mv che per assoluta
    vector<punti_massimo> punti_di_massimo_root; //solo per root
    get_maxima(campione, indici_dei_massimi, punti_di_massimo);
    lettura_parametri(parametri);
    //cout << parametri[0].coeff_c.size();

    get_maxima_root(parametri, punti_di_massimo_root);

    //for (int i = 0; i < punti_di_massimo.size(); i++)
    //{
    //    ofstream fout_temp("../Dati/punti_max_" + to_string(punti_di_massimo[i].freq_ref) + ".txt");
    //    for (int j = 0; j < punti_di_massimo[i].t_max.size(); j++)
    //    {
    //        fout_temp << punti_di_massimo[i].t_max[j] << "\t" << punti_di_massimo[i].ampl_max[j] << "\t" << punti_di_massimo[i].coeff_a[j] << "\t" << punti_di_massimo[i].coeff_b[j] << "\t" << punti_di_massimo[i].coeff_c[j] << endl;
    //    }
    //}

    vector<x_y> valutaz_offset;
    offset(punti_di_massimo, valutaz_offset);
    //for (auto d : valutaz_offset)
    //{
    //    ofstream fout_off("../Dati/offset_" + to_string(d.freq) + ".txt");
    //    for (auto c : d.offset)
    //    {
    //        fout_off << c << endl;
    //    }
    //}

    vector<x_y> picchi_picchi_assoluti;
    vector<x_y> picchi_picchi_mv;
    vector<x_y> picchi_picchi_root;

    picco_picco_ass(punti_di_massimo, picchi_picchi_assoluti);
    picco_picco_mv(punti_di_massimo, picchi_picchi_mv);
    picco_picco_chi(punti_di_massimo_root, picchi_picchi_root);

    //cout << picchi_picchi_assoluti.size() << endl;
    //for (auto d : picchi_picchi_assoluti)
    //{
    //    ofstream fout_picco("../Dati/pic_pic_mezz_ass_" + to_string(d.freq) + ".txt");
    //    for (auto c : d.picco_picco)
    //    {
    //        fout_picco << c << endl;
    //    }
    //}

    vector<punto_regime> campana_lorentziana_assoluta;
    vector<punto_regime> campana_lorentziana_mv;
    vector<punto_regime> campana_lorentziana_root;

    add_picco_medio(picchi_picchi_assoluti, campana_lorentziana_assoluta);
    add_picco_medio(picchi_picchi_mv, campana_lorentziana_mv);
    add_picco_medio(picchi_picchi_root, campana_lorentziana_root);

    add_omega_2(periodi, campana_lorentziana_assoluta);
    add_omega_2(periodi, campana_lorentziana_mv);
    add_omega_2(periodi, campana_lorentziana_root);

    //STAMPA LA LORENTZIANA

    ofstream fout_lor_mv("lore_mv.txt");
    ofstream fout_lor_ass("lore_ass.txt");
    ofstream fout_lor_root("lore_root.txt");

    for (int i = 0; i < campana_lorentziana_mv.size(); i++)
    {
        fout_lor_mv << campana_lorentziana_mv[i].omega << "\t" << campana_lorentziana_mv[i].theta << "\t" << campana_lorentziana_mv[i].err_omega << "\t" << campana_lorentziana_mv[i].err_theta << endl;
    }

    for (int i = 0; i < campana_lorentziana_assoluta.size(); i++)
    {
        fout_lor_ass << campana_lorentziana_assoluta[i].omega << "\t" << campana_lorentziana_assoluta[i].theta << "\t" << campana_lorentziana_assoluta[i].err_omega << "\t" << campana_lorentziana_assoluta[i].err_theta << endl;
    }

    for (int i = 0; i < campana_lorentziana_root.size(); i++)
    {
        fout_lor_root << campana_lorentziana_root[i].omega << "\t" << campana_lorentziana_root[i].theta * 2. * M_PI << "\t" << campana_lorentziana_root[i].err_omega << "\t" << campana_lorentziana_root[i].err_theta << endl;
    }

    return 0;
}
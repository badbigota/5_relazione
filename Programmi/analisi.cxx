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
    vector<x_y> tempi_decay;
    vector<x_y> periodi;
    vector<x_y> periodi_decadimento;
    vector<x_y> media_periodi;
    vector<x_y> media_periodi_decadimento;
    vector<x_y> omega1;
    vector<x_y> omega2;
    vector<punti_massimo> parametri;
    vector<x_y> omega_s_smorzamento;

    string map_file = "../Dati/mappa.txt";
    load_data(map_file, campione);
    get_zero_time(campione, tempi);           //calcolo di tutti tempi con interpolazione a 4 punti
    get_periods(tempi, periodi);              //calcolo di tutti i periodi
    get_periodi_medi(periodi, media_periodi); //calcola media dei periodi con errore

    //Stampa file per verifica corrispondenza punti interpolazione con effettiva intercetta con asse x
    int num_file = 0;
    for (int j = 0; j < tempi.size(); j++)
    {
        ofstream fout_time("../Dati/Zeros_Regime/zeroes_" + to_string(num_file) + ".csv");
        for (int i = 0; i < tempi[j].time.size(); i++)
        {
            fout_time << tempi[j].time[i] << "\t" << tempi[j].amplitude[i] << "\t0" << endl;
        }
        num_file++;
    }

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
    for (int k = 0; k < periodi.size(); k++)
    {
        ofstream fout("../Dati/Periodi/periodi_" + to_string(k + 1) + ".txt");
        for (int i = 0; i < periodi[k].time.size(); i++)
        {
            fout << periodi[k].time[i] << endl;
        }
    }

    //Stampa media periodi
    //for (auto d : media_periodi)
    //{
    //    cout << d.avg_time << "+/-" << d.err_avg_time << endl;
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


    for (int j = 0; j < indici_dei_massimi.size(); j++)
    {
        ofstream fout_maxima("../Dati/Maxima/maxima_" + to_string(j) + ".txt");
        for (auto c : indici_dei_massimi[j].vettore2)
        {
            fout_maxima << campione[j].time[c] << "\t" << campione[j].a[c] << endl;
        }
    }

    for (int i = 0; i < campione.size(); i++)
    {
        ofstream lout("../Dati/Range_" + to_string(i) + ".txt");
        for (int j = 0; j < indici_dei_massimi[i].vettore2.size(); j++)
        {
            lout << "," << campione[i].time[indici_dei_massimi[i].vettore2[j] - 9] << ","
                 << " \t " << campione[i].time[indici_dei_massimi[i].vettore2[j] + 9] << endl;
        }
    }
    //Genera vettore di strutture con
    vector<punti_massimo> punti_di_massimo_mv;
    vector<punti_massimo> punti_di_massimo_ass;  //usati sia per mv che per assoluta
    vector<punti_massimo> punti_di_massimo_root; //solo per root
    get_maxima_mv(campione, indici_dei_massimi, punti_di_massimo_mv);
    get_maxima_ass(campione, indici_dei_massimi, punti_di_massimo_ass);

    lettura_parametri(parametri);
    get_maxima_root(parametri, punti_di_massimo_root);
    //cout << parametri[0].coeff_c.size();

    //for (int i = 0; i < punti_di_massimo_root.size(); i++)
    //{
    //    ofstream fout_temp("../Dati/Punti_max/punti_max_" + to_string(punti_di_massimo_root[i].freq_ref) + ".txt");
    //    for (int j = 0; j < punti_di_massimo_root[i].t_max.size(); j++)
    //    {
    //        fout_temp << punti_di_massimo_root[i].t_max[j] << "\t" << punti_di_massimo_root[i].ampl_max[j] << "\t" << punti_di_massimo_root[i].coeff_a[j] << "\t" << punti_di_massimo_root[i].coeff_b[j] << "\t" << punti_di_massimo_root[i].coeff_c[j] << endl;
    //    }
    //}

    vector<x_y> valutaz_offset_ass;
    offset(punti_di_massimo_ass, valutaz_offset_ass);
    for (auto d : valutaz_offset_ass)
    {
        ofstream fout_off("../Dati/Offset_regime/offset_" + to_string(d.freq) + ".txt");
        for (auto c : d.offset)
        {
            fout_off << c << endl;
        }
    }

    vector<x_y> picchi_picchi_assoluti;
    vector<x_y> picchi_picchi_mv;
    vector<x_y> picchi_picchi_root;

    picco_picco_ass(punti_di_massimo_ass, picchi_picchi_assoluti);
    picco_picco_mv(punti_di_massimo_mv, picchi_picchi_mv);
    picco_picco_chi(punti_di_massimo_root, picchi_picchi_root);

    //cout << picchi_picchi_assoluti.size() << endl;
    for (auto d : picchi_picchi_assoluti)
    {
        ofstream fout_picco("../Dati/Picchi_ass/pic_pic_mezz_ass_" + to_string(d.freq) + ".txt");
        for (auto c : d.picco_picco)
        {
            fout_picco << c << endl;
        }
    }

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
        fout_lor_mv << campana_lorentziana_mv[i].omega << "\t" << campana_lorentziana_mv[i].theta << "\t" << campana_lorentziana_mv[i].err_omega << "\t" << campana_lorentziana_mv[i].err_theta << "\t 1" << endl;
    }

    for (int i = 0; i < campana_lorentziana_assoluta.size(); i++)
    {
        fout_lor_ass << campana_lorentziana_assoluta[i].omega << "\t" << campana_lorentziana_assoluta[i].theta << "\t" << campana_lorentziana_assoluta[i].err_omega << "\t" << campana_lorentziana_assoluta[i].err_theta << "\t 1" << endl;
    }

    for (int i = 0; i < campana_lorentziana_root.size(); i++)
    {
        fout_lor_root << campana_lorentziana_root[i].omega << "\t" << campana_lorentziana_root[i].theta << "\t" << campana_lorentziana_root[i].err_omega << "\t" << campana_lorentziana_root[i].err_theta << "\t 1" << endl;
    }

    vector<vettoredoppio> hist_root;
    vector<double> n_bin(20, 10.0); //definisce il numero di bin, cambia il secondo valore
    gauss_hist_root(punti_di_massimo_root, hist_root, n_bin);
    ofstream fout_disp("../Histogram/dispersione_sigma.txt");
    for (int j = 0; j < hist_root.size(); j++)
    {
        ofstream fout_hist("../Histogram/" + to_string(j + 10) + "-hist.txt");
        for (int k = 0; k < hist_root[j].vettore2.size(); k++)
        {
            fout_hist << hist_root[j].vettore2[k] << "\t" << hist_root[j].vettore3[k] << endl;
        }
        fout_disp << j << "\t" << hist_root[j].dispersione_amp << endl;
    }

    vector<double> compatibilty_ass;
    vector<double> compatibilty_mv;
    vector<double> compatibilty_root;
    compatibilita_omega(campana_lorentziana_assoluta, campione, compatibilty_ass);
    compatibilita_omega(campana_lorentziana_assoluta, campione, compatibilty_mv);
    compatibilita_omega(campana_lorentziana_assoluta, campione, compatibilty_root);

    //cout << "Compatiblità omega sper root con omega th";

    ofstream fout_com("../Dati/Compatib/comp.txt");
    for (int i = 0; i < campana_lorentziana_root.size(); i++) //sono tutte uguali, non cambia un cazzo con le altre
    {
        fout_com << i + 1 << "\t" << campana_lorentziana_root[i].omega << "\t" << campana_lorentziana_root[i].err_omega << "\t" << campione[i].data_file_freq * 2. * M_PI / 1000. << "\t" << sigma_dist_uni(0.001, 1) << "\t" << compatibilty_root[i] << endl;
    }

    //*******************************************************SMORZAMENTO*****************************************************************************************************
    load_data_decay(campione_decadimento);
    get_zero_time(campione_decadimento, tempi_decay);

    int num_file_d = 0;
    for (int p = 0; p < tempi_decay.size(); p++)
    {
        ofstream fout_time_smorzamento("../Dati/Zeros_smorzamento/zeroes_d_" + to_string(num_file_d) + ".csv");
        for (int i = 0; i < tempi_decay[p].time.size(); i++)
        {
            fout_time_smorzamento << tempi_decay[p].time[i] << "\t 0" << endl;
        }
        num_file_d++;
    }

    get_periods(tempi_decay, periodi_decadimento);
    get_periodi_medi(periodi_decadimento, media_periodi_decadimento); //#cavallo

    vector<vettoredoppio> indici_dei_massimi_decadimento;
    vector<punti_massimo> massimi_decadimento_ass;
    vector<punti_massimo> massimi_decadimento_mv;
    get_index_maxima(campione_decadimento, tempi_decay, indici_dei_massimi_decadimento);
    omega_s_scamorza(periodi_decadimento, omega_s_smorzamento);
    get_maxima_ass(campione_decadimento, indici_dei_massimi_decadimento, massimi_decadimento_ass);
    get_maxima_mv(campione_decadimento, indici_dei_massimi_decadimento, massimi_decadimento_mv);

    vector<punti_massimo> parametri_smorzamento_root;
    vector<punti_massimo> massimi_decadimento_root;
    lettura_parametri_smorz(parametri_smorzamento_root);
    get_maxima_root(parametri_smorzamento_root, massimi_decadimento_root);

    //Stampa i punti max di smorzamento
    for (int i = 0; i < massimi_decadimento_ass.size(); i++)
    {
        ofstream fout_smor_max("../Dati/Maxmin_smorzamento/smorz_" + to_string(i) + "_ass.txt");
        for (int j = 0; j < massimi_decadimento_ass[i].t_max.size(); j++)
        {
            fout_smor_max << massimi_decadimento_ass[i].t_max[j] << "\t" << massimi_decadimento_ass[i].ampl_max[j] << endl;
        }
    }

    for (int i = 0; i < massimi_decadimento_root.size(); i++)
    {
        ofstream fout_smor_max("../Dati/Maxmin_smorzamento/smorz_" + to_string(i) + "_root.txt");
        for (int j = 0; j < massimi_decadimento_root[i].t_max.size(); j++)
        {
            fout_smor_max << massimi_decadimento_root[i].t_max[j] << "\t" << massimi_decadimento_root[i].ampl_max[j] << endl;
        }
    }
    for (int i = 0; i < massimi_decadimento_mv.size(); i++)
    {
        ofstream fout_smor_max("../Dati/Maxmin_smorzamento/smorz_" + to_string(i) + "_mv.txt");
        for (int j = 0; j < massimi_decadimento_mv[i].t_max.size(); j++)
        {
            fout_smor_max << massimi_decadimento_mv[i].t_max[j] << "\t" << massimi_decadimento_mv[i].ampl_max[j] << endl;
        }
    }

    //Stampa/controllo righe precendenti
    /*int h = 1;
    for (int j = 0; j < periodi_decadimento.size(); j++)
    {
        cout << h << endl;
        h++;
        for (int i = 0; i < periodi_decadimento[j].time.size(); i++)
        {
            cout << periodi_decadimento[j].time[i] << endl;
        }
        cout << "Periodo medio: " << media_periodi_decadimento[j].avg_time << "+/-" << media_periodi_decadimento[j].err_avg_time << endl;
    }*/

    //TUTTA ROBA CHE SERVE A FABIO - intervalli per max root
    for (int j = 0; j < indici_dei_massimi_decadimento.size(); j++)
    {
        ofstream smorz("../Dati/periodo_smorzamento" + to_string(j) + ".txt");
    
        for (int i = 0; i < indici_dei_massimi_decadimento[j].vettore2.size(); i++)
        {
            smorz << campione_decadimento[j].time[indici_dei_massimi_decadimento[j].vettore2[i] - 9] << ",\t" << campione_decadimento[j].time[indici_dei_massimi_decadimento[j].vettore2[i] + 9] << ",\t" << endl;
        }
        smorz << endl;
    }

    /*int h = 1;
    for (auto c : omega_s_smorzamento)
    {
        cout << h << endl;
        h++;
        cout << c.omega_media2 << "+/-" << c.err_omega_media2 << endl;
    }*/

    vector<punti_massimo> punti_max_ln_root;
    vector<punti_massimo> punti_min_ln_root;
    vector<punti_massimo> punti_max_ln_ass;
    vector<punti_massimo> punti_min_ln_ass;
    vector<punti_massimo> punti_max_ln_mv;
    vector<punti_massimo> punti_min_ln_mv;

    //Prende i valori del fit di gnuplot con err a posteriori
    vector<double> par_root = {1.82653, 6.06287, -0.0558678, -0.026352};
    vector<double> par_ass = {1.84435, 6.0624, -0.0566001, 0.002301};
    vector<double> par_mv = {1.82795, 6.06275, -0.0559494, -0.0265116};
    double err_post_root = post_lor(campana_lorentziana_root, par_root);
    double err_post_ass = post_lor(campana_lorentziana_assoluta, par_ass);
    double err_post_mv = post_lor(campana_lorentziana_mv, par_mv);

    //Linearizza trovando i massimi/minimi pronti per essere plottati con gli err a posteriori precedentemente calclati e poi propagati
    linearize_max(massimi_decadimento_root, punti_max_ln_root, err_post_root);
    linearize_min(massimi_decadimento_root, punti_min_ln_root, err_post_root);
    linearize_max(massimi_decadimento_ass, punti_max_ln_ass, err_post_ass);
    linearize_min(massimi_decadimento_ass, punti_min_ln_ass, err_post_ass);
    linearize_max(massimi_decadimento_mv, punti_max_ln_mv, err_post_mv);
    linearize_min(massimi_decadimento_mv, punti_min_ln_mv, err_post_mv);

    //Stampa tutti i ln(theta) per fare i grafici con relativi errori
    for (int i = 0; i < punti_max_ln_root.size(); i++)
    {
        ofstream fout_ln_max("../Dati/Linearize/ln_max_" + to_string(i) + "_root.txt");
        ofstream fout_ln_min("../Dati/Linearize/ln_min_" + to_string(i) + "_root.txt");
        for (int j = 0; j < punti_max_ln_root[i].ampl_max.size(); j++)
        {
            fout_ln_max << punti_max_ln_root[i].t_max[j] << "\t" << punti_max_ln_root[i].ampl_max[j] << "\t" << punti_max_ln_root[i].err_ampl_max[j] << endl;
            fout_ln_min << punti_min_ln_root[i].t_max[j] << "\t" << punti_min_ln_root[i].ampl_max[j] << "\t" << punti_min_ln_root[i].err_ampl_max[j] << endl;
        }
    }
    for (int i = 0; i < punti_max_ln_ass.size(); i++)
    {
        ofstream fout_ln_max("../Dati/Linearize/ln_max_" + to_string(i) + "_ass.txt");
        ofstream fout_ln_min("../Dati/Linearize/ln_min_" + to_string(i) + "_ass.txt");
        for (int j = 0; j < punti_max_ln_ass[i].ampl_max.size(); j++)
        {
            fout_ln_max << punti_max_ln_ass[i].t_max[j] << "\t" << punti_max_ln_ass[i].ampl_max[j] << "\t" << punti_max_ln_ass[i].err_ampl_max[j] << endl;
            fout_ln_min << punti_min_ln_ass[i].t_max[j] << "\t" << punti_min_ln_ass[i].ampl_max[j] << "\t" << punti_min_ln_ass[i].err_ampl_max[j] << endl;
        }
    }
    for (int i = 0; i < punti_max_ln_mv.size(); i++)
    {
        ofstream fout_ln_max("../Dati/Linearize/ln_max_" + to_string(i) + "_mv.txt");
        ofstream fout_ln_min("../Dati/Linearize/ln_min_" + to_string(i) + "_mv.txt");
        for (int j = 0; j < punti_max_ln_mv[i].ampl_max.size(); j++)
        {
            fout_ln_max << punti_max_ln_mv[i].t_max[j] << "\t" << punti_max_ln_mv[i].ampl_max[j] << "\t" << punti_max_ln_mv[i].err_ampl_max[j] << endl;
            fout_ln_min << punti_min_ln_mv[i].t_max[j] << "\t" << punti_min_ln_mv[i].ampl_max[j] << "\t" << punti_min_ln_mv[i].err_ampl_max[j] << endl;
        }
    }

    vector<interpolazione_gamma> parametri_interpolazioni_root;
    vector<interpolazione_gamma> parametri_interpolazioni_ass;
    vector<interpolazione_gamma> parametri_interpolazioni_mv;
    return_angolari(punti_max_ln_root, punti_min_ln_root, parametri_interpolazioni_root);
    return_angolari(punti_max_ln_ass, punti_min_ln_ass, parametri_interpolazioni_ass);
    return_angolari(punti_max_ln_mv, punti_min_ln_mv, parametri_interpolazioni_mv);

    // for (int i = 0; i < parametri_interpolazioni.size(); i++)
    // {
    //     cout<<"Massimi "<<i+1<<endl;
    //     cout << parametri_interpolazioni[i].a_intercetta_gamma_max<<"\t"<< parametri_interpolazioni[i].b_angolare_gamma_max<<endl;
    //     cout<< "Minimi "<<i+1 <<endl;
    //     cout << parametri_interpolazioni[i].a_intercetta_gamma_min<<"\t"<< parametri_interpolazioni[i].b_angolare_gamma_min<<endl;
    // }

    return 0;
}
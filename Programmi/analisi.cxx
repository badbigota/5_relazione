#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <algorithm>
#include "functions.h"
using namespace std;

int main()
{
    vector<data_campione> campione;
    vector<data_campione> campione_decadimento;
    vector<x_y> tempi;
    vector<x_y> tempi_forcing;
    vector<x_y> tempi_decay;
    vector<x_y> periodi;
    vector<x_y> periodi_forcing;
    vector<x_y> periodi_decadimento;
    vector<x_y> periodi_decadimento_reiettati;
    vector<x_y> media_periodi;
    vector<x_y> media_periodi_forcing;
    vector<x_y> media_periodi_decadimento;
    vector<x_y> media_periodi_decadimento_reiettati;
    vector<x_y> omega1;
    vector<x_y> omega2;
    vector<punti_massimo> parametri;
    vector<x_y> omega_s_smorzamento;

    string map_file = "../Dati/mappa.txt";
    load_data(map_file, campione);

    /*Calcolo dei periodi con ampiezza e forcing*/
    get_zero_time(campione, tempi);           //calcolo di tutti tempi con interpolazione a 4 punti
    get_periods(tempi, periodi);              //calcolo di tutti i periodi
    get_periodi_medi(periodi, media_periodi); //calcola media dei periodi con errore
    get_zero_time_forcing(campione, tempi_forcing);
    get_periods(tempi_forcing, periodi_forcing);
    get_periodi_medi(periodi_forcing, media_periodi_forcing);

    //Stampa file per verifica corrispondenza punti interpolazione con effettiva intercetta con asse x
    int num_file = 0;
    for (int j = 0; j < tempi.size(); j++)
    {
        ofstream fout_time("../Dati/Zeros_Regime/z_" + to_string(num_file) + "_amp.txt");
        for (int i = 0; i < tempi[j].time.size(); i++)
        {
            fout_time << tempi[j].time[i] << "\t" << tempi[j].amplitude[i] << "\t0" << endl;
        }
        num_file++;
    }

    int num_file_f = 0;
    for (int j = 0; j < tempi_forcing.size(); j++)
    {
        ofstream fout_time_f("../Dati/Zeros_Regime/z_" + to_string(num_file_f) + "_forcing.txt");
        for (int i = 0; i < tempi_forcing[j].time.size(); i++)
        {
            fout_time_f << tempi_forcing[j].time[i] << "\t0" << endl;
        }
        num_file_f++;
    }

    //Stampa per prova dei periodi sia da dati puri che da forcing
    for (int k = 0; k < periodi.size(); k++)
    {
        ofstream fout_a("../Dati/Periodi/periodi_" + to_string(k) + "_amp.txt");
        for (int i = 0; i < periodi[k].time.size(); i++)
        {
            fout_a << periodi[k].time[i] << endl;
        }
    }
    for (int k = 0; k < periodi_forcing.size(); k++)
    {
        ofstream fout_f("../Dati/Periodi/periodi_" + to_string(k) + "_for.txt");
        for (int i = 0; i < periodi_forcing[k].time.size(); i++)
        {
            fout_f << periodi_forcing[k].time[i] << endl;
        }
    }

    //Stampa media periodi con errore, sia di grezzi che da forzante
    ofstream fout_media_tempi("../Dati/Periodi/periodi_medi.txt");
    ofstream fout_media_tempi_f("../Dati/Periodi/periodi_medi_forzante.txt");
    for (int i = 0; i < media_periodi.size(); i++)
    {
        fout_media_tempi << i + 1 << "\t" << media_periodi[i].avg_time << "\t" << media_periodi[i].err_avg_time << "\t" << media_periodi[i].err_time << endl;
    }
    for (int i = 0; i < media_periodi_forcing.size(); i++)
    {
        fout_media_tempi_f << i + 1 << "\t" << media_periodi_forcing[i].avg_time << "\t" << media_periodi_forcing[i].err_avg_time << "\t" << media_periodi_forcing[i].err_time << endl;
    }

    //Stampa gli istogrammi dei periodi sia grezzi che forcing
    vector<vettoredoppio> hist_tempi;
    vector<vettoredoppio> hist_tempi_f;
    vector<double> bin_tempi(20, 25.0);
    gauss_hist_tempi(periodi, hist_tempi, bin_tempi);
    gauss_hist_tempi(periodi_forcing, hist_tempi_f, bin_tempi);

    for (int j = 0; j < hist_tempi.size(); j++)
    {
        ofstream fout_hist_tempi("../Dati/Hist_tempi/" + to_string(j) + "-hist.txt");
        for (int k = 0; k < hist_tempi[j].vettore2.size(); k++)
        {
            fout_hist_tempi << hist_tempi[j].vettore2[k] << "\t" << hist_tempi[j].vettore3[k] << endl;
        }
    }
    for (int j = 0; j < hist_tempi_f.size(); j++)
    {
        ofstream fout_hist_tempi_f("../Dati/Hist_tempi/" + to_string(j) + "-hist_for.txt");
        for (int k = 0; k < hist_tempi_f[j].vettore2.size(); k++)
        {
            fout_hist_tempi_f << hist_tempi_f[j].vettore2[k] << "\t" << hist_tempi_f[j].vettore3[k] << endl;
        }
    }

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
    //Genera vettore di strutture con punti di massimo
    vector<punti_massimo> punti_di_massimo_mv;
    vector<punti_massimo> punti_di_massimo_ass;  //usati sia per mv che per assoluta
    vector<punti_massimo> punti_di_massimo_root; //solo per root
    get_maxima_mv(campione, indici_dei_massimi, punti_di_massimo_mv);
    get_maxima_ass(campione, indici_dei_massimi, punti_di_massimo_ass);
    lettura_parametri(parametri);
    get_maxima_root(parametri, punti_di_massimo_root);

    for (int i = 0; i < punti_di_massimo_root.size(); i++)
    {
        ofstream fout_temp("../Dati/Maxima_root/punti_max_" + to_string(i + 1) + "_root.txt");
        for (int j = 0; j < punti_di_massimo_root[i].t_max.size(); j++)
        {
            fout_temp << punti_di_massimo_root[i].t_max[j] << "\t" << punti_di_massimo_root[i].ampl_max[j] << endl;
        }
    }

    vector<x_y> valutaz_offset_ass;
    vector<x_y> valutaz_offset_root;
    offset(punti_di_massimo_ass, valutaz_offset_ass);
    offset(punti_di_massimo_root, valutaz_offset_root);
    for (auto d : valutaz_offset_ass)
    {
        ofstream fout_off("../Dati/Offset_regime/offset_" + to_string(d.freq) + ".txt");
        for (auto c : d.offset)
        {
            fout_off << c << "\t0" << endl;
        }
    }

    int counter = 0;
    for (auto d : valutaz_offset_root)
    {
        ofstream fout_off_ro("../Dati/Offset_regime/offset_" + to_string(counter) + "_root.txt");
        for (auto c : d.offset)
        {
            fout_off_ro << c << "\t0" << endl;
        }
        counter++;
    }

    vector<x_y> picchi_picchi_assoluti;
    vector<x_y> picchi_picchi_mv;
    vector<x_y> picchi_picchi_root;

    picco_picco_ass(punti_di_massimo_ass, picchi_picchi_assoluti);
    picco_picco_mv(punti_di_massimo_mv, picchi_picchi_mv);
    picco_picco_chi(punti_di_massimo_root, picchi_picchi_root);

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
    vector<punto_regime> campana_lor_forzante;

    add_picco_medio(picchi_picchi_assoluti, campana_lorentziana_assoluta);
    add_picco_medio(picchi_picchi_mv, campana_lorentziana_mv);
    add_picco_medio(picchi_picchi_root, campana_lorentziana_root);
    add_picco_medio(picchi_picchi_root, campana_lor_forzante);

    add_omega_2(periodi, campana_lorentziana_assoluta);
    add_omega_2(periodi, campana_lorentziana_mv);
    add_omega_2(periodi, campana_lorentziana_root);
    add_omega_2(periodi_forcing, campana_lor_forzante);

    //Prende i valori del fit di gnuplot con err a posteriori

    vector<double> par_root = {1.78062, 6.06001, -0.0542075}; 
    vector<double> par_ass = {1.85208, 6.05963, -0.0560454};  
    vector<double> par_mv = {1.78114, 6.05996, -0.0542634};   
    vector<double> par_for = {};
    double err_post_root = post_lor(campana_lorentziana_root, par_root);
    double err_post_ass = post_lor(campana_lorentziana_assoluta, par_ass);
    double err_post_mv = post_lor(campana_lorentziana_mv, par_mv);

    //STAMPA LA LORENTZIANA
    ofstream fout_lor_mv("lore_mv.txt");
    ofstream fout_lor_ass("lore_ass.txt");
    ofstream fout_lor_root("lore_root.txt");
    ofstream fout_lor_for("lore_for.txt");
    ofstream fout_omega_par("omega_par_root_ass.txt");

    for (int i = 0; i < campana_lorentziana_assoluta.size(); i++) 
    {
        fout_omega_par << i + 1 << "&&$" << campana_lorentziana_root[i].theta << " \\pm " << campana_lorentziana_root[i].err_theta << "$ && $" << campana_lorentziana_assoluta[i].theta << " \\pm " << campana_lorentziana_assoluta[i].err_theta << "$ \\\\" << endl;
    }

    //salva tutti i file la lorentiana
    for (int i = 0; i < campana_lorentziana_mv.size(); i++)
    {
        fout_lor_mv << campana_lorentziana_mv[i].omega << "\t" << campana_lorentziana_mv[i].theta << "\t" << campana_lorentziana_mv[i].err_omega << "\t" << campana_lorentziana_mv[i].err_theta << "\t" << err_post_mv << endl;
    }

    for (int i = 0; i < campana_lorentziana_assoluta.size(); i++)
    {
        fout_lor_ass << campana_lorentziana_assoluta[i].omega << "\t" << campana_lorentziana_assoluta[i].theta << "\t" << campana_lorentziana_assoluta[i].err_omega << "\t" << campana_lorentziana_assoluta[i].err_theta << "\t" << err_post_ass << endl;
    }

    for (int i = 0; i < campana_lorentziana_root.size(); i++)
    {
        fout_lor_root << campana_lorentziana_root[i].omega << "\t" << campana_lorentziana_root[i].theta << "\t" << campana_lorentziana_root[i].err_omega << "\t" << campana_lorentziana_root[i].err_theta << "\t" << err_post_root << endl;
    }

    for (int i = 0; i < campana_lor_forzante.size(); i++)
    {
        fout_lor_for << campana_lor_forzante[i].omega << "\t" << campana_lor_forzante[i].theta << "\t" << campana_lor_forzante[i].err_omega << "\t" << campana_lor_forzante[i].err_theta << "\t" << err_post_mv << endl;
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

    vector<vettoredoppio> hist_ass;
    vector<double> n_bin_ass(20, 10.0); //definisce il numero di bin, cambia il secondo valore
    gauss_hist_root(punti_di_massimo_ass, hist_ass, n_bin_ass);
    ofstream fout_disp_ass("../Histogram/dispersione_sigma_ass.txt");
    for (int j = 0; j < hist_ass.size(); j++)
    {
        ofstream fout_hist_ass("../Histogram/" + to_string(j + 10) + "-hist_ass.txt");
        for (int k = 0; k < hist_ass[j].vettore2.size(); k++)
        {
            fout_hist_ass << hist_ass[j].vettore2[k] << "\t" << hist_ass[j].vettore3[k] << endl;
        }
        fout_disp_ass << j << "\t" << hist_ass[j].dispersione_amp << endl;
    }

    vector<double> compatibilty_ass;
    vector<double> compatibilty_for;
    vector<double> compatibilty_rel;
    vector<double> compatibilty_root;
    compatibilita_omega(campana_lorentziana_assoluta, campione, compatibilty_ass);
    compatibilita_omega(campana_lor_forzante, campione, compatibilty_for);
    compatibilita_omega(campana_lorentziana_root, campione, compatibilty_root);
    compatibilita_omega_frozante_sperimentale(campana_lor_forzante, campana_lorentziana_root, compatibilty_rel);


    ofstream fout_com_ass("../Dati/Compatib/comp_ass.txt");
    ofstream fout_com_for("../Dati/Compatib/comp_for.txt");
    ofstream fout_com_rel("../Dati/Compatib/comp_rel.txt");
    ofstream fout_com_root("../Dati/Compatib/comp_root.txt");

    double sum_comp_squared = 0;
    double sum_comp_squared_for = 0;
    double sum_comp_squared_rel = 0;
    double sum_comp_squared_root = 0;

    for (int i = 0; i < campana_lorentziana_assoluta.size(); i++) 
    {
        fout_com_ass << i + 1 << "\t" << campana_lorentziana_assoluta[i].omega << "\t" << campana_lorentziana_assoluta[i].err_omega << "\t" << campione[i].data_file_freq * 2. * M_PI / 1000. << "\t" << sigma_dist_uni(0.001, 1) << "\t"
                     << "\t" << compatibilty_ass[i] << endl;
        sum_comp_squared += pow(compatibilty_ass[i], 2);
    }
    for (int i = 0; i < campana_lorentziana_root.size(); i++) 
    {
        fout_com_root << i + 1 << "\t" << campana_lorentziana_root[i].omega << "\t" << campana_lorentziana_root[i].err_omega << "\t" << campione[i].data_file_freq * 2. * M_PI / 1000. << "\t" << sigma_dist_uni(0.001, 1) << "\t"
                      << "\t" << compatibilty_root[i] << endl;
        sum_comp_squared_root += pow(compatibilty_root[i], 2);
    }
    for (int i = 0; i < campana_lor_forzante.size(); i++) 
    {
        fout_com_for << i + 1 << "\t" << campana_lor_forzante[i].omega << "\t" << campana_lor_forzante[i].err_omega << "\t" << campione[i].data_file_freq * 2. * M_PI / 1000. << "\t" << sigma_dist_uni(0.001, 1) << "\t"
                     << "\t" << compatibilty_for[i] << endl;
        sum_comp_squared_for += pow(compatibilty_for[i], 2);
    }

    for (int i = 0; i < campana_lor_forzante.size(); i++) 
    {
        fout_com_rel << i + 1 << "\t" << campana_lor_forzante[i].omega << "\t" << campana_lor_forzante[i].err_omega << "\t" << campana_lorentziana_root[i].omega << "\t" << campana_lorentziana_root[i].err_omega << "\t"
                     << "\t" << compatibilty_rel[i] << endl;
        sum_comp_squared_rel += pow(compatibilty_rel[i], 2);
    }
    cout << "Somma comp squared: " << sum_comp_squared << endl;
    cout << "Somma comp squared for: " << sum_comp_squared_for << endl;
    cout << "Somma comp squared rel: " << sum_comp_squared_rel << endl;
    cout << "Somma comp squared root: " << sum_comp_squared_root << endl;

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
    get_periodi_medi(periodi_decadimento, media_periodi_decadimento); 
    reiezione(periodi_decadimento, periodi_decadimento_reiettati, media_periodi_decadimento);
    get_periodi_medi(periodi_decadimento_reiettati, media_periodi_decadimento_reiettati);

    for (int i = 0; i < periodi_decadimento_reiettati.size(); i++)
    {
        ofstream fout_periodis("../Dati/Periodi_smorzamento/pseudoperiodi_" + to_string(i + 1) + ".txt");
        for (int j = 0; j < periodi_decadimento_reiettati[i].time.size(); j++)
        {
            fout_periodis << periodi_decadimento_reiettati[i].time[j] << endl;
        }
        fout_periodis << media_periodi_decadimento_reiettati[i].avg_time << " +/- " << media_periodi_decadimento_reiettati[i].err_avg_time << endl;
    }

    vector<vettoredoppio> indici_dei_massimi_decadimento;
    vector<punti_massimo> massimi_decadimento_ass;
    vector<punti_massimo> massimi_decadimento_mv;
    get_index_maxima(campione_decadimento, tempi_decay, indici_dei_massimi_decadimento);
    omega_s_scamorza(periodi_decadimento_reiettati, omega_s_smorzamento);
    get_maxima_ass(campione_decadimento, indici_dei_massimi_decadimento, massimi_decadimento_ass);
    get_maxima_mv(campione_decadimento, indici_dei_massimi_decadimento, massimi_decadimento_mv);

    for (int i = 0; i < omega_s_smorzamento.size(); i++)
    {
        ofstream fout_omegas("../Dati/Omega_smorzamento/omegas_" + to_string(i + 1) + ".txt");
        fout_omegas << omega_s_smorzamento[i].omega_media2 << " +/- " << omega_s_smorzamento[i].err_omega_media2 << endl;
    }

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

    for (int j = 0; j < indici_dei_massimi_decadimento.size(); j++)
    {
        ofstream smorz("../Dati/periodo_smorzamento" + to_string(j) + ".txt");

        for (int i = 0; i < indici_dei_massimi_decadimento[j].vettore2.size(); i++)
        {
            smorz << campione_decadimento[j].time[indici_dei_massimi_decadimento[j].vettore2[i] - 9] << ",\t" << campione_decadimento[j].time[indici_dei_massimi_decadimento[j].vettore2[i] + 9] << ",\t" << endl;
        }
        smorz << endl;
    }

cout<<"OMEGAS"<<endl;
    int h = 1;
    for (auto c : omega_s_smorzamento)
    {
        cout << h << endl;
        h++;
        cout << c.omega_media2 << "+/-" << c.err_omega_media2 << endl;
    }

    vector<punti_massimo> punti_max_ln_root;
    vector<punti_massimo> punti_min_ln_root;
    vector<punti_massimo> punti_max_ln_ass;
    vector<punti_massimo> punti_min_ln_ass;
    vector<punti_massimo> punti_max_ln_mv;
    vector<punti_massimo> punti_min_ln_mv;

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

    for (int i = 0; i < parametri_interpolazioni_root.size(); i++)
    {
        //cout<<"Massimi "<<i+1<<endl;
        //cout << parametri_interpolazioni_root[i].a_intercetta_gamma_max<<"+/-"<<parametri_interpolazioni_root[i].err_a_post_max<<"\t"<< parametri_interpolazioni_root[i].b_angolare_gamma_max<<"+/-"<<parametri_interpolazioni_root[i].err_b_post_max<<endl;
        //cout<< "Minimi "<<i+1 <<endl;
        //cout << parametri_interpolazioni_root[i].a_intercetta_gamma_min<<"+/-"<<parametri_interpolazioni_root[i].err_a_post_min<<"\t"<< parametri_interpolazioni_root[i].b_angolare_gamma_min<<"+/-"<<parametri_interpolazioni_root[i].err_b_post_min<<endl;
        cout << i << "Max" << parametri_interpolazioni_root[i].chi_q_max << "\t" << parametri_interpolazioni_root[i].r_max << "\t" << parametri_interpolazioni_root[i].t_max << endl;
        cout << i << "Min" << parametri_interpolazioni_root[i].chi_q_min << "\t" << parametri_interpolazioni_root[i].r_min << "\t" << parametri_interpolazioni_root[i].t_min << endl;
    }

    //calcolo di omega r e omega o a partire da omega smorzamento e gamma
    vector<double> omega_o;
    vector<double> omega_r;
    vector<double> err_omega_o;
    vector<double> err_omega_r;
    vector<double> gamma_ponderata;
    vector<double> err_gamma_ponderata;
    for (int i = 0; i < omega_s_smorzamento.size(); i++)
    {
        double omega_s = omega_s_smorzamento[i].omega_media2;
        double err_omega_s = omega_s_smorzamento[i].err_omega_media2;
        double gamma_min = parametri_interpolazioni_root[i].b_angolare_gamma_min;
        double gamma_max = parametri_interpolazioni_root[i].b_angolare_gamma_max;
        double err_gamma_min = parametri_interpolazioni_root[i].err_b_post_min;
        double err_gamma_max = parametri_interpolazioni_root[i].err_b_post_max;
        vector<double> gamma_pon = {abs(gamma_min), abs(gamma_max)};
        vector<double> err_gamma_pon = {err_gamma_min, err_gamma_max};
        gamma_ponderata.push_back(media_ponderata(gamma_pon, err_gamma_pon));
        err_gamma_ponderata.push_back(errore_media_ponderata(err_gamma_pon));
        cout << "Gamma ponderata:\t" << gamma_ponderata[i] << "+/-" << err_gamma_ponderata[i] << endl;

        double dp_s_o = omega_s / sqrt(pow(omega_s, 2) + pow(gamma_ponderata[i], 2));
        double dp_g_o = gamma_ponderata[i] / sqrt(pow(omega_s, 2) + pow(gamma_ponderata[i], 2));
        double dp_s_r = omega_s / sqrt(pow(omega_s, 2) - pow(gamma_ponderata[i], 2));
        double dp_g_r = gamma_ponderata[i] / sqrt(pow(omega_s, 2) - pow(gamma_ponderata[i], 2));
        omega_o.push_back(sqrt(pow(omega_s, 2) + pow(gamma_ponderata[i], 2))); 
        omega_r.push_back(sqrt(pow(omega_s, 2) - pow(gamma_ponderata[i], 2)));
        err_omega_o.push_back(sqrt(pow(dp_s_o * err_omega_s, 2) + pow(dp_g_o * err_gamma_ponderata[i], 2)));
        err_omega_r.push_back(sqrt(pow(dp_s_r * err_omega_s, 2) + pow(dp_g_r * err_gamma_ponderata[i], 2)));
        cout << "Omega_o:\t" << omega_o[i] << "+/-" << err_omega_o[i] << endl;
        cout << "Omega_r:\t" << omega_r[i] << "+/-" << err_omega_r[i] << endl;
    }
    double mp_omega_r_riass = media_ponderata(omega_r, err_omega_r);
    double err_mp_omega_r_riass = errore_media_ponderata(err_omega_r);
    cout << "Gamma_pond_riassuntiva:\t" << media_ponderata(gamma_ponderata, err_gamma_ponderata) << "+/-" << errore_media_ponderata(err_gamma_ponderata) << endl;
    cout << "Omega_o_riassuntivo:\t" << media_ponderata(omega_o, err_omega_o) << "+/-" << errore_media_ponderata(err_omega_o) << endl;
    cout << "Omega_r_riassuntivo:\t" << mp_omega_r_riass << "+/-" << err_mp_omega_r_riass << endl;
    cout << "Comp_r_fit_decad:\t" << comp_3(mp_omega_r_riass, par_root[1], err_mp_omega_r_riass, 0.001032) << endl;

    cout << endl;
    cout << "Comp gamma 1 2: " << comp(gamma_ponderata[0], gamma_ponderata[1], err_gamma_ponderata[0], err_gamma_ponderata[1]) << endl;
    cout << "Comp gamma 1 3: " << comp(gamma_ponderata[0], gamma_ponderata[2], err_gamma_ponderata[0], err_gamma_ponderata[2]) << endl;
    cout << "Comp gamma 2 3: " << comp(gamma_ponderata[1], gamma_ponderata[2], err_gamma_ponderata[1], err_gamma_ponderata[2]) << endl;
    cout << "Comp omega_o 1 2: " << comp(omega_o[0], omega_o[1], err_omega_o[0], err_omega_o[1]) << endl;
    cout << "Comp omega_o 1 3: " << comp(omega_o[0], omega_o[2], err_omega_o[0], err_omega_o[2]) << endl;
    cout << "Comp omega_o 2 3: " << comp(omega_o[1], omega_o[2], err_omega_o[1], err_omega_o[2]) << endl;
    cout << "Comp omega_r 1 2: " << comp(omega_r[0], omega_r[1], err_omega_r[0], err_omega_r[1]) << endl;
    cout << "Comp omega_r 1 3: " << comp(omega_r[0], omega_r[2], err_omega_r[0], err_omega_r[2]) << endl;
    cout << "Comp omega_r 2 3: " << comp(omega_r[1], omega_r[2], err_omega_r[1], err_omega_r[2]) << endl;

    return 0;
}
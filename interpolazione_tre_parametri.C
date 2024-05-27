// Interpolazione a più parametri (max 3)
#include <cstdlib> // Necessario per std::exit() per gestire correttamente gli errori durante l'esecuzione del programma

//Variabili globali
vector<double> x, y, sigma_y;                                                        //Dati iniziali
int num_parametri;                                                                   //Numero parametri
double chi_piu_uno_val;                                                              //chi = 1, 2.3, 3.5, 4.7, ecc (in base a num_parametri)
double chi_quadro_min;                                                               //Chi quadro minimo assoluto
double par_1i, par_1f, par_2i, par_2f, par_3i, par_3f;                               //Range variazione parametri (per graficare)
vector<double> griglia_par1, griglia_par2, griglia_par3, griglia_chi_quadri;         //Griglia parametri e chi quadri relativi (per graficare)
double par1_def, par2_def, par3_def, sigma_par1_def, sigma_par2_def, sigma_par3_def; //parametri ed errori definitivi da interpolazione
string nome_file = "dati_da_interpolare.txt";                                        //File input dati


double funzione_interpolante(vector<double> x, double& par1, double par2, double par3, int i_esimo_valore) {
    return M_PI/2 - atan((par1 * par1 - pow(2 * M_PI * x[i_esimo_valore] * 1000, 2)) / (2 * M_PI * x[i_esimo_valore] * 1000 * par2 * 2));
}

void ricerca_parametri() {

    double numero_punti_iniziali = 100;                // MODIFICABILE
    double range_min1 = 1050000;                       //
    double range_max1 = 1150000;                       //
    double range_min2 = 20000;                         //
    double range_max2 = 40000;                         //
    double range_min3 = 0;                             //
    double range_max3 = 1;                             //
    double precisione_bisezione = 0.0000001;           //


    if ((range_max1 <= range_min1 || range_max2 <= range_min2) || (range_max3 <= range_min3))
    {
        cout << endl << "------------------------------";
        cout << endl << "| Range parametri non validi |";
        cout << endl << "------------------------------";
        cout << endl << endl; exit(1);
    }

    double passo_iniziale1 = abs(range_max1 - range_min1) / numero_punti_iniziali;      //step avanzamento cicli for sottostanti
    double passo_iniziale2 = abs(range_max2 - range_min2) / numero_punti_iniziali;
    double passo_iniziale3 = abs(range_max3 - range_min3) / numero_punti_iniziali;

    if (num_parametri < 2)                                                              //in questo modo entra nei cicli for sottostanti una sola volta per iterazione se quei parametri non sono necessari
    {
        double range_min2 = 0;
        double range_max2 = passo_iniziale2;
    }
    if (num_parametri < 3)
    {
        double range_min3 = 0;
        double range_max3 = passo_iniziale3;
    }

    double chi_min_ass_iniziale = 1000000;
    double par_i_1 = 0;
    double par_i_2 = 0;
    double par_i_3 = 0;

    cout << endl << endl;

    for (double par1 = range_min1; par1 < range_max1; par1 += passo_iniziale1)                      //par1
    {
        // percentuale di completamento
        double progress = (par1 + passo_iniziale1 - range_min1) / (range_max1 - range_min1) * 100;
        if (progress > 100) { progress = 100; }
        cout << "\rStima iniziale parametri in corso: " << setprecision(2) << progress << "%   ";

        for (double par2 = range_min2; par2 < range_max2; par2 += passo_iniziale2)                  //par2
        {
            for (double par3 = range_min3; par3 < range_max3; par3 += passo_iniziale3)              //par3
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1, par2, par3, i)) / sigma_y[i], 2);
                }
                if (sum_chi < chi_min_ass_iniziale)
                {
                    chi_min_ass_iniziale = sum_chi;
                    par_i_1 = par1;
                    par_i_2 = par2;
                    par_i_3 = par3;

                }
            }
        }
        // Pulisci la riga precedente
        cout << flush;
    }

    //range di perfezionamento della stima attraverso un metodo di bisezione
    double range_min_par1 = par_i_1 - 1.5 * passo_iniziale1;
    if (range_min_par1 < 0) { range_min_par1 = 0; }
    double range_max_par1 = par_i_1 + 2 * passo_iniziale1;
    double range_min_par2 = par_i_2 - 1.5 * passo_iniziale2;
    if (range_min_par2 < 0) { range_min_par2 = 0; }
    double range_max_par2 = par_i_2 + 2 * passo_iniziale2;
    double range_min_par3 = par_i_3 - 1.5 * passo_iniziale3;
    if (range_min_par3 < 0) { range_min_par3 = 0; }
    double range_max_par3 = par_i_3 + 2 * passo_iniziale3;


    //ricerca rapida con algoritmo di bisezione
    double par1 = par_i_1;
    double par2 = par_i_2;
    double par3 = par_i_3;

    // ----------------------------------------- PARAMETRO 1 ---------------------------------------------------------------------------------
    double precisione_par1 = precisione_bisezione * par1;
    int controllo_par1 = 0;
    vector<double> chi_par1, par_chi1;
    double min_chi_par1 = 100000;
    double sec_min_chi_par1 = 1000000;
    double posizione_min_chi_par1 = 0;
    double posizione_sec_min_chi_par1 = 1;

    for (double p = range_min_par1; p < range_max_par1 * 1.1; p += (range_max_par1 - range_min_par1))
    {
        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, p, par2, par3, i)) / sigma_y[i], 2);
        }
        chi_par1.push_back(sum_chi);
        par_chi1.push_back(p);
    }


    while (abs(par_chi1[posizione_min_chi_par1] - par_chi1[posizione_sec_min_chi_par1]) > precisione_par1)
    {

        if (par_chi1[posizione_min_chi_par1] > par_chi1[posizione_sec_min_chi_par1])
        {
            par1 = par_chi1[posizione_sec_min_chi_par1] + abs(par_chi1[posizione_min_chi_par1] - par_chi1[posizione_sec_min_chi_par1]) / 2;
        }
        else
        {
            par1 = par_chi1[posizione_min_chi_par1] + abs(par_chi1[posizione_min_chi_par1] - par_chi1[posizione_sec_min_chi_par1]) / 2;
        }

        if (num_parametri > 1)
        {
            // ----------------------------------------- PARAMETRO 2 ---------------------------------------------------------------------------------
            double precisione_par2 = precisione_bisezione * par2;
            int controllo_par2 = 0;
            vector<double> chi_par2, par_chi2;
            double min_chi_par2 = 100000;
            double sec_min_chi_par2 = 1000000;
            double posizione_min_chi_par2 = 0;
            double posizione_sec_min_chi_par2 = 1;

            for (double p = range_min_par2; p < range_max_par2 * 1.1; p += (range_max_par2 - range_min_par2))
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1, p, par3, i)) / sigma_y[i], 2);
                }
                chi_par2.push_back(sum_chi);
                par_chi2.push_back(p);
            }


            while (abs(par_chi2[posizione_min_chi_par2] - par_chi2[posizione_sec_min_chi_par2]) > precisione_par2)
            {
                if (par_chi2[posizione_min_chi_par2] > par_chi2[posizione_sec_min_chi_par2])
                {
                    par2 = par_chi2[posizione_sec_min_chi_par2] + abs(par_chi2[posizione_min_chi_par2] - par_chi2[posizione_sec_min_chi_par2]) / 2;
                }
                else
                {
                    par2 = par_chi2[posizione_min_chi_par2] + abs(par_chi2[posizione_min_chi_par2] - par_chi2[posizione_sec_min_chi_par2]) / 2;
                }

                if (num_parametri > 2)
                {
                    // ----------------------------------------- PARAMETRO 3 ---------------------------------------------------------------------------------
                    double precisione_par3 = precisione_bisezione * par3;
                    int controllo_par3 = 0;
                    vector<double> chi_par3, par_chi3;
                    double min_chi_par3 = 100000;
                    double sec_min_chi_par3 = 1000000;
                    double posizione_min_chi_par3 = 0;
                    double posizione_sec_min_chi_par3 = 1;

                    for (double p = range_min_par3; p < range_max_par3 * 1.1; p += (range_max_par3 - range_min_par3))
                    {
                        double sum_chi = 0;
                        for (int i = 0; i < x.size(); i++)
                        {
                            sum_chi += pow((y[i] - funzione_interpolante(x, par1, par2, p, i)) / sigma_y[i], 2);
                        }
                        chi_par3.push_back(sum_chi);
                        par_chi3.push_back(p);
                    }

                    while (abs(par_chi3[posizione_min_chi_par3] - par_chi3[posizione_sec_min_chi_par3]) > precisione_par3)
                    {
                        if (par_chi3[posizione_min_chi_par3] > par_chi3[posizione_sec_min_chi_par3])
                        {
                            par3 = par_chi3[posizione_sec_min_chi_par3] + abs(par_chi3[posizione_min_chi_par3] - par_chi3[posizione_sec_min_chi_par3]) / 2;
                        }
                        else
                        {
                            par3 = par_chi3[posizione_min_chi_par3] + abs(par_chi3[posizione_min_chi_par3] - par_chi3[posizione_sec_min_chi_par3]) / 2;
                        }


                        double sum_chi = 0;
                        for (int i = 0; i < x.size(); i++)
                        {
                            sum_chi += pow((y[i] - funzione_interpolante(x, par1, par2, par3, i)) / sigma_y[i], 2);
                        }
                        par_chi3.push_back(par3);
                        chi_par3.push_back(sum_chi);


                        if (sum_chi < min_chi_par3)
                        {
                            sec_min_chi_par3 = min_chi_par3;
                            min_chi_par3 = sum_chi;
                            posizione_sec_min_chi_par3 = posizione_min_chi_par3;
                            posizione_min_chi_par3 = chi_par3.size() - 1;
                        }
                        else {
                            sec_min_chi_par3 = sum_chi;
                            posizione_sec_min_chi_par3 = chi_par3.size() - 1;
                        }

                        controllo_par3++;
                        if (controllo_par3 > 100)
                        {
                            cout << endl << "------------------------------------";
                            cout << endl << "| ERRORE bisezione iniziale (par3) |";
                            cout << endl << "------------------------------------";
                            cout << endl << endl; exit(1);
                            break;
                        }
                    }
                }
                // ----------- ritorno a parametro 2 ----------------

                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1, par2, par3, i)) / sigma_y[i], 2);
                }
                par_chi2.push_back(par2);
                chi_par2.push_back(sum_chi);


                if (sum_chi < min_chi_par2)
                {
                    sec_min_chi_par2 = min_chi_par2;
                    min_chi_par2 = sum_chi;
                    posizione_sec_min_chi_par2 = posizione_min_chi_par2;
                    posizione_min_chi_par2 = chi_par2.size() - 1;
                }
                else {
                    sec_min_chi_par2 = sum_chi;
                    posizione_sec_min_chi_par2 = chi_par2.size() - 1;
                }

                controllo_par2++;
                if (controllo_par2 > 100)
                {
                    cout << endl << "------------------------------------";
                    cout << endl << "| ERRORE bisezione iniziale (par2) |";
                    cout << endl << "------------------------------------";
                    cout << endl << endl; exit(1);
                    break;
                }
            }
        }
        // ----------- ritorno a parametro 1 ----------------

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par1, par2, par3, i)) / sigma_y[i], 2);
        }
        par_chi1.push_back(par1);
        chi_par1.push_back(sum_chi);
        //cout << "sum_chi = " << sum_chi << endl;


        if (sum_chi < min_chi_par1)
        {
            sec_min_chi_par1 = min_chi_par1;
            min_chi_par1 = sum_chi;
            posizione_sec_min_chi_par1 = posizione_min_chi_par1;
            posizione_min_chi_par1 = chi_par1.size() - 1;
        }
        else {
            sec_min_chi_par1 = sum_chi;
            posizione_sec_min_chi_par1 = chi_par1.size() - 1;
        }

        controllo_par1++;
        if (controllo_par1 > 100)
        {
            cout << endl << "------------------------------------";
            cout << endl << "| ERRORE bisezione iniziale (par1) |";
            cout << endl << "------------------------------------";
            cout << endl << endl; exit(1);
            break;
        }
    }

    // fine bisezione dei 3 parametri e aggiornamento valore parametri
    par1_def = par1;
    par2_def = par2;
    par3_def = par3;

    //Trovo chi_quadro minimo
    chi_quadro_min = 0;
    for (int i = 0; i < x.size(); i++)
    {
        chi_quadro_min += pow((y[i] - funzione_interpolante(x, par1_def, par2_def, par3_def, i)) / sigma_y[i], 2);
    }

}


void numero_parametri() {
    int errore = 0;
    int p2_check = 0;
    int p3_check = 0;
    int i_esimo_val_es = 3;
    double es1 = 3.48;              //3.48 in modo che non sia un numero particolare (1 o 0) che annulli ecc la funzione
    double es2 = 3.48;
    double es3 = 3.48;
    double a = funzione_interpolante(x, es1, es2, es3, i_esimo_val_es) * 1000000000000;      //Valori molto alti affinché sia facilmente individuabile la differenza anche con "code" molto basse a valori quasi costanti

    es1 = 10000;
    es2 = 3.48;
    es3 = 3.48;
    double b = funzione_interpolante(x, es1, es2, es3, i_esimo_val_es) * 1000000000000;

    if (abs(a - b) < 0.00000001)
    {
        cout << endl << "------------------------------------------";
        cout << endl << "| ERRORE manca il primo parametro (par1) |";
        cout << endl << "------------------------------------------";
        cout << endl << endl; exit(1);
    }
    es1 = 3.48;
    es2 = 20000;
    es3 = 3.48;
    double c = funzione_interpolante(x, es1, es2, es3, i_esimo_val_es) * 1000000000000;
    if (abs(a - c) < 0.00000001)
    {
        num_parametri = 1;
        p2_check = 1;
    }

    es1 = 3.48;
    es2 = 3.48;
    es3 = 30000;
    double d = funzione_interpolante(x, es1, es2, es3, i_esimo_val_es) * 1000000000000;
    if (abs(a - d) < 0.00000001)
    {
        if (p2_check == 0)
        {
            num_parametri = 2;
        }
        p3_check = 1;
    }
    else
    {
        num_parametri = 3;
    }
    if ((p2_check == 1) && (p3_check == 0))
    {
        cout << endl << "------------------------------------";
        cout << endl << "| ERRORE presente par3 ma non par2 |";
        cout << endl << "------------------------------------";
        cout << endl << endl; exit(1);
    }
}

void chi_quadro_piu_uno() {
    if (num_parametri == 1)
    {
        chi_piu_uno_val = 1;
    }
    else
    {
        double chi_critico = 0;
        double cum_prob = 0;
        double step_size = 0.0001;
        double x = 0;

        while (cum_prob < 0.682689492137086) {
            double pdf = exp(-0.5 * x) * pow(x, (num_parametri - 2) / 2.0) / (pow(2, num_parametri / 2.0) * tgamma(num_parametri / 2.0));
            chi_critico = x;
            cum_prob += pdf * step_size;
            x += step_size;
        }
        chi_piu_uno_val = chi_critico;
    }
}


void Popolamento_vettori(vector<double>& x, vector<double>& y, vector<double>& sigma_y, string file_dati) {

    ifstream infile(file_dati);

    if (!infile.is_open()) {
        cout << endl << "------------------------------";
        cout << endl << "| Impossibile aprire il file |";
        cout << endl << "------------------------------";
        cout << endl << endl; exit(1);
    }

    // lettura dei dati dal file
    double val_x, val_y, val_sigma_y;

    while (infile >> val_x >> val_y >> val_sigma_y) {
        x.push_back(val_x);
        y.push_back(val_y);
        sigma_y.push_back(val_sigma_y);
    }
    infile.close();         // chiusura del file
}


void chi_piu_uno_par1(double& sigma_chi_piu_uno) {

    //precisione
    double precisione_f = par1_def * 0.00001;

    //Aggiorno max_dx e max_sx per indicare alla funzione chi_piu_uno gli estremi del parametro di dx e sx per cui il chi quadro è minore della variabile: double chi_massimo
    double max_val_dx = par_1f;
    double max_val_sx = par_1i;

    //Cerco il chi + 1 di destra

    vector<double> dx_f_chi, dx_chi;
    double min_chi_dx = par1_def;
    double sec_min_chi_dx = max_val_dx;
    int posizione_min_chi_dx = 0;
    int posizione_sec_min_chi_dx = 1;

    double f_t = par1_def;

    vector<double> primi_chi_dx = { par1_def, max_val_dx };
    for (int t = 0; t < primi_chi_dx.size(); t++)
    {
        double sum_chi = 0;
        f_t = primi_chi_dx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, f_t, par2_def, par3_def, i)) / sigma_y[i], 2);
        }
        dx_f_chi.push_back(f_t);
        dx_chi.push_back(sum_chi);
    }

    double precisione_f_t = precisione_f;
    int controllo = 0;

    while (abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) > precisione_f_t)
    {

        f_t = dx_f_chi[posizione_min_chi_dx] + abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2;

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, f_t, par2_def, par3_def, i)) / sigma_y[i], 2);
        }
        dx_f_chi.push_back(f_t);
        dx_chi.push_back(sum_chi);

        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_dx = sum_chi;
            posizione_sec_min_chi_dx = dx_chi.size() - 1;
        }
        else {
            min_chi_dx = sum_chi;
            posizione_min_chi_dx = dx_chi.size() - 1;
        }

        controllo++;
        if (controllo > 100)
        {
            cout << endl << "---------------------------------";
            cout << endl << "| ERRORE chi+1 di destra (par1) |";
            cout << endl << "---------------------------------";
            cout << endl << endl; exit(1);
            break;
        }

    }

    //Cerco il chi + 1 di sinistra

    vector<double> sx_f_chi, sx_chi;
    double min_chi_sx = par1_def * 0.999;
    double sec_min_chi_sx = max_val_sx;
    int posizione_min_chi_sx = 0;
    int posizione_sec_min_chi_sx = 1;

    double f_t_sx = par1_def;

    vector<double> primi_chi_sx = { par1_def * 0.999, max_val_sx };
    for (int t = 0; t < primi_chi_sx.size(); t++)
    {
        double sum_chi = 0;
        f_t_sx = primi_chi_sx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, f_t_sx, par2_def, par3_def, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(f_t_sx);
        sx_chi.push_back(sum_chi);
    }

    double precisione_f_t_sx = precisione_f;
    int controllo_sx = 0;


    while (abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) > precisione_f_t_sx)
    {

        f_t_sx = sx_f_chi[posizione_sec_min_chi_sx] + abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2;
        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, f_t_sx, par2_def, par3_def, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(f_t_sx);
        sx_chi.push_back(sum_chi);


        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_sx = sum_chi;
            posizione_sec_min_chi_sx = sx_chi.size() - 1;
        }
        else {
            min_chi_sx = sum_chi;
            posizione_min_chi_sx = sx_chi.size() - 1;
        }

        controllo_sx++;
        if (controllo_sx > 100)
        {
            cout << endl << "-----------------------------------";
            cout << endl << "| ERRORE chi+1 di sinistra (par1) |";
            cout << endl << "-----------------------------------";
            cout << endl << endl; exit(1);
            break;
        }

    }

    // Definisco il sigma finale, media dei due chi+1 trovati
    double a = dx_f_chi[posizione_min_chi_dx] + abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2; //cout << endl << "-->" << a << endl;
    double b = sx_f_chi[posizione_min_chi_sx] + abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2; //cout << endl << "-->" << b << endl;
    sigma_chi_piu_uno = (a - b) / 2;
}

void chi_piu_uno_par2(double& sigma_chi_piu_uno) {

    //precisione
    double precisione_f = par2_def * 0.00001;

    //Aggiorno max_dx e max_sx per indicare alla funzione chi_piu_uno gli estremi del parametro di dx e sx per cui il chi quadro è minore della variabile: double chi_massimo
    double max_val_dx = par_2f;
    double max_val_sx = par_2i;

    //Cerco il chi + 1 di destra

    vector<double> dx_f_chi, dx_chi;
    double min_chi_dx = par2_def;
    double sec_min_chi_dx = max_val_dx;
    int posizione_min_chi_dx = 0;
    int posizione_sec_min_chi_dx = 1;

    double f_t = par2_def;

    vector<double> primi_chi_dx = { par2_def, max_val_dx };
    for (int t = 0; t < primi_chi_dx.size(); t++)
    {
        double sum_chi = 0;
        f_t = primi_chi_dx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par3_def, f_t, par3_def, i)) / sigma_y[i], 2);
        }
        dx_f_chi.push_back(f_t);
        dx_chi.push_back(sum_chi);
    }

    double precisione_f_t = precisione_f;
    int controllo = 0;

    while (abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) > precisione_f_t)
    {

        f_t = dx_f_chi[posizione_min_chi_dx] + abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2;

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par3_def, f_t, par3_def, i)) / sigma_y[i], 2);
        }
        dx_f_chi.push_back(f_t);
        dx_chi.push_back(sum_chi);


        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_dx = sum_chi;
            posizione_sec_min_chi_dx = dx_chi.size() - 1;
        }
        else {
            min_chi_dx = sum_chi;
            posizione_min_chi_dx = dx_chi.size() - 1;
        }

        controllo++;
        if (controllo > 100)
        {
            cout << endl << "--------------------------------";
            cout << endl << "| ERRORE chi+1 di destra (par2) |";
            cout << endl << "--------------------------------";
            cout << endl << endl; exit(1);
            break;
        }

    }

    //Cerco il chi + 1 di sinistra    

    vector<double> sx_f_chi, sx_chi;
    double min_chi_sx = par2_def * 0.999;
    double sec_min_chi_sx = max_val_sx;
    int posizione_min_chi_sx = 0;
    int posizione_sec_min_chi_sx = 1;

    double f_t_sx = par2_def;

    vector<double> primi_chi_sx = { par2_def * 0.999, max_val_sx };
    for (int t = 0; t < primi_chi_sx.size(); t++)
    {
        double sum_chi = 0;
        f_t_sx = primi_chi_sx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par1_def, f_t_sx, par3_def, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(f_t_sx);
        sx_chi.push_back(sum_chi);
    }

    double precisione_f_t_sx = precisione_f;
    int controllo_sx = 0;

    while (abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) > precisione_f_t_sx)
    {

        f_t_sx = sx_f_chi[posizione_sec_min_chi_sx] + abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2;

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par1_def, f_t_sx, par3_def, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(f_t_sx);
        sx_chi.push_back(sum_chi);


        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_sx = sum_chi;
            posizione_sec_min_chi_sx = sx_chi.size() - 1;
        }
        else {
            min_chi_sx = sum_chi;
            posizione_min_chi_sx = sx_chi.size() - 1;
        }

        controllo_sx++;
        if (controllo_sx > 100)
        {
            cout << endl << "-----------------------------------";
            cout << endl << "| ERRORE chi+1 di sinistra (par2) |";
            cout << endl << "-----------------------------------";
            cout << endl << endl; exit(1);
            break;
        }

    }

    // Definisco il sigma finale, media dei due chi+1 trovati
    double a = dx_f_chi[posizione_min_chi_dx] + abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2; //cout << endl << "-->" << a << endl;
    double b = sx_f_chi[posizione_min_chi_sx] + abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2; //cout << endl << "-->" << b << endl;
    sigma_chi_piu_uno = (a - b) / 2;
}

void chi_piu_uno_par3(double& sigma_chi_piu_uno) {

    //precisione
    double precisione_f = par3_def * 0.00001;

    //Aggiorno max_dx e max_sx per indicare alla funzione chi_piu_uno gli estremi del parametro di dx e sx per cui il chi quadro è minore della variabile: double chi_massimo
    double max_val_dx = par_3f;
    double max_val_sx = par_3i;

    //Cerco il chi + 1 di destra

    vector<double> dx_f_chi, dx_chi;
    double min_chi_dx = par3_def;
    double sec_min_chi_dx = max_val_dx;
    int posizione_min_chi_dx = 0;
    int posizione_sec_min_chi_dx = 1;

    double f_t = par3_def;

    vector<double> primi_chi_dx = { par3_def, max_val_dx };
    for (int t = 0; t < primi_chi_dx.size(); t++)
    {
        double sum_chi = 0;
        f_t = primi_chi_dx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par3_def, par2_def, f_t, i)) / sigma_y[i], 2);
        }
        dx_f_chi.push_back(f_t);
        dx_chi.push_back(sum_chi);
    }


    double precisione_f_t = precisione_f;
    int controllo = 0;

    while (abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) > precisione_f_t)
    {

        f_t = dx_f_chi[posizione_min_chi_dx] + abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2;

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par3_def, par2_def, f_t, i)) / sigma_y[i], 2);
        }
        dx_f_chi.push_back(f_t);
        dx_chi.push_back(sum_chi);


        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_dx = sum_chi;
            posizione_sec_min_chi_dx = dx_chi.size() - 1;
        }
        else {
            min_chi_dx = sum_chi;
            posizione_min_chi_dx = dx_chi.size() - 1;
        }

        controllo++;
        if (controllo > 100)
        {
            cout << endl << "---------------------------------";
            cout << endl << "| ERRORE chi+1 di destra (par3) |";
            cout << endl << "---------------------------------";
            cout << endl << endl; exit(1);
            break;
        }

    }

    //Cerco il chi + 1 di sinistra    

    vector<double> sx_f_chi, sx_chi;
    double min_chi_sx = par3_def * 0.999;
    double sec_min_chi_sx = max_val_sx;
    int posizione_min_chi_sx = 0;
    int posizione_sec_min_chi_sx = 1;

    double f_t_sx = par3_def;

    vector<double> primi_chi_sx = { par3_def * 0.999, max_val_sx };
    for (int t = 0; t < primi_chi_sx.size(); t++)
    {
        double sum_chi = 0;
        f_t_sx = primi_chi_sx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par1_def, par2_def, f_t_sx, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(f_t_sx);
        sx_chi.push_back(sum_chi);
    }

    double precisione_f_t_sx = precisione_f;
    int controllo_sx = 0;

    while (abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) > precisione_f_t_sx)
    {

        f_t_sx = sx_f_chi[posizione_sec_min_chi_sx] + abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2;

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par1_def, par2_def, f_t_sx, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(f_t_sx);
        sx_chi.push_back(sum_chi);


        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_sx = sum_chi;
            posizione_sec_min_chi_sx = sx_chi.size() - 1;
        }
        else {
            min_chi_sx = sum_chi;
            posizione_min_chi_sx = sx_chi.size() - 1;
        }

        controllo_sx++;
        if (controllo_sx > 100)
        {
            cout << endl << "-----------------------------------";
            cout << endl << "| ERRORE chi+1 di sinistra (par3) |";
            cout << endl << "-----------------------------------";
            cout << endl << endl; exit(1);
            break;
        }

    }

    // Definisco il sigma finale, media dei due chi+1 trovati
    double a = dx_f_chi[posizione_min_chi_dx] + abs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2; //cout << endl << "-->" << a << endl;
    double b = sx_f_chi[posizione_min_chi_sx] + abs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2; //cout << endl << "-->" << b << endl;
    sigma_chi_piu_uno = (a - b) / 2;
}


void risultato(double valore, double errore, string nome) {
    cout << endl << nome << " = ( " << fixed << setprecision(10) << valore << " +- " << setprecision(10) << errore << " )" << endl;

    // Includere:   #include <iomanip>
    int digits_val = 0 - floor(log10(errore));
    if (errore >= 1) { digits_val = 0; }
    cout << endl << nome << " = ( " << fixed << setprecision(digits_val) << valore << " +- " << setprecision(digits_val) << errore << " )" << endl;
}

void calcolo_matrice_chi_quadri(double fattore_dimensione_griglia, double punti_per_parametro, int coppia_parametri) {

    
    
    if (fattore_dimensione_griglia >= 5)
    {
        double chi_massimo = fattore_dimensione_griglia;

        double par1 = par1_def;
        double par2 = par2_def;
        double par3 = par3_def;

        //PAR1 - Ricerca iniziale regione con chi_quadro < 100
        double par1_chi_massimo_iniziale1 = 0;
        double par1_chi_massimo_iniziale2 = 0;
        if (1 == 1) {

            double sum_chi_0 = 0;
            double par_zero = 0;
            for (int i = 0; i < x.size(); i++)
            {
                sum_chi_0 += pow((y[i] - funzione_interpolante(x, par_zero, par2, par3, i)) / sigma_y[i], 2);
            }

            double passo_iniziale;              //passo di scansione per trovare regione minore del chi_massimo necessario (es. |chi|<100)
            if (par1 > 3)
            {
                passo_iniziale = 0.01;
            }
            else
            {
                passo_iniziale = 0.001 * par1;
            }

            if (sum_chi_0 < chi_massimo)
            {
                par1_chi_massimo_iniziale1 = 0;
            }
            else
            {
                vector<double> par1_chi_massimo_iniz, chi_chi_massimo;
                double s = 0;
                chi_chi_massimo.push_back(chi_massimo * 10);
                par1_chi_massimo_iniz.push_back(s);
                s++;
                int indice = 0;
                while (chi_chi_massimo[indice] > chi_massimo)
                {
                    double sum_chi_s = 0;
                    for (int i = 0; i < x.size(); i++)
                    {
                        sum_chi_s += pow((y[i] - funzione_interpolante(x, s, par2, par3, i)) / sigma_y[i], 2);
                    }
                    par1_chi_massimo_iniz.push_back(s);
                    chi_chi_massimo.push_back(sum_chi_s);
                    s += passo_iniziale;
                    indice++;
                }
                par1_chi_massimo_iniziale1 = par1_chi_massimo_iniz[par1_chi_massimo_iniz.size() - 2];           //trovato estremo range di SX
            }

            //Ricerca fine dove di nuovo il chi quadro comincia a valere più del range necessario deciso (es. |chi|<100)
            for (double k = par1_chi_massimo_iniziale1; k < 100; k += passo_iniziale)
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, k, par2, par3, i)) / sigma_y[i], 2);
                }

                if (sum_chi > chi_massimo && k > (par1_chi_massimo_iniziale1 + 2 * passo_iniziale)) {
                    par1_chi_massimo_iniziale2 = k;                                                             //trovato estremo range di DX
                    break;
                }
            }





        }

        //PAR2 - Ricerca iniziale regione con chi_quadro < 100
        double par2_chi_massimo_iniziale1 = 0;
        double par2_chi_massimo_iniziale2 = 0;
        if (num_parametri > 1) {

            double sum_chi_0 = 0;
            double par_zero = 0;
            for (int i = 0; i < x.size(); i++)
            {
                sum_chi_0 += pow((y[i] - funzione_interpolante(x, par1, par_zero, par3, i)) / sigma_y[i], 2);
            }

            double passo_iniziale;              //passo di scansione per trovare regione minore del chi_massimo necessario (es. |chi|<100)
            if (par1 > 3)
            {
                passo_iniziale = 0.01;
            }
            else
            {
                passo_iniziale = 0.001 * par1;
            }

            if (sum_chi_0 < chi_massimo)
            {
                par2_chi_massimo_iniziale1 = 0;
            }
            else
            {
                vector<double> par1_chi_massimo_iniz, chi_chi_massimo;
                double s = 0;
                chi_chi_massimo.push_back(chi_massimo * 10);
                par1_chi_massimo_iniz.push_back(s);
                s++;
                int indice = 0;
                while (chi_chi_massimo[indice] > chi_massimo)
                {
                    double sum_chi_s = 0;
                    for (int i = 0; i < x.size(); i++)
                    {
                        sum_chi_s += pow((y[i] - funzione_interpolante(x, par1, s, par3, i)) / sigma_y[i], 2);
                    }
                    par1_chi_massimo_iniz.push_back(s);
                    chi_chi_massimo.push_back(sum_chi_s);
                    s += passo_iniziale;
                    indice++;
                }
                par2_chi_massimo_iniziale1 = par1_chi_massimo_iniz[par1_chi_massimo_iniz.size() - 2];           //trovato estremo range di SX
            }

            //Ricerca fine dove di nuovo il chi quadro comincia a valere più del range necessario deciso (es. |chi|<100)
            for (double k = par1_chi_massimo_iniziale1; k < 100; k += passo_iniziale)
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1, k, par3, i)) / sigma_y[i], 2);
                }

                if (sum_chi > chi_massimo && k > (par1_chi_massimo_iniziale1 + 2 * passo_iniziale)) {
                    par2_chi_massimo_iniziale2 = k;                                                             //trovato estremo range di DX
                    break;
                }
            }





        }


        //PAR3 - Ricerca iniziale regione con chi_quadro < 100
        double par3_chi_massimo_iniziale1 = 0;
        double par3_chi_massimo_iniziale2 = 0;
        if (num_parametri > 2) {

            double sum_chi_0 = 0;
            double par_zero = 0;
            for (int i = 0; i < x.size(); i++)
            {
                sum_chi_0 += pow((y[i] - funzione_interpolante(x, par1, par2, par_zero, i)) / sigma_y[i], 2);
            }

            double passo_iniziale;              //passo di scansione per trovare regione minore del chi_massimo necessario (es. |chi|<100)
            if (par1 > 3)
            {
                passo_iniziale = 0.01;
            }
            else
            {
                passo_iniziale = 0.001 * par1;
            }

            if (sum_chi_0 < chi_massimo)
            {
                par3_chi_massimo_iniziale1 = 0;
            }
            else
            {
                vector<double> par1_chi_massimo_iniz, chi_chi_massimo;
                double s = 0;
                chi_chi_massimo.push_back(chi_massimo * 10);
                par1_chi_massimo_iniz.push_back(s);
                s++;
                int indice = 0;
                while (chi_chi_massimo[indice] > chi_massimo)
                {
                    double sum_chi_s = 0;
                    for (int i = 0; i < x.size(); i++)
                    {
                        sum_chi_s += pow((y[i] - funzione_interpolante(x, par1, par2, s, i)) / sigma_y[i], 2);
                    }
                    par1_chi_massimo_iniz.push_back(s);
                    chi_chi_massimo.push_back(sum_chi_s);
                    s += passo_iniziale;
                    indice++;
                }
                par3_chi_massimo_iniziale1 = par1_chi_massimo_iniz[par1_chi_massimo_iniz.size() - 2];           //trovato estremo range di SX
            }

            //Ricerca fine dove di nuovo il chi quadro comincia a valere più del range necessario deciso (es. |chi|<100)
            for (double k = par1_chi_massimo_iniziale1; k < 100; k += passo_iniziale)
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1, par2, k, i)) / sigma_y[i], 2);
                }

                if (sum_chi > chi_massimo && k > (par1_chi_massimo_iniziale1 + 2 * passo_iniziale)) {
                    par3_chi_massimo_iniziale2 = k;                                                             //trovato estremo range di DX
                    break;
                }
            }
        }

        //Trovati i range li aggiorno e per andare meglio li porto al livello globale per le altre funzioni (senza dover passare un botto di roba in ogni funzione)
        par_1i = par1_chi_massimo_iniziale1;
        par_1f = par1_chi_massimo_iniziale2;
        par_2i = par2_chi_massimo_iniziale1;
        par_2f = par2_chi_massimo_iniziale2;
        par_3i = par3_chi_massimo_iniziale1;
        par_3f = par3_chi_massimo_iniziale2;
    }

    else if (fattore_dimensione_griglia == 0)
    {
        cout << endl << "-----------------------------------------------------------";
        cout << endl << "| ERRORE fattore_dimensione_griglia valore non valido (0) |";
        cout << endl << "-----------------------------------------------------------";
        cout << endl << endl; exit(1);
    }

    else
    {
        //Faccio come prima ma aggiorno in base al numero di sigma (relazionati ai corrispondenti parametri)
        par_1i = par1_def - fattore_dimensione_griglia * sigma_par1_def;
        par_1f = par1_def + fattore_dimensione_griglia * sigma_par1_def;
        par_2i = par2_def - fattore_dimensione_griglia * sigma_par2_def;
        par_2f = par2_def + fattore_dimensione_griglia * sigma_par2_def;
        par_3i = par3_def - fattore_dimensione_griglia * sigma_par3_def;
        par_3f = par3_def + fattore_dimensione_griglia * sigma_par3_def;
    }

    //Creazione griglia scansione completa definitiva da graficare --------------------------------------------------------------------------

    double chi_min_ass_iniziale = 1000000000;

    double passo_par1 = (par_1f - par_1i) / punti_per_parametro;
    double passo_par2 = (par_2f - par_2i) / punti_per_parametro;
    double passo_par3 = (par_3f - par_3i) / punti_per_parametro;

    if (num_parametri < 2)
    {
        par_2f = par_2i + passo_par2;
    }
    if (num_parametri < 3)
    {
        par_3f = par_3i + passo_par3;
    }


    if (coppia_parametri == 12 || coppia_parametri == 21)
    {
        double par3_g = par3_def;
        for (double par1_g = par_1i; par1_g < par_1f; par1_g += passo_par1)                  //par1
        {
            // percentuale di completamento
            float progress = static_cast<float>(par1_g - par_1i + passo_par1) / (par_1f - par_1i) * 100;
            if (progress > 100) { progress = 100; }
            cout << "\rCalcolo scansione completa per graficare in corso: " << setprecision(2) << progress << "%";

            for (double par2_g = par_2i; par2_g < par_2f; par2_g += passo_par2)              //par2
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1_g, par2_g, par3_g, i)) / sigma_y[i], 2);
                }
                if (sum_chi < 1000000000)                       //MODIFICABILE: per avere solo chi quadri minori di un certo valore indipendentemente dalle combinazioni lineari dei parametri
                {
                    griglia_par1.push_back(par1_g);
                    if (num_parametri > 1)
                    {
                        griglia_par2.push_back(par2_g);
                    }
                    else
                    {
                        griglia_par2.push_back(0);
                    }

                    if (num_parametri > 2)
                    {
                        griglia_par3.push_back(par3_g);
                    }
                    else
                    {
                        griglia_par3.push_back(0);
                    }
                    griglia_chi_quadri.push_back(sum_chi);
                    //cout << par1_g << "\t" << par2_g << "\t" << par3_g << endl;
                }
            }
            // Pulisci la riga precedente
            cout << flush;
        }
    }

    else if (coppia_parametri == 23 || coppia_parametri == 32)
    {
        double par1_g = par1_def; //cout << "Ti si entra qui dentro?" << endl;
        for (double par2_g = par_2i; par2_g < par_2f; par2_g += passo_par2)              //par2
        {
            // percentuale di completamento
            float progress = static_cast<float>(par2_g - par_2i + passo_par2) / (par_2f - par_2i) * 100;
            if (progress > 100) { progress = 100; }
            cout << "\rCalcolo scansione completa per graficare in corso: " << setprecision(2) << progress << "%";

            for (double par3_g = par_3i; par3_g < par_3f; par3_g += passo_par3)          //par3
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1_g, par2_g, par3_g, i)) / sigma_y[i], 2);
                }
                if (sum_chi < 1000000000)                       //MODIFICABILE: per avere solo chi quadri minori di un certo valore indipendentemente dalle combinazioni lineari dei parametri
                {
                    griglia_par1.push_back(par1_g);
                    if (num_parametri > 1)
                    {
                        griglia_par2.push_back(par2_g);
                    }
                    else
                    {
                        griglia_par2.push_back(0);
                    }

                    if (num_parametri > 2)
                    {
                        griglia_par3.push_back(par3_g);
                    }
                    else
                    {
                        griglia_par3.push_back(0);
                    }
                    griglia_chi_quadri.push_back(sum_chi);
                }
            }
            // Pulisci la riga precedente
            cout << flush;
        }
    }

    else if (coppia_parametri == 31 || coppia_parametri == 13)
    {
        double par2_g = par2_def;
        for (double par1_g = par_1i; par1_g < par_1f; par1_g += passo_par1)                  //par1
        {
            // percentuale di completamento
            float progress = static_cast<float>(par1_g - par_1i + passo_par1) / (par_1f - par_1i) * 100;
            if (progress > 100) { progress = 100; }
            cout << "\rCalcolo scansione completa per graficare in corso: " << setprecision(2) << progress << "%";

            for (double par3_g = par_3i; par3_g < par_3f; par3_g += passo_par3)          //par3
            {
                double sum_chi = 0;
                for (int i = 0; i < x.size(); i++)
                {
                    sum_chi += pow((y[i] - funzione_interpolante(x, par1_g, par2_g, par3_g, i)) / sigma_y[i], 2);
                }
                if (sum_chi < 1000000000)                       //MODIFICABILE: per avere solo chi quadri minori di un certo valore indipendentemente dalle combinazioni lineari dei parametri
                {
                    griglia_par1.push_back(par1_g);
                    if (num_parametri > 1)
                    {
                        griglia_par2.push_back(par2_g);
                    }
                    else
                    {
                        griglia_par2.push_back(0);
                    }

                    if (num_parametri > 2)
                    {
                        griglia_par3.push_back(par3_g);
                    }
                    else
                    {
                        griglia_par3.push_back(0);
                    }
                    griglia_chi_quadri.push_back(sum_chi);
                }
            }
            // Pulisci la riga precedente
            cout << flush;
        }
    }

    else
    {
        for (double par1_g = par_1i; par1_g < par_1f; par1_g += passo_par1)                  //par1
        {
            // percentuale di completamento
            float progress = static_cast<float>(par1_g - par_1i + passo_par1) / (par_1f - par_1i) * 100;
            if (progress > 100) { progress = 100; }
            cout << "\rCalcolo scansione completa per graficare in corso: " << setprecision(2) << progress << "%";

            for (double par2_g = par_2i; par2_g < par_2f; par2_g += passo_par2)              //par2
            {
                for (double par3_g = par_3i; par3_g < par_3f; par3_g += passo_par3)          //par3
                {
                    double sum_chi = 0;
                    for (int i = 0; i < x.size(); i++)
                    {
                        sum_chi += pow((y[i] - funzione_interpolante(x, par1_g, par2_g, par3_g, i)) / sigma_y[i], 2);
                    }
                    if (sum_chi < 1000000000)                       //MODIFICABILE: per avere solo chi quadri minori di un certo valore indipendentemente dalle combinazioni lineari dei parametri
                    {
                        griglia_par1.push_back(par1_g);
                        if (num_parametri > 1)
                        {
                            griglia_par2.push_back(par2_g);
                        }
                        else
                        {
                            griglia_par2.push_back(0);
                        }

                        if (num_parametri > 2)
                        {
                            griglia_par3.push_back(par3_g);
                        }
                        else
                        {
                            griglia_par3.push_back(0);
                        }
                        griglia_chi_quadri.push_back(sum_chi);
                    }

                    if (sum_chi < chi_min_ass_iniziale)
                    {
                        chi_min_ass_iniziale = sum_chi;
                        par1_def = par1_g;
                        par2_def = par2_g;
                        par3_def = par3_g;
                    }
                }
            }
            // Pulisci la riga precedente
            cout << flush;
        }
    }
    cout << endl;
}



void Programma_interpolazione_tre_parametri()
{
    // accedo ai dati contenuti nel file di input
    Popolamento_vettori(x, y, sigma_y, nome_file);

    // Cerco di capire quanti parametri sono stati immessi nella funzione interpolante iniziale
    numero_parametri();

    // Calcolo il chi quadro = 1, 2.3, 3.5, 4.7, ecc (in base a num_parametri)
    chi_quadro_piu_uno(); //chi_piu_uno_val = 1;  //Voglio il chi+1 pari a 1 anche se ci sono più parametri

    // Stima dei parametri d'interpolazione
    ricerca_parametri();

    // Stima degli errori legati ai parametri d'interpolazione
                             chi_piu_uno_par1(sigma_par1_def);
    if (num_parametri > 1) { chi_piu_uno_par2(sigma_par2_def); }
    if (num_parametri > 2) { chi_piu_uno_par3(sigma_par3_def); }

    // Stampo a schermo dei risultati ottenuti
    cout << endl << endl;
    cout << "----------------------------------------------" << endl;
    cout << "Risultati parametri (con errore composto per piu' parametri):" << endl;
                             risultato(par1_def, sigma_par1_def, "par1"); cout << endl;
    if (num_parametri > 1) { risultato(par2_def, sigma_par2_def, "par2"); cout << endl; }
    if (num_parametri > 2) { risultato(par3_def, sigma_par3_def, "par3"); cout << endl; }
    cout << "----------------------------------------------" << endl;
    cout << "Chi quadro minimo: " << setprecision(10) << chi_quadro_min << endl << endl;
}



void interpolazione_tre_parametri() {

    //Lancia il programma d'interpolazione per stimare parametri ed errori accessibili da variabili globali
    Programma_interpolazione_tre_parametri();

    gStyle->SetOptStat(0);

    if (num_parametri == 1)
    {
        TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 1400, 500);
        c1->Divide(2, 1, 0.01, 0.01);

        //Grafico 1 ------------------------------------------------------------------------------

        vector<double> griglia_par1_1, griglia_chi_quadri_1;

        for (double i = par1_def - 3 * sigma_par1_def; i < par1_def + 3 * sigma_par1_def; i += sigma_par1_def * 6 / 40) {
            double sum_chi = 0;
            double par2_s = 0;
            double par3_s = 0;
            for (int s = 0; s < x.size(); s++)
            {
                sum_chi += pow((y[s] - funzione_interpolante(x, i, par2_s, par3_s, s)) / sigma_y[s], 2);
            }
            griglia_par1_1.push_back(i);
            griglia_chi_quadri_1.push_back(sum_chi);
        }


        TGraph* g1_chi_quadri = new TGraph();
        for (int i = 0; i < griglia_par1_1.size(); i++) {
            g1_chi_quadri->SetPoint(i, griglia_par1_1[i], griglia_chi_quadri_1[i]);
        }

        g1_chi_quadri->SetTitle("Grafico interpolazione");
        g1_chi_quadri->GetXaxis()->SetTitle("Frequenza [kHz]");
        g1_chi_quadri->GetYaxis()->SetTitle("V_in/V_out");
        g1_chi_quadri->SetMarkerStyle(22);
        g1_chi_quadri->SetMarkerSize(0.7);
        g1_chi_quadri->SetLineWidth(1);


        TMultiGraph* grafico_1 = new TMultiGraph("Grafico interpolazione", "Grafico interpolazione");
        grafico_1->Add(g1_chi_quadri, "pc");



        //Grafico 2 ------------------------------------------------------------------------------

        TGraphErrors* dati_plot = new TGraphErrors();
        for (int i = 0; i < x.size(); i++) {
            dati_plot->SetPoint(i, x[i], y[i]);
            dati_plot->SetPointError(i, 0, sigma_y[i]);
        }

        dati_plot->SetTitle("Grafico interpolazione");
        dati_plot->GetXaxis()->SetTitle("Frequenza [kHz]");
        dati_plot->GetYaxis()->SetTitle("V_in/V_out");
        dati_plot->SetMarkerStyle(22);
        dati_plot->SetMarkerSize(0.7);
        dati_plot->SetLineWidth(1);


        TGraph* funzione_interpolante_dati = new TGraph();

        vector<double> funz_interpolante;
        for (double i = x[0]; i < x[x.size() - 1]; i += abs(x[0] - x[x.size() - 1]) / 50) {
            funz_interpolante.push_back(i);
        }

        for (int i = 0; i < funz_interpolante.size(); i++)
        {
            funzione_interpolante_dati->SetPoint(i, funz_interpolante[i], funzione_interpolante(funz_interpolante, par1_def, par2_def, par3_def, i));
        }


        auto legend = new TLegend(0.7, 0.85, 0.95, 0.95);
        legend->AddEntry(dati_plot, "Dati sperimentali", "p");
        legend->AddEntry(funzione_interpolante_dati, "fit globale", "l");


        TMultiGraph* grafico_2 = new TMultiGraph("Grafico interpolazione", "Grafico interpolazione");
        grafico_2->Add(dati_plot, "p");
        grafico_2->Add(funzione_interpolante_dati, "c");


        //fine grafici ---------------------------------------------------------------------------

        c1->cd(1);
        grafico_1->Draw("a");
        c1->cd(2);
        //gPad->SetLogx();              //Imposta scala logaritmica sulle x
        grafico_2->Draw("a");
        legend->Draw();
        c1->cd(1);
        gPad->Pop();

        c1->Update();

    }


    if (num_parametri == 2)
    {
        TCanvas* c2 = new TCanvas("c2", "c2", 0, 0, 1400, 500);
        c2->Divide(2, 0, 0.01, 0.01);

        //Grafico 1 ------------------------------------------------------------------------------

        double fattore_dimensione_griglia1 = 3;      //se n<=5 vengono salvati parametri con n sigma di differenza dalla migliore stima, se n>5 vengono salvati parametri che garantiscono lungo le 3 rette principali (perpendicolari tra loro passanti per le migliori stime dei parametri) valori che restituiscono chi quadri minori di n
        double punti_per_parametro1 = 40;
        int coppia_parametri1 = 12;                  //scegliere tra 12 (o 21), 23 (o 32), 31 (o 13) oppure un valore differente per prenderli tutti e tre (per grafico 3D)
        calcolo_matrice_chi_quadri(fattore_dimensione_griglia1, punti_per_parametro1, coppia_parametri1);

        TH2F* grafico_1 = new TH2F("Parametri 1-2", "Parametri 1-2", punti_per_parametro1, par_1i, par_1f, punti_per_parametro1, par_2i, par_2f); // Definisci i bin e l'intervallo

        grafico_1->SetStats(0);	//Rimuovo il box

        // Riempi l'istogramma con i parametri e i chi_quadri associati
        for (int i = 0; i < griglia_par1.size(); i++) {
            int bin = grafico_1->FindBin(griglia_par1[i], griglia_par2[i]); // Trova il bin corrispondente alle coordinate (x, y)
            grafico_1->SetBinContent(bin, griglia_chi_quadri[i]); // Imposta il valore associato al bin
        }

        grafico_1->SetContour(1000);

        grafico_1->GetXaxis()->SetTitle("Asse x");
        grafico_1->GetYaxis()->SetTitle("Asse y");
        grafico_1->GetZaxis()->SetTitle("Chi Quadri");


        //Grafico 2 ------------------------------------------------------------------------------

        TGraphErrors* dati_plot = new TGraphErrors();
        for (int i = 0; i < x.size(); i++) {
            dati_plot->SetPoint(i, x[i], y[i]);
            dati_plot->SetPointError(i, 0, sigma_y[i]);
        }

        dati_plot->SetTitle("Grafico interpolazione");
        dati_plot->GetXaxis()->SetTitle("Frequenza [kHz]");
        dati_plot->GetYaxis()->SetTitle("V_in/V_out");
        dati_plot->SetMarkerStyle(22);
        dati_plot->SetMarkerSize(0.7);
        dati_plot->SetLineWidth(1);



        TGraph* funzione_interpolante_dati = new TGraph();

        vector<double> funz_interpolante;
        for (double i = x[0]; i < x[x.size() - 1]; i += abs(x[0] - x[x.size() - 1]) / 50) {
            funz_interpolante.push_back(i);
        }

        for (int i = 0; i < funz_interpolante.size(); i++)
        {
            funzione_interpolante_dati->SetPoint(i, funz_interpolante[i], funzione_interpolante(funz_interpolante, par1_def, par2_def, par3_def, i));
        }


        auto legend = new TLegend(0.7, 0.85, 0.95, 0.95);
        legend->AddEntry(dati_plot, "Dati sperimentali", "p");
        legend->AddEntry(funzione_interpolante_dati, "fit globale", "l");


        TMultiGraph* grafico_2 = new TMultiGraph("Grafico interpolazione", "Grafico interpolazione");
        grafico_2->Add(dati_plot, "p");
        grafico_2->Add(funzione_interpolante_dati, "c");


        //fine grafici ---------------------------------------------------------------------------

        c2->cd(1);
        grafico_1->Draw("colz");
        c2->cd(2);
        grafico_2->Draw("a");
        legend->Draw();
        c2->cd(1);
        gPad->Pop();

        c2->Update();
    }


    if (num_parametri == 3)
    {
        TCanvas* c3 = new TCanvas("c3", "c3", 0, 0, 1500, 1000);
        c3->Divide(2, 2, 0.02, 0.01);

        //Grafico 1 ------------------------------------------------------------------------------

        double fattore_dimensione_griglia1 = 3;      //se n<=5 vengono salvati parametri con n sigma di differenza dalla migliore stima, se n>5 vengono salvati parametri che garantiscono lungo le 3 rette principali (perpendicolari tra loro passanti per le migliori stime) valori che restituiscono chi quadri minori di n
        double punti_per_parametro1 = 40;
        int coppia_parametri1 = 12;                  //scegliere tra 12 (o 21), 23 (o 32), 31 (o 13) oppure un valore differente per prenderli tutti e tre (per grafico 3D)
        calcolo_matrice_chi_quadri(fattore_dimensione_griglia1, punti_per_parametro1, coppia_parametri1);

        TH2F* grafico_1 = new TH2F("Parametri 1-2", "Parametri 1-2", punti_per_parametro1, par_1i, par_1f, punti_per_parametro1, par_2i, par_2f); // Definisci i bin e l'intervallo

        grafico_1->SetStats(0);	//Rimuovo il box

        // Riempi l'istogramma con i parametri e i chi_quadri associati
        for (int i = 0; i < griglia_par1.size(); i++) {
            int bin = grafico_1->FindBin(griglia_par1[i], griglia_par2[i]); // Trova il bin corrispondente alle coordinate (x, y)
            grafico_1->SetBinContent(bin, griglia_chi_quadri[i]); // Imposta il valore associato al bin
        }

        grafico_1->SetContour(1000);

        grafico_1->GetXaxis()->SetTitle("Asse x");
        grafico_1->GetYaxis()->SetTitle("Asse y");
        grafico_1->GetZaxis()->SetTitle("Chi Quadri");


        //Grafico 2 ------------------------------------------------------------------------------

        double fattore_dimensione_griglia2 = 3;      //se n<=5 vengono salvati parametri con n sigma di differenza dalla migliore stima, se n>5 vengono salvati parametri che garantiscono lungo le 3 rette principali (perpendicolari tra loro passanti per le migliori stime) valori che restituiscono chi quadri minori di n
        double punti_per_parametro2 = 40;
        int coppia_parametri2 = 23;                  //scegliere tra 12 (o 21), 23 (o 32), 31 (o 13) oppure un valore differente per prenderli tutti e tre (per grafico 3D)
        calcolo_matrice_chi_quadri(fattore_dimensione_griglia1, punti_per_parametro1, coppia_parametri1);

        TH2F* grafico_2 = new TH2F("Parametri 2-3", "Parametri 2-3", punti_per_parametro2, par_2i, par_2f, punti_per_parametro2, par_3i, par_3f); // Definisci i bin e l'intervallo

        grafico_2->SetStats(0);	//Rimuovo il box

        // Riempi l'istogramma con i parametri e i chi_quadri associati
        for (int i = 0; i < griglia_par1.size(); i++) {
            int bin = grafico_2->FindBin(griglia_par2[i], griglia_par3[i]); // Trova il bin corrispondente alle coordinate (x, y)
            grafico_2->SetBinContent(bin, griglia_chi_quadri[i]); // Imposta il valore associato al bin
        }

        grafico_2->SetContour(1000);

        grafico_2->GetXaxis()->SetTitle("Asse x");
        grafico_2->GetYaxis()->SetTitle("Asse y");
        grafico_2->GetZaxis()->SetTitle("Chi Quadri");

        //Grafico 3 ------------------------------------------------------------------------------

        double fattore_dimensione_griglia3 = 3;      //se n<=5 vengono salvati parametri con n sigma di differenza dalla migliore stima, se n>5 vengono salvati parametri che garantiscono lungo le 3 rette principali (perpendicolari tra loro passanti per le migliori stime) valori che restituiscono chi quadri minori di n
        double punti_per_parametro3 = 40;
        int coppia_parametri3 = 31;                  //scegliere tra 12 (o 21), 23 (o 32), 31 (o 13) oppure un valore differente per prenderli tutti e tre (per grafico 3D)
        calcolo_matrice_chi_quadri(fattore_dimensione_griglia1, punti_per_parametro1, coppia_parametri1);

        TH2F* grafico_3 = new TH2F("Parametri 3-1", "Parametri 3-1", punti_per_parametro3, par_3i, par_3f, punti_per_parametro3, par_1i, par_1f); // Definisci i bin e l'intervallo

        grafico_3->SetStats(0);	//Rimuovo il box

        // Riempi l'istogramma con i parametri e i chi_quadri associati
        for (int i = 0; i < griglia_par1.size(); i++) {
            int bin = grafico_3->FindBin(griglia_par3[i], griglia_par1[i]); // Trova il bin corrispondente alle coordinate (x, y)
            grafico_3->SetBinContent(bin, griglia_chi_quadri[i]); // Imposta il valore associato al bin
        }

        grafico_3->SetContour(1000);

        grafico_3->GetXaxis()->SetTitle("Asse x");
        grafico_3->GetYaxis()->SetTitle("Asse y");
        grafico_3->GetZaxis()->SetTitle("Chi Quadri");

        //Grafico 4 ------------------------------------------------------------------------------

        TGraphErrors* dati_plot = new TGraphErrors();
        for (int i = 0; i < x.size(); i++) {
            dati_plot->SetPoint(i, x[i], y[i]);
            dati_plot->SetPointError(i, 0, sigma_y[i]);
        }

        dati_plot->SetTitle("Grafico interpolazione");
        dati_plot->GetXaxis()->SetTitle("Frequenza [kHz]");
        dati_plot->GetYaxis()->SetTitle("V_in/V_out");
        dati_plot->SetMarkerStyle(22);
        dati_plot->SetMarkerSize(0.7);
        dati_plot->SetLineWidth(1);


        TGraph* funzione_interpolante_dati = new TGraph();

        vector<double> funz_interpolante;
        for (double i = x[0]; i < x[x.size() - 1]; i += abs(x[0] - x[x.size() - 1]) / 50) {
            funz_interpolante.push_back(i);
        }

        for (int i = 0; i < funz_interpolante.size(); i++)
        {
            funzione_interpolante_dati->SetPoint(i, funz_interpolante[i], funzione_interpolante(funz_interpolante, par1_def, par2_def, par3_def, i));
        }


        auto legend = new TLegend(0.7, 0.85, 0.95, 0.95);
        legend->AddEntry(dati_plot, "Dati sperimentali", "p");
        legend->AddEntry(funzione_interpolante_dati, "fit globale", "l");


        TMultiGraph* grafico_4 = new TMultiGraph("Grafico interpolazione", "Grafico interpolazione");
        grafico_4->Add(dati_plot, "p");
        grafico_4->Add(funzione_interpolante_dati, "c");



        //fine grafici ---------------------------------------------------------------------------


        c3->cd(1);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_1->Draw("colz");
        c3->cd(2);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_2->Draw("colz");
        c3->cd(3);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_3->Draw("colz");
        c3->cd(4);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_4->Draw("a");
        legend->Draw();
        c3->cd(1);
        gPad->Pop();


        c3->Update();
    }
}


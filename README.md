# Interpolazione a più parametri

Questo programma implementa un algoritmo di interpolazione a più parametri (fino a 3). 

## Funzionalità

* **Lettura dati:** Il programma legge i dati da un file di testo, in cui ogni riga rappresenta un punto dati con le seguenti colonne:
    * `val_x`: valore della variabile indipendente
    * `val_y`: valore della variabile dipendente
    * `val_sigma_y`: errore sulla variabile dipendente
* **Determinazione del numero di parametri:** Il programma identifica automaticamente il numero di parametri necessari per l'interpolazione, analizzando la funzione interpolante.
* **Ricerca dei parametri:** Il programma utilizza un algoritmo di bisezione per trovare i parametri che minimizzano il chi-quadro.
* **Calcolo degli errori:** Il programma calcola gli errori sui parametri utilizzando la defnizione di standard error.
* **Grafici:** Il programma genera grafici che mostrano i dati, la funzione interpolante e il chi-quadro in funzione dei parametri.

## Come usare il programma

1. **Compilare il programma:** Utilizzare un compilatore Root / C++ (con librerie root) per compilare il codice sorgente.
2. **Creare un file di dati:** Creare un file di testo chiamato `dati_da_interpolare.txt` con i dati da interpolare, nel formato descritto sopra (ogni riga x, y, sigma_y distanziati da un tab e con separatore decimale ".").
3. **Definire la funzione interpolante:** Definire, nell'apposita riga di codice, la funzione interpolante i dati utilizzando come variabile indipendente x[i_esimo_valore], mentre come parametri (in numero necessario non maggiore a 3) par1, par2 e par3.
4. **Definire il range dei parametri:** Definire, nelle apposite righe di codice, il range di valori in cui devono essere ricercati i parametri (definire numero_punti_iniziali=100 per un impegno computazionale medio).
5. **Eseguire il programma:** Eseguire il programma compilato.

## Requisiti

1. **Root:** Utilizzare Root per compilare il codice.

## Esempio di file di dati
160	0.281486702	0.010260399

161	0.293361922	0.010324526

162	0.325720326	0.010388654
163	0.34821413	0.010452781
164	0.370959261	0.010516909
165	0.41469023	0.010581036
166	0.458923855	0.010645164
167	0.503660134	0.010709291
168	0.55945482	0.010773419
169	0.63711499	0.010837546
170	0.737017637	0.010901674 
171	0.848795503	0.010965801
172	0.994251243	0.011029929
173	1.173950343	0.011094056
174	1.377525547	0.011158184
175	1.594358272	0.011222311
176	1.813578607	0.011286439
177	2.001822839	0.022701132
178	2.16970955	0.022829387
179	2.294367947	0.022957642
180	2.397663513	0.023085897
181	2.479219259	0.023214152
182	2.561528986	0.023342407
183	2.621596238	0.023470662
184	2.659044022	0.023598917
185	2.696743134	0.023727172
186	2.734693573	0.023855427
187	2.77289534	0.023983682
188	2.811348434	0.024111937
189	2.826302415	0.024240192
190	2.841256396	0.024368447

## Esempio di funzione interpolante i dati
return M_PI/2 - atan((par1 * par1 - pow(2 * M_PI * x[i_esimo_valore] * 1000, 2)) / (2 * M_PI * x[i_esimo_valore] * 1000 * par2 * 2));

## Esempio range parametri
double numero_punti_iniziali = 100;                // MODIFICABILE
double range_min1 = 1050000;                       //
double range_max1 = 1150000;                       //
double range_min2 = 20000;                         //
double range_max2 = 40000;                         //
double range_min3 = 0;                             //
double range_max3 = 1;                             //
double precisione_bisezione = 0.0000001;           //

## Note

* Il programma è stato progettato per interpolare funzioni con un massimo di tre parametri.
* Il file di dati deve essere chiamato `dati_da_interpolare.txt` e deve essere situato nella stessa directory del programma.
* Il programma utilizza la libreria ROOT per la generazione dei grafici.

## Autori

* [Nicola Lavarda]

## Licenza

Questo programma è rilasciato sotto la licenza [Apache License].

Apache License
Version 2.0, January 2004
http://www.apache.org/licenses/

TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

1. Definitions.
...

(Il testo completo della licenza può essere trovato su http://www.apache.org/licenses/LICENSE-2.0)


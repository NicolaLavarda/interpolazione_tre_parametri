# Interpolazione a più parametri

Questo programma implementa un algoritmo di interpolazione a più parametri (fino a 3). 

## Funzionalità

* **Lettura dati:** Il programma legge i dati da un file di testo, in cui ogni riga rappresenta un punto dati con le seguenti colonne:
    * `val_x`: valore della variabile indipendente
    * `val_y`: valore della variabile dipendente
    * `val_sigma_y`: errore sulla variabile dipendente
* **Determinazione del numero di parametri:** Il programma identifica automaticamente il numero di parametri necessari per l'interpolazione, analizzando la funzione interpolante.
* **Ricerca dei parametri:** Il programma utilizza una serie di cicli for e un rapido algoritmo di bisezione per trovare i parametri che minimizzano il chi-quadro.
* **Calcolo degli errori:** Il programma calcola gli errori sui parametri utilizzando la defnizione di standard error.
* **Grafici:** Il programma genera grafici che mostrano i dati, la funzione interpolante e il chi-quadro in funzione dei parametri.

## Come usare il programma

1. **Includere le librerie:** Utilizzando direttamente Root non è necessario, altrimenti includere tutte le librerie necessarie.
2. **Creare un file di dati:** Creare un file di testo chiamato `dati_da_interpolare.txt` con i dati da interpolare, nel formato descritto sopra (ogni riga `x`, `y`, `sigma_y` distanziati da un tab e con separatore decimale un punto `.`).
3. **Definire la funzione interpolante:** Definire, nell'apposita riga di codice, la funzione interpolante i dati utilizzando come variabile indipendente `x[i_esimo_valore]`, mentre come parametri (in numero necessario non maggiore a 3) `par1`, `par2` e `par3`.
4. **Definire il range dei parametri:** Definire, nelle apposite righe di codice, il range di valori in cui devono essere ricercati i parametri (definire `numero_punti_iniziali=100` per un impegno computazionale medio).
5. **Eseguire il programma:** Compilare ed seguire il programma con Root.

## Requisiti

1. **Root:** Utilizzare Root per compilare il codice.

## Esempio di file di dati
160   0.281486702   0.010260399<br>
161	0.293361922   0.010324526<br>
162	0.325720326   0.010388654<br>
163	0.348214130   0.010452781<br>
164	0.370959261   0.010516909<br>
165	0.414690230   0.010581036<br>
166   0.458923855   0.010645164<br>
167   0.503660134   0.010709291<br>
168	0.559454820   0.010773419<br>
169	0.637114990   0.010837546<br>
170	0.737017637   0.010901674<br>
171	0.848795503   0.010965801<br>
172	0.994251243   0.011029929<br>
173	1.173950343   0.011094056<br>
174	1.377525547   0.011158184<br>
175	1.594358272   0.011222311<br>
176	1.813578607   0.011286439<br>
177	2.001822839   0.022701132<br>
178	2.169709550   0.022829387<br>
179	2.294367947   0.022957642<br>
180	2.397663513   0.023085897<br>
181	2.479219259   0.023214152<br>
182	2.561528986   0.023342407<br>
183	2.621596238   0.023470662<br>
184	2.659044022   0.023598917<br>
185	2.696743134   0.023727172<br>
186	2.734693573   0.023855427<br>
187	2.772895340   0.023983682<br>
188	2.811348434   0.024111937<br>
189	2.826302415   0.024240192<br>
190	2.841256396   0.024368447<br>

## Esempio di funzione interpolante i dati
`return M_PI/2 - atan((par1 * par1 - pow(2 * M_PI * x[i_esimo_valore] * 1000, 2)) / (2 * M_PI * x[i_esimo_valore] * 1000 * par2 * 2));`

## Esempio range parametri
double numero_punti_iniziali = 100;                <br>
double range_min1 = 1050000;                       <br>
double range_max1 = 1150000;                       <br>
double range_min2 = 20000;                         <br>
double range_max2 = 40000;                         <br>
double range_min3 = 0;                             <br>
double range_max3 = 1;                             <br>
double precisione_bisezione = 0.0000001;           <br>

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


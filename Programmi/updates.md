# 5 Relazione


```
struct data_torque
{
    int data_file_number;  //numero progressivo
    double data_file_freq; //freq fornita
    vector<double> time;
    vector<double> a;
    vector<double> forcing;
};
```


## Lettura file esterni
- Salva per ogni torq in vettore di strutture come sopra riportato
- Cerca in modo approssimato la presenza degli zeri, calcola l'interolazione lineare e salva la radice
- Calcola tutti i periodi
- Calcolo omega in due modi, da media omega o media tmepi
- Fix intercette
- Metodo per trovare indice dei massimi
- Interpolazione con MV parabola sui 18 punti vicino ai picchi
- Stima di ampiezza max da coeff derivati da MV

## Stima degli zeri

## Stima delgi omega
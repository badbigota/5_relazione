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

## Stima degli zeri

## Stima delgi omega
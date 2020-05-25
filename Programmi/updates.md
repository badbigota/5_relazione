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
- Idem con root
- Idem con max ass
- Stima periodo con intercette di ampiezza
- stima di ampiezza media da tutti e 3 i metodi
- Plot lorenziana
- Stima di err a post su lorentz
- Check offset
- Check gaussiane ampiezze theta


- Decadimento, trova zeri
- Trova max
- Plot max  elinearizzzione
- Uso err di theta post di lorenz
- Fit stima di gamma

## Stima degli zeri

## Stima delgi omega
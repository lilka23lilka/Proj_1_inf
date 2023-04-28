
READ ME - Projekt 1 Informatyka geodezyjna II

1) Do czego służy program ?
Program został napisany, aby w szybki sposób przeprowadzić transformacje wspołrzędnych między układami. 

2) Jakie transformacje współrzędnych on objemuje?
Zawiera on cztery przekształcenia wspolrzednych:
- XYZ2BLH -> jest to transformacja obejmująca przekształcenie współrzędnych ortokartograficznych na układ współrzędnych geodezyjnych. Powszechnie nazwyna Algorytmem Hirvonena. Początkowymi danymi są X,Y,Z podanego punktu. Następnie poprzez obliczenia uzyskujemy ostetcznie szerokości gedozyjnej (fi), długości geodezyjnej (lambda), a także wysokości. Otrzymane wartości zostaną zapisane do nowo utworzonego pliku o rozszerzeniu txt. Uzyskane dane bedą przedstawione w odpowiednich jednostkach: fi i lamda w stopniach dziesietnych, natomiast wysokość w metrach. 
- BLH2XYZ -> jest do przekształcenie odwrotne do powyższego Algorytmu Hirvonena. Podając szrokość, długość geodezyjną w stopniach dziesiętnych oraz wysokość w metrach otrzymamy X,Y,Z w metrach. 
- XYZ2neu -> przekształca współrzędne XYZ na wektor przestrzenny neu
- BL21992 -> przekształca współrzędne geodezyjne na współrzędne X,Y w układzie 1992. Dane początkowe oraz końcowe powinny byc w metrach
- BL22000 -> polega na przekształceniu współrzędnych geodezyjnych na współrzędne X,Y w układzie 2000. Zarówno dane początkowe jak i końcowe powinny być w jednostce jaką jest metr.

3) Na jakich modelach elipsoid program umożliwia obliczenia?
- GRS80
- WGS84
- Krasowski

4) Na jakim systemie program został napisany oraz jakie są wymagania sprzętowe?
-System operacyjny Windows 10,
-64 - biotwy system operacyjny
-Python 3.10 - proponowane : Spyder (3.4.5) - wymagane działanie komendy 'input'
-Biblioteki użyte w programie: Math, Numpy, Argparse

5) Przykładowa transformacja XYZ2BLH:

Przykładowe dane początkowe jakie należy podać w pliku tekstowym:
```sh
1   3689561.567     1096783.987     5039874.085
3   3681034.736     1102865.847     5049573.678
4   3664987.476     1117895.076     5025767.983
6   3679837.285     1073684.373     5073443.437
```

WAZNE! Liczby zmiennoprzecinkowe wymagają użycia kropki, nie przecinka. Dane należy oddzielic od siebie tab.

Wywołanie:
Odpowiedzieć na podane komendy:
- podać nazwę pliku tekstowego z danymi w formacie podanym powyżej, bądź cieżkę do tego pliku,
- podać nazwę transformacji, którą chcemy uruchomić,
- podać nazwę elipsoidy, dla której mają byc wykonane obliczenia,

Przykładowe wywołanie:
- Dane.txt
- XYZ2BLH
- Krasowski

Po uruchomieniu programu zostanie utworzony raport z transformowanymi współrzędnymi
przykład otrzymanych wspolrzednych geodezyjnych - transformacja XYZ -> BLH :
(format : nr B L H)
```sh
1   54.2926     22.2912     16
2   53.2826     19.2929     256
3   53.2929     19.2785     85
4   51.2823     16.2685     164
5   52.2959     19.2816     248                                  
```
Po wpisaniu 'KONTYNUUJ', program zostanie ponownie uruchomiony w pętli gotowy do ponownych obliczeń. w przeciwnym przypadku zakonczy on dzialanie.

6) Błedy
po wpisaniu błędnych parametrów wyświetlany jest następujący komunikat:
  PODANY PLIK NIE ISTNIEJE.

po podaniu pliku z niepoprawnym formatem danych  wyświetlany jest następujący komunikat:  
  NIEPOPRAWNY FORMAT DANYCH W PLIKU.
  
wszelkie inne błędy na chwilę obecną nie zostały wykryte 
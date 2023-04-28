# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 17:05:32 2023

@author: lilja
"""
import numpy as np
from math import *
from argparse import ArgumentParser


class Transormacje_współrzędnych :
    
    def __init__(self):
        '''
        Parametry elipsoid:
            a - duża półoś elipsoidy 
            e^2 - mimośród^2
        podane w postaci słownika list w formacie:
        {'elipsoida':[a,e^2]}
        '''
        self.elipsoidy = {'GRS80':[ 6378137, 0.00669438002290],
                           'WGS84': [ 6378137,  0.00669437999014],
                           'Krasowski': [ 6378245, 0.00669342162297]}
    
    
    def odczyt_txt(self, txt):
    
        with open(txt,'r') as txt:
            wiersze = txt.readlines()
            wsp = []
            for linijka in wiersze:
                dane = linijka.split()
                a = [float(i) for i in dane]
                a[0] = int(a[0])
                wsp.append(a)
        return(wsp)

    
    def XYZ2BLH(self,txt,elipsoida) :
        '''
        Parameters: plik tekstowy ze wspołrzędnymi  XYZ  w metrach oraz parametry wybranej elipsoidy
        działanie: ALgorytm Hirvonena - przekształca współrzędne XYZ (ortokartezjańskie) na BLH( współrzędne geodezyjne + wysokosc). Funkcja zawiera iteracje poprzez którą otrzymujemy odpowiednie przybliżenie B. 
        Returns: plik tekstowy ze wspołrzędnymi  BLH, współrzędne geodezyjne zwracane są w stopniach dziesiętnych i przyblizeniu do 5 miejsca po przecinku, natomiasy wysokoć przedstawiona jest w metrach z dokladnoscia do 3 miejsc.
            
        '''
        a = self.elipsoidy[elipsoida][0]
        e2 = self.elipsoidy[elipsoida][1]
        poczat_dane = self.odczyt_txt(txt)
        ost_dane = []
        dane_raport=[]
        for element in poczat_dane:
            Pkt,X,Y,Z = element
            #zmienna do obliczen
            p = np.sqrt((X**2) + (Y**2))
            #fi
            B = np.arctan(Z/(p * (1-e2)))
            while True:
                N = a / np.sqrt(1- e2 * np.sin(B)**2)
                H = (p / np.cos(B)) - N
                Bp = B
                B = np.arctan(Z / (p * (1 - e2 * (N / (N + H)))))
                if np.abs(Bp - B) <( 0.000001/206265):
                    break
            #zaokraglenie H
            H = round(H,3)
            #lambda
            L = np.arctan2(Y,X)
            #zamiana na stopnie, minuty,sekundy fi
            B = B * 180 / pi
            # Bst=B
            # Bd = int(B)
            # Bm = int(60*(B-Bd))
            # Bs = (B-Bd -Bm/60)*3600
            # Bs = round(Bs,5)
            # B=(Bd,Bm,Bs)
            #zamiana na stopnie, minuty,sekundy fi
            L = L * 180 / pi
            # Lst =L
            # Ld = int(L)
            # Lm = int(60*(L-Ld))
            # Ls = (L-Ld-Lm/60)*3600
            # Ls = round(Ls,5)
            # L = (Ld,Lm,Ls)
            ost_dane.append([Pkt,B,L,H])
            #raport
        with open('Raport_transformacja_XYZ2BLH.txt', 'w') as plik:
            for el in ost_dane:
                plik.write('{:^10d} {:^20.5f} {:^20.5f} {:^20.3f}\n'.format(el[0], el[1], el[2], el[3]))
        return (ost_dane)    
    
    
    def BLH2XYZ(self,txt,elipsoida):
        '''
        Parameters: plik tekstowy ze wspołrzędnymi  BLH  w stopniach dziesietnych oraz parametry wybranej elipsoidy
        działanie: funkcja odwrotma do ALgorytmu Hirvonena - przekształca współrzędne BLH ( współrzędne geodezyjne + wysokosc) na XYZ( ortokartezjańskie). Wzory na XYZ zawierją promień wodzący oraz wektor normalny.
        Returns: plik tekstowy ze wspołrzędnymi  XYZ, współrzędne sa w metrach z dokladnoscia do 3 miejsc po przecinku.
        '''
        a = self.elipsoidy[elipsoida][0]
        e2 = self.elipsoidy[elipsoida][1]
        poczat_dane = self.odczyt_txt(txt)
        ost_dane = []
        for element in poczat_dane:
            Pkt,B,L,H = element
            N = a / np.sqrt(1- e2 * np.sin(B)**2)
            B = B * pi / 180
            L = L * pi / 180
            #uzyskanie XYZ i zaokraglenie
            X = (N + H) * np.cos(B) * np.cos(L)
            X=round(X,3)
            Y = (N + H) * np.cos(B) * np.sin(L)
            Y = round(Y,3)
            Z = (N * (1 - e2) + H) * np.sin(B)
            Z = round(Z,3)
            ost_dane.append([Pkt,X,Y,Z])
            #raport
        with open('Raport_transformacja_BLH2XYZ.txt', 'w') as plik:
            for el in ost_dane:
                plik.write('{:^10d} {:^20.3f} {:^20.3f} {:^20.3f}\n'.format(el[0], el[1], el[2], el[3]))
        return(ost_dane)      
           
            
    def XYZ2NEU(self,txt,elipsoida):
        '''
        Parameters: plik tekstowy ze wspołrzędnymi  XYZ  dla punktu poczatkowego oraz koncowego w metrach oraz parametry wybranej elipsoidy
        działanie: Zamiana wspolrzednych XYZ punktu poczatkowego na wspolrzedne geodezyjne BL. Nastpenie tworzone są 3 macierze, na których podstawie stworzono macierz obrotu NEU. Nastpenie otrzymana macierz transponowano, przemnozono przez wektor przestrzenny XYZ. Poprzez to otrzymano szukany wektor w ukladzie wspolrzednych topocentrycznych NEU. 
        Returns: plik tekstowy ze wartosciami wektora NEU.
        '''
        a = self.elipsoidy[elipsoida][0]
        e2 = self.elipsoidy[elipsoida][1]
        wsp = self.odczyt_txt(txt)
        wek_neu = []
        for a in wsp:
            nr,xp,yp,zp,xp,yp,zp = a
            p =np.sqrt(xp ** 2 + yp **2)
            B =np.arctan(zp /( p * (1 - e2)))
            while True:
                N =a / np.sqrt(1- e2 * np.sin(B)**2)
                h = (p / np.cos(B)) - N
                Bw = B
                B = np.arctan( zp / (p * (1 - e2 * (N / (N + h)))))
                if np.abs(Bw - B) < (  0.000001/206265):
                    brak
            L =np.arctan2(yp,xp)
            R = np.array([[-np.sin(B)*np.cos(L), -np.sin(L), np.cos(B)*np.cos(L)],
                          [-np.sin(B)*np.sin(L),  np.cos(L), np.cos(B)*np.sin(L)],
                          [np.cos(B),                   0,             np.sin(B)]])
            
            XYZ = np.array([[xk-xp],[yk-yp],[zk-zp]])
            neu = R.T @ XYZ
            b = [nr,neu[0][0], neu[1][0],neu[2][0]]
            wek_neu.append(b)
        with open('xyz_to_neu.txt','w') as plik:
            for c in wek_neu:
                plik.write('{:10} {:15.3f} {:15.3f} {:15.3f}\n'.format(c[0], c[1], c[2], c[3]))
        return (wek_neu)
           
           
    def BL22000(self,txt,elipsoida):
        '''
        Parameters: plik tekstowy ze wspołrzędnymi BLH  w stopniach dziesietnych oraz parametry wybranej elipsoidy
        działanie: Zamiana wspolrzednych BL na wspolrzedne XYZ 2000. Dla otrzymanej dlugosci geodezyjnej zostaje dopasowana jedna z 4 stref(5,6,7,8) oraz odpowiedni poludnik osiowy (15,18,21,24[stop]). Nastepnie poprzez liczne obliczenia dazymy do uzyskania wspolrzednych XY na plaszczyznie Gaussa-Krugera. Nastepnie otrzymane wartosci przeskalowywujemy (skala ukladu 2000 = 0.99923). Dodjac odpwiednia wartosc do Y, uwzgldniajaca numer przynaleznej strefy otrzymujemy ostateczne wspolrzednie XY w ukladzie 2000.
        Returns: plik tekstowy ze ze wspołrzędnymi XY w układzie 2000, współrzędne sa w metrach z dokladnoscia do 3 miejsc po przecinku.
        '''
        a = self.elipsoidy[elipsoida][0]
        e2 = self.elipsoidy[elipsoida][1]
        poczat_dane = self.odczyt_txt(txt)
        ost_dane = []
        for element in poczat_dane:
            Pkt,B,L = element
            B = B * pi / 180
            L = L * pi / 180
            #Dopasowanie strefy
            L_stref = 0
            n = 0
            if L > 13.5 * pi / 180 and L < 16.5 * pi / 180:
                L_stref = L_stref + (15 * pi / 180)
                n = n + 5
            if L > 16.5 * pi / 180 and L < 19.5 * pi / 180: 
                L_stref = L_stref + (18 * pi / 180)
                n = n + 6
            if L > 19.5 * pi / 180 and L < 22.5 * pi / 180: 
                L_stref = L_stref + (21 * pi / 180)
                n = n + 7
            if L > 22.5 * pi / 180 and L < 25.5 * pi / 180: 
                L_stref = L_stref + (24 * pi / 180)
                n = n + 8
            #dane posrednie
            b2 = (a ** 2) * (1 - e2)
            ep2 = (a ** 2 - b2) / b2
            dL = L - L_stref
            t = np.tan(B)
            n2 = ep2 * (np.cos(B) ** 2)
            N = a / np.sqrt((1- e2) * (np.sin(B)**2))
            #sigma
            A0 = 1 - (e2 / 4) - ((3 * e2 ** 2) / 64) - ((5 * e2 ** 3) / 256)
            A2 = (3 / 8) * (e2 + (e2 ** 2) / 4 + (15 * e2 ** 3) / 128)
            A4 = (15 / 256) * (e2 ** 2 + (3 * e2 ** 3) / 4)
            A6 = (35 * e2 ** 3) / 3072
            sig = a * ((A0 * B) - (A2 * np.sin(2 * B)) + (A4 * np.sin(4 * B)) - (A6 * np.sin(6 * B)))
            # Wspolrzedne Gaussa-Krugera
            Xgk = sig + ((dL ** 2 / 2) * N * np.sin(B) * np.cos(B) * (1 + (((dL ** 2)/12) * (np.cos(B) ** 2) * (5 - t **2 + 9 * n2 + 4 * n2 ** 2)) + (((dL ** 4) / 360) * (np.cos(B) ** 4 ) * (61 - 58 * (t ** 2) + t ** 4 + 270 * n2 - 330 * n2 * (t ** 2)))))
            Ygk = dL * N * np.cos(B) * (1 + (((dL ** 2)/6) * (np.cos(B) ** 2) * (1 - t ** 2 + n2)) + (((dL ** 4 ) / 120) * (np.cos(B) ** 4) * (5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2)))   
            # Przeskalowanie na uklad 2000
            X2000 = Xgk * 0.999923
            X2000 = round(X2000,3)
            Y2000 = Ygk * 0.999923 + n * 1000000 + 500000
            Y2000 = round(Y2000,3)
            ost_dane.append([Pkt,X2000,Y2000])
            #raport
        with open('Raport_transformacja_BL2XY2000.txt', 'w') as plik:
            for el in ost_dane:
                plik.write('{:^10d} {:^20.3f} {:^20.3f}\n'.format(el[0], el[1], el[2]))
        return(ost_dane)             
    
    def BL21992(self,txt,elipsoida):
        '''
        Parameters: plik tekstowy ze wspołrzędnymi BLH w stopniach dziesietnych oraz parametry wybranej elipsoidy
        działanie: Zamiana wspolrzednych BL na wspolrzedne XYZ 1992. Działa analogicznie jak transformajca do ukladu 2000, jednak w przypadku ukladu 1992 przyjmujemy wartosc poludnika osiowego rowna 19 stopni. Skala tego ukladu wynosi 1992. Rowniez dazymy poczatkowo do otrzymania wspolrzednych na plaszczyznie Gaussa-Krugera. Nastepnie poprzez przeskalowanie i dalsze obliczenia uzyskujemy wspolrzedne w ukladzie 1992.
        Returns: plik tekstowy ze ze wspołrzędnymi XY w układzie 1992, współrzędne sa w metrach z dokladnoscia do 3 miejsc po przecinku.
        '''
        a = self.elipsoidy[elipsoida][0]
        e2 = self.elipsoidy[elipsoida][1]
        wsp = self.odczyt_txt(txt)
        wsp_92 = []
        for a in wsp:
            nr, B, L = a
            Brad = B * pi/180
            Lrad = L* pi/180
            lam0 = 19* pi/180
            b2 = (a**2) * (1-e2)
            ep2 = (a**2-b2)/b2
            dL = L - lam0
            t = tan(B)
            n2 = ep2 * (cos(fi)**2)
            N = a / np.sqrt(1- e2 * np.sin(B)**2)
            #
            A0 = 1 - (e2 / 4) - ((3 * e2 ** 2)/ 64) - ((5 * e2 ** 3) / 256)
            A2 = (3 / 8) *(e2 + (e2 ** 2) / 4 + (15 * e2**3) /128)
            A4 = (15 / 256)*(e2 ** 2 + (3 * e2 **3) / 4)
            A6 = (35 * e2 ** 3) / 3072
            sigma = a * ((A0 * B) - (A2 * sin(2 *B)) + (A4* sin(4 * B)) - (A6 * sin(6 * B)))
            #
            xgk = sigma + ((dL **2 / 2) * N * sin(B) * cos(B) * (1+ (((dL**2)/12) *(cos(B) ** 2) * (5 - t **2 + 9 * n2 + 4 * n2 ** 2)) + (((dL ** 4) / 360) * (cos(B)**4) * (61 - 58 * (t ** 2) + t ** 4 + 270 * n2 - 330 * n2 * (t ** 2)))))
            ygk = dl * N * cos(B) * (1 + (((dL**2)/6) * (cos(B) ** 2) * (1 - t ** 2 + n2)) + (((dL ** 4 ) / 120) * (cos(B) ** 4) * (5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2)))   
            m = 0.9993 
            x92 = xgk * m - 5300000
            y92 = ygk * m +500000
            b = [nr,x92,y92]
            wsp_92.append(b)
        with open('BL_to_1992','w') as plik:
            for c in wsp_92:
                plik.write('{:10} {:15.3f} {:15.3f}\n'.format(c[0],c[1],c[2]))   
        return(wsp_92) 

if __name__ == '__main__':
    transformacje = {'XYZ2BLH':'XYZ2BLH', 'BLH2XYZ':'BLH2XYZ', 'XYZ2NEU':'XYZ2NEU', 'BL22000':'BL22000', 'BL21992':'BL21992'}
    elipsoidy = {'KRASOWSKI':'KRASOWSKI' , 'GRS80':'GRS80' , 'WGS84':'WGS84'}
    parser = ArgumentParser(description='Obliczanie wspolrzednych')
    parser.add_argument('-p', type=str, help='Podaj nazwę pliku z danymi z rozszerzeniem. Jes;o plik znajduje sie w innym folderze nalezy podac jego sciezke.')
    parser.add_argument('-t', type=str, help='Precyzuje nazwe wybranego typu transformacji (XYZ2BLH, BLH2XYZ, XYZ2NEU, BL22000, BL21992)')
    parser.add_argument('-e', type=str, help='Precyzuje nazwe modelu elipsoidy (WGS84/ GRS80/ KRASOWSKI)')
    args = parser.parse_args()
    KONTYNUUJ = "KONTYNUUJ"
    try:
        while KONTYNUUJ == "KONTYNUUJ" :
            if args.p==None:
                args.p = input(str('Wklej sciezke do pliku txt z danymi: '))
            if args.t==None:
                args.t = input(str('Nazwa transformacji:')).upper()
            if args.e==None:
                args.e = input(str('Model elipsoidy:')).upper()
            klasa = Transformacje_wspolrzednych()
            trans = transformacje[args.t]
            elip = elipsoidy[args.e]
            #dopasowanie argumentu z input'u do odpowiedniej transformacji z podanymi argumentami odnosnie pliku i elipsoidy
            if trans == 'XYZ2BLH':
                poczat_dane = klasa.XYZ2BLH(args.p, elip)
            if trans == 'BLH2XYZ':
                poczat_dane = klasa.BLH2XYZ(args.p, elip)
            if trans == 'XYZ2NEU':
                poczat_dane = klasa.XYZ2NEU(args.p, elip)
            if trans == 'BL22000':
                poczat_dane = klasa.BL22000(args.p, elip)
            if trans == 'BL21992':
                poczat_dane = klasa.BL21992(args.p, elip)
            KONTYNUUJ = input(str("Aby kontynuowac program wpisz KONTYNUUJ  - w przeciwnym przypadku program zakonczy dzialanie : ")).upper()
            args.e = None
            args.p= None
            args.t= None      
    # błedy
    except FileNotFoundError:
        print('PODANY PLIK NIE ISTNIEJE.')
    except KeyError:
        print('WSPISANE PARAMETRY SĄ NIEPRAWIDŁOWE .')
    except IndexError:
        print('NIEPOPRAWNY FORMAT DANYCH W PLIKU.')
    except ValueError:
        print('NIEPOPRAWNY FORMAT DANYCH W PLIKU.')
    finally:
        print('PROGRAM ZAKOŃCZYŁ PRACĘ.')            
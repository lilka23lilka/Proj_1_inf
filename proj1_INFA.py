# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 17:05:32 2023

@author: lilja
"""
import numpy as np
from math import *

class Transormacje_współrzędnych :
    
    def __init__(self,elipsoida):
        self.elipsoida = {'GRS80':[ 6378137, 0.00669438002290],
                           'WGS84': [ 6378137,  0.00669437999014],
                           'Krasowski': [ 6378245, 0.00669342162297]}
    
    
    def odczyt_txt(txt):
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
        a = self.elipsoida[elipsoida][0]
        e2 = self.elipsoida[elipsoida][1]
        poczat_dane = self.odczyt_txt(txt)
        ost_dane = []
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
            Bst=B
            Bd = int(B)
            Bm = int(60*(B-Bd))
            Bs = (B-Bd -Bm/60)*3600
            Bs = round(Bs,5)
            B=(Bd,Bm,Bs)
            #zamiana na stopnie, minuty,sekundy fi
            L = L * 180 / pi
            Lst =L
            Ld = int(L)
            Lm = int(60*(L-Ld))
            Ls = (L-Ld-Lm/60)*3600
            Ls = round(Ls,5)
            L = (Ld,Lm,Ls)
            ost_dane.append([Pkt,B,L,H])
            #raport
        with open('Raport_transformacja_XYZ2BLH.txt', 'w') as plik:
            plik.write('Otrzymane wspolrzedne geodezyjne - transformacja XYZ -> BLH \n')
            plik.write('{:^10s} {:^20s} {:^20s} {:^20s}\n'.format('Punkt','B ','L','H [m]'))
            for el in ost_dane:
                plik.write('{:^10d} {:^20.5f} {:^20.5f} {:^20.3f}\n'.format(el[0], el[1], el[2], el[3]))
        return (ost_dane)
    
    
    
    def BLH2XYZ(self,txt,elipsoida):
           a = self.elipsoida[elipsoida][0]
           e2 = self.elipsoida[elipsoida][1]
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
               ost_dane.append([Nr_pkt,X,Y,Z])
               #raport
           with open('Raport_transformacja_BLH2XYZ.txt', 'w') as plik:
               plik.write('Otrzymane wspolrzedne kartzejanskie - transformacja BLH -> XYZ \n')
               plik.write('{:^10s} {:^20s} {:^20s} {:^20s}\n'.format('Punkt','X [m]','Y [m]','Z [m]'))
               for el in ost_dane:
                   plik.write('{:^10d} {:^20.3f} {:^20.3f} {:^20.3f}\n'.format(el[0], el[1], el[2], el[3]))
           return(ost_dane)      
           
            
    def XYZ2neu(self,txt,elipsoida):
        a = self.elipsoida[elipsoida][0]
        e2 = self.elipsoida[elipsoida][1]
        wsp = odczyt_txt(txt)
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
           
           
    def BL2XY2000(self,txt,elipsoida):
           a = self.elipsoida[elipsoida][0]
           e2 = self.elipsoida[elipsoida][1]
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
               plik.write('Otrzymane wspolrzedne w ukladzie 2000 - transformacja BL -> XY 2000 \n')
               plik.write('{:^10s} {:^20s} {:^20s}\n'.format('Punkt','X [m]','Y [m]'))
               for el in ost_dane:
                   plik.write('{:^10d} {:^20.3f} {:^20.3f}\n'.format(el[0], el[1], el[2]))
           return(ost_dane)             
    
    def BL21992(self,txt,elipsoida):
        a = self.elipsoida[elipsoida][0]
        e2 = self.elipsoida[elipsoida][1]
        wsp = odczyt_txt(txt)
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
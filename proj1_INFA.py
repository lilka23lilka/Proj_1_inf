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
                   plik.write('{:^10s} {:^20s} {:^20s} {:^20s}\n'.format('Pkt','X [m]','Y [m]','Z [m]'))
                   for el in ost_dane:
                       plik.write('{:^10d} {:^20.3f} {:^20.3f} {:^20.3f}\n'.format(el[0], el[1], el[2], el[3]))
               return(ost_dane)
           
            
           #miejsce na neu
           
           
    
    
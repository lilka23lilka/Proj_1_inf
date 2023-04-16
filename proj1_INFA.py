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
        with open('Raport_transformacja_XYZ2BLH.txt', 'w') as plik:
            plik.write('{:^10s} {:^20s} {:^20s} {:^20s}\n'.format('Punkt','B','L','H'))
            for element in ost_dane:
                plik.write('{:^10} {:^20.5f} {:^20.5f} {:^20.3f}\n'.format(a[0], a[1], a[2], a[3]))
        return (ost_dane)
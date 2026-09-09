"""Intervalles fermés binary64 : opérations puis voisins vers l'extérieur.

Les entrées finies sont exactes. Un dépassement ou une division par un
intervalle contenant zéro refuse le calcul. Aucun mode global n'est changé.
Le compteur porte sur les opérations arithmétiques des extrémités, pas
sur les appels nextafter. Prototype de recherche, hors API publique.
"""
import math
import sys
import struct

from vinkulum._ports.inverse_selectionnee import BudgetInverseSelectionnee


class IntervallesBinaires:
    zero = (0.0, 0.0)
    un = (1.0, 1.0)

    def __init__(self, budget):
        if (sys.float_info.radix, sys.float_info.mant_dig,
                sys.float_info.max_exp) != (2, 53, 1024):
            raise ArithmeticError("binary64 IEEE requis")
        # Refuser FTZ/DAZ : l'encadrement au sous-flux nécessite les subnormaux.
        tiny = float.fromhex('0x1p-1074')
        if (struct.unpack('=Q',struct.pack('=d',tiny+tiny))[0] != 2
                or struct.unpack('=Q',struct.pack('=d',float.fromhex('0x1p-1022')*.5))[0] != 0x0008000000000000):
            raise ArithmeticError("arithmétique avec sous-flux graduel requise")
        self.budget = budget
        self.operations = 0

    def compter(self, n):
        if self.operations + n > self.budget:
            raise BudgetInverseSelectionnee("budget d'opérations binary64 dépassé",
                dict(phase="certification", operations_demandees=self.operations+n,
                     budget_operations=self.budget))
        self.operations += n

    @staticmethod
    def point(x):
        x = float(x)
        if not math.isfinite(x):
            raise ValueError("donnée finie requise")
        return x, x

    @staticmethod
    def enveloppe(lo, hi):
        lo = math.nextafter(lo, -math.inf)
        hi = math.nextafter(hi, math.inf)
        if not math.isfinite(lo) or not math.isfinite(hi):
            raise ArithmeticError("encadrement binary64 fini indisponible")
        return lo, hi

    def add(self, a, b):
        if a == self.zero:
            return b
        if b == self.zero:
            return a
        self.compter(2)
        return self.enveloppe(a[0]+b[0], a[1]+b[1])

    def sub(self, a, b):
        if b == self.zero:
            return a
        if a == self.zero:
            return self.neg(b)
        self.compter(2)
        return self.enveloppe(a[0]-b[1], a[1]-b[0])

    @staticmethod
    def neg(a):
        return -a[1], -a[0]

    def mul(self, a, b):
        if a == self.zero or b == self.zero:
            return self.zero
        if a == self.un:
            return b
        if b == self.un:
            return a
        if a[0] >= 0:
            if b[0] >= 0:
                lo, hi = a[0]*b[0], a[1]*b[1]
            elif b[1] <= 0:
                lo, hi = a[1]*b[0], a[0]*b[1]
            else:
                lo, hi = a[1]*b[0], a[1]*b[1]
        elif a[1] <= 0:
            if b[0] >= 0:
                lo, hi = a[0]*b[1], a[1]*b[0]
            elif b[1] <= 0:
                lo, hi = a[1]*b[1], a[0]*b[0]
            else:
                lo, hi = a[0]*b[1], a[0]*b[0]
        elif b[0] >= 0:
            lo, hi = a[0]*b[1], a[1]*b[1]
        elif b[1] <= 0:
            lo, hi = a[1]*b[0], a[0]*b[0]
        else:
            self.compter(2)
            lo = min(a[0]*b[1], a[1]*b[0])
            hi = max(a[0]*b[0], a[1]*b[1])
        self.compter(2)
        return self.enveloppe(lo, hi)

    def carre(self, a):
        if a == self.zero:
            return self.zero
        self.compter(2)
        if a[0] >= 0:
            lo, hi = a[0]*a[0], a[1]*a[1]
        elif a[1] <= 0:
            lo, hi = a[1]*a[1], a[0]*a[0]
        else:
            _, hi = self.enveloppe(0., max(a[0]*a[0], a[1]*a[1]))
            return 0., hi
        lo, hi = self.enveloppe(lo, hi)
        return max(0., lo), hi

    def div(self, a, b):
        if b[0] <= 0 <= b[1]:
            raise ArithmeticError("diviseur non séparé de zéro")
        if a == self.zero:
            return self.zero
        if b == self.un:
            return a
        self.compter(4)
        products = (a[0]/b[0], a[0]/b[1], a[1]/b[0], a[1]/b[1])
        return self.enveloppe(min(products), max(products))

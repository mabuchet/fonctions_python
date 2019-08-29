# -*- coding: utf-8 -*-

""" 
Auteur : Marc-Antoine BUCHET
Date : 28/06/2019

Script de développement d'une fonction appelée format_scientifique et qui a pour
but d'écrire un nombre x affecté de son incertitude dx en écriture scientifique
en respectant le bon nombre de chiffres significatifs sur x.

Si x est plus grand que son incertitude, alors on écrit x en écriture 
scientifique et on rajoute le nombre de zéro qu'il faut sur l'incertitude.
"""

import numpy as np
import random

################################################################################
# Fonctions : 
################################################################################
def puissance_de_dix(x):
    """Donne la puissance de 10 du réel x"""
    return int(np.floor(np.log10(abs(x))))

def format_scientifique(x,dx,CS=2,fmt='std'):
    """Utilise la valeur de x et de son incertitude dx pour générer une chaine
    de caractères contenant x et dx avec le bon nombre de chiffres significatifs
    pour une notation scientifique.
    
    Dans tous les cas, on ne garde de dx que le nombre de chiffres significatifs
    choisi par CS (2 par défaut).
    
    Si x>dx (cas le plus usuel) alors on note x en écriture scientifique et 
    c'est dx qui impose le nombre de décimales à x. On a alors deux choix de 
    formats : 
    - le choix standard ('std') : on note l'incertitude avec le même nombre de
      chiffres significatifs que x et la même puissance, la puissance noté une
      fois à la fin en utilisant des parenthèses autout de x et dx.
      ex : x = 10,777777 et dx=0,33 alors on note x = (1,0778 \pm 0,033)x10^1
    - le choix métrologique ('NIST') : on note l'incertitude entre parenthèses
      à la fin des décimales de x. Cette notation n'est adaptée que si dx
      est suffisamment petit devant x.
     
    Si x<=dx alors on note dx en écriture scietifique et on ne garde de x
    que le nombre de chiffres significatifs correspondant. Ici, il n'y a pas 
    de choix de format."""
    # On récupère les puissances de x et dx :
    p_x = puissance_de_dix(x)
    p_dx = puissance_de_dix(dx)
    # On gère la notation selon le cas x>dx ou x<=dx
    if x > dx : 
        # Le nombre de décimales dans l'écriture de x est ici imposé par la  
        # valeur de dx et le nombre de chiffres significatifs choisi :
        nbre_decimales = p_x-p_dx+CS-1
        # On récupère la mantisse de x :
        m_x=x/10**p_x
        if fmt in ['NIST']:
            # On récupère la mantisse de dx :
            m_dx = dx/10**p_dx
            # On récupère la valeur de l'incertitude sur les derniers chiffres 
            # de x avec le nombre de chiffres significatifs voulu (on ajoute 0.5 
            # pour arrondir à l'entier le plus proche et pas simplement à la
            # partie entière) :
            incertitude = int(m_dx*10**(CS-1)+0.5)
            # On formate la chaine de caractère :
            s = '{0:.'+'{}'.format(nbre_decimales)+'f}'
            s = s.format(m_x)
            s += '({0:2d})'.format(incertitude)
            s += 'e{0:+03d}'.format(p_x)
        elif fmt in ['std'] :
            # On récupère la "mantisse" de dx (à la même puissance que x) :
            m_dx = dx/10**p_x
            s = '({0:.'+'{}'.format(nbre_decimales)+'f} \pm {1:.'
            s += '{}'.format(nbre_decimales)
            s += 'f})e{2:+03d}'
            s = s.format(m_x,m_dx,p_x)
        else : 
            raise ValueError("Format choisi inconnu")
    elif x<=dx :
        # Le nombre de décimales dans l'écriture de dx est ici imposé par le
        # nombre de chiffres significatifs choisi :
        nbre_decimales = CS-1
        # On récupère la mantisse de dx :
        m_dx=dx/10**p_dx
        # On ramène la "mantisse" de x à la même puissance que dx :
        m_x = x/10**p_dx
        s = '({0:.'+'{}'.format(nbre_decimales)+'f} \pm {1:.'
        s += '{}'.format(nbre_decimales)
        s += 'f})e{2:+03d}'
        s = s.format(m_x,m_dx,p_dx)
    return s

################################################################################
# Tests : 
################################################################################
exit_values = ['Non','N','No','','0']
print("Tests des fonctions 'puissance_de_dix' et 'format_scientifique'.")
print("les valeurs pour quitter le test sont :")
for elmt in exit_values : print(elmt)

# Tests de puissance_de_10 :
print('\n\nTests de puissance_de_dix :\n')
b = True
while b :
    try : 
        x=float(input('Entrez un nombre decimal :'))
    except ValueError :
        print("Attention, ce n'est pas un nombre décimal !") 
        continue
    print(x)
    print('puissance de x :',puissance_de_dix(x))
    val=input('On continue ?')
    if val in exit_values : b=False

# Tests de format_scientifique :
print('\n\nTests de format_scientifique :\n')
b = True
while b :
    try : 
        x=float(input('Entrez un nombre decimal pour x :'))
        dx=float(input('Entrez un nombre decimal pour dx :'))
    except ValueError :
        print("Attention, ce n'est pas un nombre décimal !") 
        continue
    print(x)
    print(dx)
    print(format_scientifique(x,dx))
    print(format_scientifique(x,dx,fmt='NIST'))
    val=input('On continue ?')
    if val in exit_values : b=False

print('\nIci, des nombres sont générés aléatoirement entre 1e-8 et 1e9 :\n')
b = True
while b :
    print('\n')
    m_x = random.uniform(1,10)
    p_x = random.randint(-8,8)
    x = m_x*10**p_x
    m_dx = random.uniform(1,10)
    p_dx = random.randint(-8,8)
    dx = m_dx*10**p_dx
    print(x)
    print(dx)
    if x<dx : print('Attention, x<dx')
    print(format_scientifique(x,dx))
    print(format_scientifique(x,dx,fmt='NIST'))
    val=input('On continue ?')
    if val in exit_values : b=False
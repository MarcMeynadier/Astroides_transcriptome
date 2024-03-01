

#DDCC77 / #003272 / #B090DA / #88CCEE / #CC6677

import matplotlib.pyplot as plt
import numpy as np


def barplot(title,expList,valuesPos,valuesNeg):
    # Générer des données aléatoires pour les barres
    #values = np.random.randint(1, 10, size=(3, 5))
    print(valuesPos)

    # Couleurs de la palette
    palette = ['#DDCC77', '#003272', '#CC6677', '#88CCEE', '#B090DA']

    # Coordonnées x pour les barres
    x = np.arange(5)

    # Création du graphique
    for i in range(len(expList)):
        for j in range(5):
            plt.bar(x[j] + i * 6, valuesPos[i][j], color=palette[j], label=f'Barre {j+1}',edgecolor='black')
            plt.bar(x[j] + i * 6, -valuesNeg[i][j], color=palette[j],edgecolor='black')

    # Affichage de la ligne à la valeur 0
    plt.axhline(y=0, color='black', linewidth=0.5)

    # Ajouter des lignes verticales pour séparer les sous-parties
    for i in range(1, 3):
        plt.axvline(x=(i * 6) - 1, color='black', linestyle='--')  # Ajustement de la position

    # Ajouter les noms "exp 1", "exp 2" et "exp 3" sous l'abscisse
    plt.text(2.5, -11,expList[0], ha='center', fontsize=12)
    plt.text(8.5, -11,expList[1], ha='center', fontsize=12)
    plt.text(14.5, -11,expList[2], ha='center', fontsize=12)
    if len(expList) == 4:
        plt.text(20.5, -11,expList[2], ha='center', fontsize=12)

    # Personnalisation du graphique
    #plt.xlabel('Barres')
    plt.ylabel('Genes number - Expression relative')
    plt.title(title)

    # Utiliser une liste d'étiquettes d'axe x pour toutes les barres dans le graphique
    total_labels = [str(i) for i in range(1, 16)]
    plt.xticks(np.arange(15), total_labels)

    # Affichage de la légende uniquement pour la première sous-partie
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[:5], labels[:5])

    # Affichage du graphique
    plt.show()


# Calcif, Metabo, Immune, Homeostasis, Other
valuesPos = [[0,0,0,0,0],[2,0,0.5,0,1],[3,0,2.5,0,0]]
valuesNeg = [[1,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]] 
expList = ["sp_sp_bck","sp_sp_tro","gm_sp_trt"]
title = "May 2018 - SP"
barplot(title,expList,valuesPos,valuesNeg)

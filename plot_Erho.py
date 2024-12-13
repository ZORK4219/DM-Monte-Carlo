import pandas as pd
import matplotlib.pyplot as plt

def graph_data_multiple(file_path1, file_path2):
    # Fonction pour lire et préparer les données
    def read_data(file_path):
        data = pd.read_csv(file_path, sep="|", skipinitialspace=True)
        data.columns = [col.strip() for col in data.columns]
        data = data.dropna()
        data = data.apply(pd.to_numeric, errors='coerce')
        return data

    # Lire les deux fichiers
    data1 = read_data(file_path1)
    data2 = read_data(file_path2)

    # Vérification des données
    if data1.isnull().values.any() or data2.isnull().values.any():
        print("Erreur : certaines données dans l'un des fichiers contiennent des valeurs invalides.")
        return

    # Tracer les données du premier fichier
    plt.figure(figsize=(10, 6))

    plt.scatter(data1["Paramètre"], data1["Fonction"], label="E($\\rho$) simulation", linestyle='-', color="navy", s=4)
    plt.axhline(y=0, color="darkorange", linestyle="--", label="E($\\rho$) gaz parfait")
    #plt.scatter(data2["Paramètre"], data1["Fonction"], label="E($\\rho$) 5K", linestyle='-', color="orange", s=4)
    # Ajouter les détails au graphe
    plt.title("Énergie en fonction de la densité")
    plt.xlabel("$\\rho$")
    plt.ylabel("E($\\rho$)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Afficher le graphe
    plt.show()

# Exemple d'appel de la fonction
# Remplacez "data1.txt" et "data2.txt" par les chemins réels vers vos fichiers texte
graph_data_multiple("/home/adrien34490/WorkDir/Codes/C/Monte Carlo/DM-Monte-Carlo/E_rho/Energie_rho.txt",
                    "/home/adrien34490/WorkDir/Codes/C/Monte Carlo/DM-Monte-Carlo/E_rho/Energie_rho.txt")

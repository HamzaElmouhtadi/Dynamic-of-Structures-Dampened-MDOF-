import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig

# Fonction Calc prenant en entrée Nx (nombre de degrés de liberté), mx (matrice de masse), kx (matrice de raideur), fx (vecteur de force), bx (matrice d'amortissement)
def Calc(Nx, mx, kx, fx, bx):
    N = Nx  # Assignation du nombre de degrés de liberté
    M = mx  # Assignation de la matrice de masse
    K = kx  # Assignation de la matrice de raideur
    Force = fx  # Assignation du vecteur de force
    B = bx  # Assignation de la matrice d'amortissement

    # Affichage des matrices de raideur et de masse pour vérification
    print(K)
    print(M)
    
    # Calcul des valeurs propres et des vecteurs propres de l'équation aux valeurs propres (inv(M) @ K)
    valeurs_propres, Modes = eig(K,M)
    
    # Affichage des modes propres pour vérification
    print(Modes)
    
    # Calcul des fréquences propres à partir des valeurs propres
    freq_libre = np.sqrt(valeurs_propres) / (2 * np.pi)
    
    # Visualisation des modes propres
    for i in range(N):
        plt.subplot(N, 1, i+1)
        plt.plot(range(N), Modes[:, i], color='blue', linestyle='-', marker='o')
        plt.title('mode {} : {:.2f} Hz'.format(i+1, np.real(freq_libre[i])))
    plt.show()

    # Normalisation des modes propres
    mui = np.diag(Modes.T @ M @ Modes)
    print("mui", mui)
    Modes = (1 / np.sqrt(mui)) * Modes
    
    # Calcul des fréquences naturelles non amorties
    wi = np.sqrt(np.diag(Modes.T @ K @ Modes))
    print("wi", wi)
    
    # Calcul des coefficients d'amortissement
    xi = (1 / (2 * wi)) * np.diag(Modes.T @ B @ Modes)
    print("xi", xi)
    
    # Calcul des fréquences amorties
    wd = wi * (np.sqrt(1 - xi ** 2))
    
    # Calcul du vecteur de forces modales
    F = Modes.T @ Force

    # Réponse impulsionnelle
    t = np.linspace(0, 20, 1000)  # Temps de simulation
    u = np.zeros(len(t))  # Initialisation de la réponse en déplacement
    for i in range(N):
        u += Modes[0, i] * (F[i] * np.exp(-xi[i] * wi[i] * t) / (mui[i] * wd[i])) * np.sin(wd[i] * t)
    
    # Affichage de la réponse impulsionnelle
    plt.figure()
    plt.plot(t, u)
    plt.xlabel('(sec)')
    plt.ylabel('(m)')
    plt.title('Reponse d\'une impulsion - masse 1')
    freqstart = float(input("Definir la plage frequencielle : debut = ")) 
    freqend= float(input("fin = "))
    # Diagramme de Bode
    w = np.linspace(freqstart*(2*np.pi), freqend*(2*np.pi), 100000)  # Plage de fréquence
    U = np.zeros(len(w), dtype=complex)  # Initialisation de la réponse en fréquence
    for i in range(N):
        U += (F[i] * Modes[0, i]) / (wi[i] ** 2 - w ** 2 + 2j * xi[i] * wi[i] * w)
    
    # Affichage du diagramme de Bode
    plt.figure()
    f = w / (2 * np.pi)
    plt.plot(f, 20 * np.log10(np.abs(U)), 'k', linewidth=2)
    plt.xlabel('F(rad/s)')
    plt.ylabel('Amplitude (dB)')
    plt.title('Diagramme bode déplacement masse 1')
    
    # Diagramme de Nyquist
    realU = np.real(U)
    imageU = np.imag(U)
    plt.figure()
    plt.plot(realU, imageU)
    plt.xlabel('Partie réelle')
    plt.ylabel('Partie imaginaire')
    plt.title('Diagramme nyquist déplacement masse 1')
    # Diagramme de phase
    phase = np.angle(U)
    phase_multiple_pi = phase / np.pi  # Conversion des angles en multiples de pi
    plt.figure()
    plt.plot(f, phase_multiple_pi)
    plt.xlabel('Fréquence')
    plt.ylabel('Angle (π)')
    plt.yticks(np.arange(-2, 3, 1), ['-2π', '-π', '0', 'π', '2π'])
    plt.title('Diagramme phase déplacement masse 1')
    plt.show()

# Données d'exemple
#N = 3  # Nombre de degrés de liberté
#M = np.array(np.zeros((N, N), dtype=float))  # Matrice de masse initialisée à zéro
#K = np.array(np.zeros((N, N), dtype=float))  # Matrice de raideur initialisée à zéro
#B = np.array(np.zeros((N, N), dtype=float))  # Matrice d'amortissement initialisée à zéro
#F = np.array(np.zeros((N, 1), dtype=float))  # Vecteur de force initialisé à zéro

# Assignation des valeurs de la matrice de masse, de raideur, d'amortissement et du vecteur de force
#M = np.array([[6.0000, 0, 0], [0, 0.5000, 0], [0, 0, 0.0800]])
#K = np.array([[1013000, -13000, 0], [-13000, 44000, -11000], [0, -11000, 11500]])
#B = np.array([[0.4, -0.4, 0], [-0.4, 1.1, 0], [0, 0, 0.7]])
#F = np.array([0, 0, 1000])

# Appel de la fonction Calc avec les valeurs assignées
#Calc(3, M, K, F, B)

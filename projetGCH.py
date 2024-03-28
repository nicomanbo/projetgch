import numpy as np

class ailette:
    Tw = 100
    T_inf = 0
    k = 10
    L = 0.2
    R = 0.2

def TempAilette(ailette,nz,nr, h):
    """ Fonction assemblant la matrice A et le vecteur b

    Entrées:
        - ailette : Classe contenant les paramètres de l'ailette
        - Z : Bornes du domaine en z, Z = [z_min, z_max]
        - R : Bornes du domaine en r, R = [r_min, r_max]
        - nz : Discrétisation de l'espace en z (nombre de points)
        - nr : Discrétisation de l'espace en r (nombre de points)
        - h : Coefficient de convection
    Sorties (dans l'ordre énuméré ci-bas):
        - A : Matrice (array)
        - b : Vecteur (array)
    """
    Z = [0, ailette.L]
    R = [0, ailette.R]
    vr = np.linspace(R[0], R[1], nr)
    dz = (Z[1] - Z[0]) / (nz-1)
    dr = (R[1] - R[0]) / (nr-1)
    n = nr*nz
    b = np.zeros(n)
    A = np.zeros((n, n))
    for i in range(1, nz-1):
        for j in range(1, nr-1):
            p = j + i * nr
            A[p,p] = (-2/dr**2) + (-2/dz**2)
            A[p, p+nr] = 1/dz**2
            A[p, p-nr] = 1/dz**2
            A[p, p+1] = (1/(2*dr*vr[j]) + 1/dr**2)
            A[p, p-1] = (-1/(2*dr*vr[j]) + 1/dr**2)
    
    for i in range(nz):
        for j in range(nr):
            p = i + j * nr
            au_extremite = (p + 1) % nr == 0
            debut_ailette = j == 0
            fin_ailette = p >= (nr * (nz-1))
            au_centre = i == 0
            if au_extremite:
                A[p][p] = -3
                A[p][p-1] = 4
                A[p][p-2] = -1
                b[p] = h * ailette.T_inf
            if au_centre:
                A[p][p] = (h + ((3 * ailette.k)/(2 * dr)))
                A[p][p-1] = (-2 * ailette.k) / dr
                A[p][p-2] = -ailette.k / (2 * dr)
            if debut_ailette:
                A[p][p] = 1 
                b[p] = ailette.Tw
            if fin_ailette:
                A[p][p] = 1
                b[p] = ailette.T_inf
    return A, b

print(TempAilette(ailette, 3, 3, 10 ))
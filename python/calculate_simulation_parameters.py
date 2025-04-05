import numpy as np
from dataclasses import dataclass

@dataclass
class Molecule:
    d_ref: float
    T_ref: float
    m: float


def mean_free_path_VHS(mol: Molecule, n: float, T: float, nu:float):
    return 1/(np.sqrt(2)*np.pi*mol.d_ref**2*n*(mol.T_ref/T)**nu)

def density_from_lambda_VHS(mol, lambda_, T, nu):
    return 1/(lambda_*np.sqrt(2)*np.pi*mol.d_ref**2*(mol.T_ref/T)**nu)



if __name__ == "__main__":

    N2 = Molecule(m=4.6518e-27, d_ref=4.17e-10, T_ref=273.15)

    n = 1e20
    T = 1
    nu = .24

    print(mean_free_path_VHS(N2, n, T, nu))

    print(density_from_lambda_VHS(N2, 10, T=300, nu=nu))
    print(density_from_lambda_VHS(N2, 1, T=100, nu=nu))

    print("nu", mean_free_path_VHS(N2, n, T, nu)*100)
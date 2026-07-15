import math
import numpy as np
import sympy as sp
from scipy.special import eval_hermite

import config as cfg
from numerical_methods import solve_hamiltonian, INTEGRATION_DISPATCH

# Global variables to cache the compiled SymPy functions invisibly
_CALC_V_DIAG = None
_CALC_V_OFF_DIAG = None

def _compile_analytical_matrix_elements():
    """Compiles the SymPy matrix elements exactly once and caches them."""
    global _CALC_V_DIAG, _CALC_V_OFF_DIAG
    
    # If they are already compiled, skip the math
    if _CALC_V_DIAG is not None and _CALC_V_OFF_DIAG is not None:
        return
        
    print("Initializing SymPy Engine: Compiling analytical matrix elements... (This takes ~10s once)")
    
    # 1. Define Symbolic Variables (Using 'k' as the index so it doesn't collide with mass 'm')
    u, L_sym, mass_sym, omega_sym = sp.symbols('u L m omega', real=True, positive=True)
    n_sym, k_sym = sp.symbols('n k', integer=True, positive=True)
    
    prefactor = (mass_sym * omega_sym**2) / L_sym
    position_squared = (u - L_sym/2)**2
    
    integrand_diag = prefactor * position_squared * sp.sin(n_sym * sp.pi * u / L_sym)**2
    V_diag_expr = sp.integrate(integrand_diag, (u, 0, L_sym)).simplify()
    
    integrand_off = prefactor * position_squared * sp.sin(n_sym * sp.pi * u / L_sym) * sp.sin(k_sym * sp.pi * u / L_sym)
    V_off_expr = sp.integrate(integrand_off, (u, 0, L_sym)).simplify()
    
    # 2. Compile with perfectly unique variable names
    _CALC_V_DIAG = sp.lambdify((n_sym, L_sym, mass_sym, omega_sym), V_diag_expr, 'numpy')
    _CALC_V_OFF_DIAG = sp.lambdify((n_sym, k_sym, L_sym, mass_sym, omega_sym), V_off_expr, 'numpy')
    
    print("Compilation complete. Engine ready!")


def basis_function(x, n, L) -> float:
    """Evaluates the n-th infinite square well basis function at position x."""
    normalization = np.sqrt(2.0 / L)
    argument = (n * np.pi * (x + L / 2.0)) / L
    return normalization * np.sin(argument)

def integrand_V(x, n, k, L, mass, omega) -> float:
    """The potential energy integrand: chi_n(x) * V(x) * chi_k(x)."""
    chi_n = basis_function(x, n, L)
    chi_k = basis_function(x, k, L)
    
    V_x = 0.5 * mass * ((omega * x)**2)
    return chi_n * V_x * chi_k

def build_hamiltonian_matrix(N_BASIS, L, mass, omega, integration_method='trapezoidal', N_GRID=1000):
    """Builds the N_BASIS x N_BASIS Hamiltonian matrix"""
    H = np.zeros((N_BASIS, N_BASIS))
    
    for i in range(N_BASIS):
        for j in range(i, N_BASIS):
            # Shift indices to match basis functions (k instead of m)
            n = i + 1
            k = j + 1

            # --- Kinetic Energy (Exact analytical result) ---
            T = 0.0
            if n == k:
                momentum = (n * np.pi) / L
                T = ((cfg.HBAR * momentum)**2) / (2.0 * mass)

            # --- Potential Energy (Analytical or Numerical) ---
            V = 0.0
            if integration_method == 'analytical':
                if _CALC_V_DIAG is None or _CALC_V_OFF_DIAG is None:
                    _compile_analytical_matrix_elements()
                
                if n == k:
                    V = _CALC_V_DIAG(n, L, mass, omega)
                else:
                    V = _CALC_V_OFF_DIAG(n, k, L, mass, omega)

            elif integration_method in INTEGRATION_DISPATCH:
                integration_method_func = INTEGRATION_DISPATCH[integration_method]
                V = integration_method_func(integrand_V, -L/2, L/2, N_GRID, n, k, L, mass, omega)

            else:
                raise ValueError("Invalid integration method selected!")
            
            H[i, j] = T + V
            if i != j:
                H[j, i] = H[i, j] 
            
    return H

def exact_qho_wavefunction(x, n, mass=cfg.MASS, omega=cfg.OMEGA, hbar=cfg.HBAR):
    """Evaluates the exact analytical QHO wavefunction for state n at position(s) x."""
    xi = np.sqrt((mass * omega) / hbar) * x
    prefactor = (mass * omega / (np.pi * hbar)) ** 0.25
    
    # Use float to avoid 64-bit integer overflow for large n
    N_n = prefactor / np.sqrt(float(2**n) * float(math.factorial(n))) 
    gaussian = np.exp(-(xi**2) / 2.0)
    hermite = eval_hermite(n, xi)
    
    return N_n * gaussian * hermite

def reconstruct_numerical_wavefunction(x_array, L, N_basis, eigenvector):
    """Reconstructs the spatial wavefunction by multiplying the expansion coefficients by basis functions."""
    psi_num = np.zeros_like(x_array)
    for i in range(N_basis):
        n = i + 1
        psi_num += eigenvector[i] * basis_function(x_array, n, L)
    return psi_num

def calculate_qho_energies_and_states(L_val, N_BASIS_val, N_STATES_val,  mass=cfg.MASS, omega=cfg.OMEGA):
    """Master pipeline: Builds the matrix and solves it for given parameters."""
    H_matrix = build_hamiltonian_matrix(
        N_BASIS=N_BASIS_val, 
        L=L_val, 
        mass=mass, 
        omega=omega, 
        integration_method=cfg.integration_method, 
        N_GRID=cfg.N_GRID
    )
    
    energies, states = solve_hamiltonian(
        H=H_matrix, 
        N_STATES=N_STATES_val, 
        eigen_method=cfg.eigenvalue_method
    )
    
    return energies, states
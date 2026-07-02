import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh

def integral_solve_using_trapezoidal_rule(f, a, b, N, *args) -> float:
    """
    Solves the definite integral of a function using the trapezoidal rule.

    Parameters
    ----------
    f : callable
        The objective function to be integrated. Must accept the integration variable
        as its first argument, followed by any additional arguments (*args).
    a : float
        The lower limit of integration.
    b : float
        The upper limit of integration.
    N : int
        The number of sub-intervals to divide the integration range into.
    *args : tuple, optional
        Additional arguments to pass to the objective function `f`.

    Returns
    -------
    float
        The estimated numerical value of the integral. The error scales as O(h^2).
    """
    h = (b - a) / N
    s = f(a, *args) / 2.0
    
    for i in range(1, N):
        s += f(a + h * i, *args)
        
    s += f(b, *args) / 2.0
    s *= h
    
    return s

def integral_solve_using_romberg_method(f, a, b, N, *args) -> float:
    """
    Solves the definite integral of a function using Romberg's method.
    
    This method extrapolates a more accurate estimate from two trapezoidal 
    approximations, significantly reducing the approximation error.

    Parameters
    ----------
    f : callable
        The objective function to be integrated. Must accept the integration variable
        as its first argument, followed by any additional arguments (*args).
    a : float
        The lower limit of integration.
    b : float
        The upper limit of integration.
    N : int
        The initial number of sub-intervals. The method will also evaluate 
        the function at 2*N intervals internally.
    *args : tuple, optional
        Additional arguments to pass to the objective function `f`.

    Returns
    -------
    float
        The estimated numerical value of the integral. The error scales as O(h^4).
    """
    # Calculate S_n using the trapezoidal rule
    s_n = integral_solve_using_trapezoidal_rule(f, a, b, N, *args)

    # Calculate S_2n using the trapezoidal rule with double the intervals
    s_2n = integral_solve_using_trapezoidal_rule(f, a, b, 2 * N, *args)

    # Romberg extrapolation formula
    solution = (4.0 / 3.0) * s_2n - (1.0 / 3.0) * s_n
    
    return solution

def integral_solve_using_simpsons_rule(f, a, b, N, *args) -> float:
    """
    Solves the definite integral of a function using Simpson's 1/3 rule.

    Parameters
    ----------
    f : callable
        The objective function to be integrated. 
    a : float
        The lower limit of integration.
    b : float
        The upper limit of integration.
    N : int
        The number of sub-intervals. If an odd integer is provided, 
        it will be automatically incremented to the next even integer.
    *args : tuple, optional
        Additional arguments to pass to the objective function `f`.

    Returns
    -------
    float
        The estimated numerical value of the integral. The error scales as O(h^4).
    """
    # Force N to be an integer first, just in case a float was passed
    N = int(N)
    
    # Automatically adjust odd numbers to the nearest even number (rounding up)
    if N % 2 != 0:
        N += 1
        print(f"Warning: Simpson's rule requires an even N. Automatically adjusted N to {N}.")

    h = (b - a) / N
    s = f(a, *args) + f(b, *args)

    # Sum the odd-indexed terms (multiplied by 4)
    for i in range(1, N, 2):
        s += 4.0 * f(a + i * h, *args)
        
    # Sum the even-indexed terms (multiplied by 2)
    for i in range(2, N - 1, 2):
        s += 2.0 * f(a + i * h, *args)

    return (h / 3.0) * s

def scipy_quad_wrapper(f, a, b, N_GRID, *args):
    """
    Wraps SciPy's quad function to perfectly match our custom integration interface.
    (We accept N_GRID to satisfy the interface, but ignore it internally).
    """
    val, _ = quad(f, a, b, args=args, limit=100)
    return val

def custom_qr_algorithm(H, max_iter=2000, tol=1e-10):
    """
    Diagonalizes a symmetric matrix using the iterative QR algorithm.
    
    Parameters
    ----------
    H : np.ndarray
        The symmetric Hamiltonian matrix.
    max_iter : int
        The maximum number of QR iterations.
    tol : float
        The tolerance for convergence (off-diagonal elements near zero).
        
    Returns
    -------
    energies : np.ndarray
        1D array of eigenvalues.
    eigenvectors : np.ndarray
        2D array of eigenvectors (columns).
    """
    # Create a copy of H to avoid modifying the original matrix
    H_k = np.copy(H)
    n = H_k.shape[0]
    
    # Initialize the eigenvector matrix as the Identity matrix
    V_k = np.eye(n)
    
    for i in range(max_iter):
        # 1. QR Decomposition
        Q, R = np.linalg.qr(H_k)
        
        # 2. Reverse Multiplication (Similarity Transformation)
        H_k = R @ Q
        
        # 3. Accumulate the eigenvectors
        V_k = V_k @ Q
        
        # 4. Convergence Check (Are the off-diagonal elements practically zero?)
        # For a symmetric matrix, we can just check the lower triangle (excluding diagonal)
        off_diagonal_sum = np.sum(np.abs(np.tril(H_k, k=-1)))
        
        if off_diagonal_sum < tol:
            print(f"QR Algorithm converged in {i} iterations.")
            break
    else:
        print(f"Warning: QR Algorithm did not fully converge within {max_iter} iterations.")

    # The eigenvalues are the diagonal elements of the final H_k matrix
    energies = np.diag(H_k)
    
    return energies, V_k

def custom_jacobi_iteration(H, max_iter=5000, tol=1e-10):
    """
    Diagonalizes a symmetric matrix using the Jacobi Eigenvalue Algorithm.
    
    Parameters
    ----------
    H : np.ndarray
        The symmetric Hamiltonian matrix.
    max_iter : int
        The maximum number of rotations.
    tol : float
        The tolerance for convergence (off-diagonal elements near zero).
        
    Returns
    -------
    energies : np.ndarray
        1D array of eigenvalues.
    eigenvectors : np.ndarray
        2D array of eigenvectors (columns).
    """
    H_k = np.copy(H)
    n = H_k.shape[0]
    V_k = np.eye(n)
    
    for i in range(max_iter):
        # 1. Find the largest off-diagonal element
        # We temporarily set the diagonal to zero so np.argmax only searches off-diagonals
        H_off_diagonal = H_k - np.diag(np.diag(H_k))
        
        # np.unravel_index converts the flattened 1D max index back into (row, col) coordinates
        p, q = np.unravel_index(np.argmax(np.abs(H_off_diagonal)), H_k.shape)
        
        # 2. Convergence Check
        if np.abs(H_k[p, q]) < tol:
            print(f"Jacobi Algorithm converged in {i} iterations.")
            break
            
        # 3. Calculate the rotation angle theta
        if H_k[p, p] == H_k[q, q]:
            # If diagonal elements are equal, tan(2*theta) approaches infinity
            theta = np.pi / 4.0
        else:
            theta = 0.5 * np.arctan(2.0 * H_k[p, q] / (H_k[p, p] - H_k[q, q]))
            
        # 4. Construct the Orthogonal Rotation Matrix (P)
        P = np.eye(n)
        P[p, p] = np.cos(theta)
        P[q, q] = np.cos(theta)
        P[p, q] = np.sin(theta)
        P[q, p] = -np.sin(theta)
        
        # 5. Apply the Similarity Transformation (H = P^T * H * P)
        H_k = P.T @ H_k @ P
        
        # 6. Accumulate the Eigenvectors (V = V * P)
        V_k = V_k @ P
        
    else:
        print(f"Warning: Jacobi Algorithm did not fully converge within {max_iter} iterations.")

    # The eigenvalues are the isolated elements on the main diagonal
    energies = np.diag(H_k)
    
    return energies, V_k

def solve_hamiltonian(H, N_STATES, eigen_method):
    """Diagonalizes the Hamiltonian and returns the lowest N_STATES.
    
    Parameters
    ----------
    H : np.ndarray
        The Hamiltonian matrix to diagonalize.
    N_STATES : int
        The number of lowest energy states to return.
    eigen_method : str
        The method to use for eigenvalue decomposition ('scipy', 'QR', 'Jacobi').

    Returns
    -------
    tuple
        A tuple containing:
        - sorted_energies: np.ndarray of the lowest N_STATES eigenvalues (energies).
        - sorted_wavefunctions: np.ndarray of the corresponding eigenvectors (wavefunctions).
    """
    # Solve secular equation, which is just the eigenvalue problem for the Hamiltonian matrix H.
    if eigen_method in EIGEN_PROBLEM_DISPATCH:
        eigen_method_func = EIGEN_PROBLEM_DISPATCH[eigen_method]

        energies, wavefunctions = eigen_method_func(H)
    
    else:
        raise ValueError(f"Invalid eigenvalue method: {eigen_method}")
        
    # Sort eigenvalues from lowest to highest
    sort_indices = np.argsort(energies)
    sorted_energies = energies[sort_indices]
    sorted_wavefunctions = wavefunctions[:, sort_indices]
    
    return sorted_energies[:N_STATES], sorted_wavefunctions[:, :N_STATES]

# The "Strategy" Dispatcher (OOP paradigm)
INTEGRATION_DISPATCH = {
    'trapezoidal': integral_solve_using_trapezoidal_rule,
    'romberg': integral_solve_using_romberg_method,
    'simpsons': integral_solve_using_simpsons_rule,
    'scipy': scipy_quad_wrapper
}

EIGEN_PROBLEM_DISPATCH = {
    'Jacobi': custom_jacobi_iteration,
    'QR': custom_qr_algorithm,
    'scipy': eigh
}

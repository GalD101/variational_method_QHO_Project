from scipy.integrate import quad

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
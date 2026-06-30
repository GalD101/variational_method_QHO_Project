# --- Physical Parameters ---
# We parameterize the physical constants so the equations remain exact.
# They default to 1.0 (Natural Units) to prevent floating-point underflow.
HBAR = 1.0            # Reduced Planck's constant
MASS = 1.0            # Particle mass
OMEGA = 1.0           # Oscillator frequency

# --- Baseline Defaults (For single-run testing) ---
# We label them as defaults so we know they are meant to be overwritten in loops
DEFAULT_L = 30.0#10.0      
DEFAULT_N_BASIS = 50  
DEFAULT_N_STATES = 5

# --- Execution Flags ---
# These control the logic flow, making debugging much easier.
DEBUG_PRINT = True       # Set to True to print matrix values (keep N_BASIS small if True!)
PLOT_ERROR = True        # Flag to trigger the matplotlib convergence graph
PLOT_WAVEFUNCTION = True # Flag to trigger visualizing the final wavefunctions

# --- Numerical Methods ---
integration_method = 'trapezoidal' # Options: 'analytical', 'scipy', 'trapezoidal', 'romberg', 'simpsons'
eigenvalue_method = 'scipy'        # Options: 'scipy', 'QR', 'Jacobi'
N_GRID = 1000                      # Number of grid intervals for custom integration methods

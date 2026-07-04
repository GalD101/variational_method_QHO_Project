# --- Physical Parameters ---
# We parameterize the physical constants so the equations remain exact.
# They default to 1.0 (Natural Units) to prevent floating-point underflow.
HBAR = 1.0  # Reduced Planck's constant
MASS = 1.0  # Particle mass
OMEGA = 1.0 # Oscillator frequency

# --- Baseline Defaults (For single-run testing) ---
# We label them as defaults so we know they are meant to be overwritten in loops
DEFAULT_L = 25.0      
DEFAULT_N_BASIS = 100  
DEFAULT_N_STATES = 38

# --- Numerical Methods ---
integration_method = 'analytical' # Options: 'analytical', 'scipy', 'trapezoidal', 'romberg', 'simpsons'
eigenvalue_method = 'scipy'       # Options: 'scipy', 'QR', 'Jacobi'
N_GRID = 1500                     # Number of grid intervals for custom integration methods

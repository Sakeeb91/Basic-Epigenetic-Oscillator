#!/usr/bin/env python3
# Epigenetic Oscillator Simulator based on Goodwin Equations (Eq. 14)

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
from datetime import datetime
import shutil

# Create plots directory if it doesn't exist
if not os.path.exists('plots'):
    os.makedirs('plots')

# Create a unique timestamped folder for this run
run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
run_folder = f"plots/run_{run_timestamp}"
os.makedirs(run_folder, exist_ok=True)

# Create sample_outputs directory if it doesn't exist
if not os.path.exists('sample_outputs'):
    os.makedirs('sample_outputs')

def calculate_steady_state(params):
    """
    Calculate the steady state values (p_i, q_i) for (X, Y) based on the parameters.
    
    For a biologically meaningful steady state (p_i > 0, q_i > 0), we require a_i/b_i > A_i.
    
    p_i = (β_i/α_i)q_i
    q_i = (a_i/b_i - A_i) / k_i
    """
    alpha_i, beta_i, a_i, b_i, k_i, A_i = params
    
    # Check if a_i/b_i > A_i for biologically meaningful steady state
    if a_i/b_i <= A_i:
        raise ValueError(f"Parameters do not satisfy the condition a_i/b_i > A_i for meaningful steady state. a_i/b_i = {a_i/b_i}, A_i = {A_i}")
    
    # Calculate steady state values
    q_i = (a_i/b_i - A_i) / k_i
    p_i = (beta_i/alpha_i) * q_i
    
    return p_i, q_i

def goodwin_oscillator_ode(t, state, params):
    """
    Goodwin oscillator ODE system (Equation 14)
    dX/dt = a_i / (A_i + k_i*Y) - b_i
    dY/dt = α_i*X - β_i*Y
    
    Parameters:
    - t: time point
    - state: [X, Y] array of state variables
    - params: [alpha_i, beta_i, a_i, b_i, k_i, A_i] parameter list
    """
    X, Y = state
    alpha_i, beta_i, a_i, b_i, k_i, A_i = params
    
    # Ensure denominator is positive (handle potential numerical issues)
    denom = A_i + k_i * max(Y, 0)
    
    dX_dt = (a_i / denom) - b_i
    dY_dt = alpha_i * X - beta_i * Y
    
    return [dX_dt, dY_dt]

def calculate_talandic_energy(X, Y, params):
    """
    Calculate the Talandic Energy G (Equation 15) for the Goodwin system
    G(X, Y) = (α_i/2)*X^2 - β_i*X + b_i*Y + (a_i/k_i)*log(A_i + k_i*Y)
    
    This should be conserved for a given trajectory.
    """
    alpha_i, beta_i, a_i, b_i, k_i, A_i = params
    
    # Ensure argument of log is positive (handle potential numerical issues)
    log_term_arg = A_i + k_i * Y
    safe_log_term_arg = np.maximum(log_term_arg, 1e-9)  # Prevent log(<=0)
    
    G = (alpha_i / 2.0) * X**2 - beta_i * X + b_i * Y + (a_i / k_i) * np.log(safe_log_term_arg)
    
    return G

def simulate_oscillator(params=None, t_span=(0, 500), t_points=1000, initial_state=None, steady_state_factor=None):
    """
    Simulate the Goodwin oscillator model
    
    Parameters:
    - params: [alpha_i, beta_i, a_i, b_i, k_i, A_i] parameter list
    - t_span: (t_start, t_end) tuple for simulation timespan
    - t_points: number of time points to evaluate
    - initial_state: [X0, Y0] initial conditions, if None calculated from steady state and factors
    - steady_state_factor: [X_factor, Y_factor] to set initial state relative to steady state
    
    Returns:
    - t: time points
    - states: array of state variables at each time point
    - G: Talandic Energy at each time point
    - steady_state: (p_i, q_i) steady state values
    """
    # Default parameters if none provided
    if params is None:
        # Example parameters as suggested in the plan
        alpha_i = 0.2   # protein/mRNA/min
        beta_i = 0.02   # 1/min (~50 min protein half-life)
        a_i = 5.0       # mRNA/min (max synthesis rate)
        b_i = 0.1       # mRNA/min (degradation rate)
        k_i = 0.05      # 1/protein_conc
        A_i = 1.0       # dimensionless
        
        params = [alpha_i, beta_i, a_i, b_i, k_i, A_i]
    
    # Calculate steady state
    p_i, q_i = calculate_steady_state(params)
    steady_state = (p_i, q_i)
    
    # Set initial conditions
    if initial_state is None:
        if steady_state_factor is None:
            # Default: X0 = 1.5*p_i, Y0 = 0.8*q_i
            steady_state_factor = [1.5, 0.8]
        
        X0 = p_i * steady_state_factor[0]
        Y0 = q_i * steady_state_factor[1]
        initial_state = [X0, Y0]
    
    # Create time points for evaluation
    t_eval = np.linspace(t_span[0], t_span[1], t_points)
    
    # Solve the ODE system
    solution = solve_ivp(
        goodwin_oscillator_ode, t_span, initial_state, args=(params,),
        method='RK45', t_eval=t_eval, rtol=1e-6, atol=1e-9
    )
    
    # Calculate Talandic Energy for each point
    X_sol = solution.y[0]
    Y_sol = solution.y[1]
    G = np.array([calculate_talandic_energy(x, y, params) for x, y in zip(X_sol, Y_sol)])
    
    return solution.t, solution.y, G, steady_state

def plot_results(t, states, G, steady_state, params):
    """
    Plot the simulation results and save to files
    
    Parameters:
    - t: time points
    - states: array of state variables at each time point
    - G: Talandic Energy at each time point
    - steady_state: (p_i, q_i) steady state values
    - params: [alpha_i, beta_i, a_i, b_i, k_i, A_i] parameter list
    """
    X_sol, Y_sol = states
    p_i, q_i = steady_state
    alpha_i, beta_i, a_i, b_i, k_i, A_i = params
    
    # 1. Time Series plot
    plt.figure(figsize=(10, 6))
    plt.plot(t, X_sol, label='mRNA (X)')
    plt.plot(t, Y_sol, label='Protein (Y)')
    
    # Add horizontal lines for steady state values
    plt.axhline(y=p_i, color='r', linestyle='--', alpha=0.5, label=f'X steady state ({p_i:.2f})')
    plt.axhline(y=q_i, color='g', linestyle='--', alpha=0.5, label=f'Y steady state ({q_i:.2f})')
    
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.title('Goodwin Oscillator Time Series (Eq. 14)')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    filename = 'time_series.png'
    filepath = os.path.join(run_folder, filename)
    plt.savefig(filepath)
    
    # Create markdown documentation
    md_filename = 'time_series.md'
    md_filepath = os.path.join(run_folder, md_filename)
    with open(md_filepath, 'w') as f:
        f.write(f'# Time Series Plot for Goodwin Oscillator (Eq. 14)\n\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        f.write('## Parameters:\n')
        f.write(f'- α_i: {alpha_i}\n')
        f.write(f'- β_i: {beta_i}\n')
        f.write(f'- a_i: {a_i}\n')
        f.write(f'- b_i: {b_i}\n')
        f.write(f'- k_i: {k_i}\n')
        f.write(f'- A_i: {A_i}\n\n')
        f.write('## Steady State Values:\n')
        f.write(f'- p_i (X steady state): {p_i:.4f}\n')
        f.write(f'- q_i (Y steady state): {q_i:.4f}\n\n')
        f.write('## Description:\n')
        f.write('This plot shows the temporal evolution of mRNA (X) and protein (Y) concentrations in the Goodwin oscillator model.\n')
        f.write('The oscillations represent the cyclic behavior resulting from the negative feedback loop.\n')
        f.write('Note the phase lag between mRNA (X) and protein (Y) peaks, which is characteristic of this type of oscillator.\n')
        f.write(f'\n![Time Series Plot](./{filename})\n')
    
    # Close the plot to free memory
    plt.close()
    
    # 2. Phase Portrait
    plt.figure(figsize=(8, 8))
    plt.plot(X_sol, Y_sol)
    
    # Mark the initial condition
    plt.scatter(X_sol[0], Y_sol[0], color='blue', s=100, label='Initial state')
    
    # Mark the steady state point
    plt.scatter(p_i, q_i, color='red', marker='x', s=100, label='Steady state (p_i, q_i)')
    
    plt.xlabel('mRNA (X)')
    plt.ylabel('Protein (Y)')
    plt.title('Phase Portrait of Goodwin Oscillator')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    filename = 'phase_portrait.png'
    filepath = os.path.join(run_folder, filename)
    plt.savefig(filepath)
    
    # Create markdown documentation
    md_filename = 'phase_portrait.md'
    md_filepath = os.path.join(run_folder, md_filename)
    with open(md_filepath, 'w') as f:
        f.write(f'# Phase Portrait for Goodwin Oscillator (Eq. 14)\n\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        f.write('## Parameters:\n')
        f.write(f'- α_i: {alpha_i}\n')
        f.write(f'- β_i: {beta_i}\n')
        f.write(f'- a_i: {a_i}\n')
        f.write(f'- b_i: {b_i}\n')
        f.write(f'- k_i: {k_i}\n')
        f.write(f'- A_i: {A_i}\n\n')
        f.write('## Steady State Values:\n')
        f.write(f'- p_i (X steady state): {p_i:.4f}\n')
        f.write(f'- q_i (Y steady state): {q_i:.4f}\n\n')
        f.write('## Description:\n')
        f.write('This phase portrait shows the trajectory of the system in the mRNA-Protein state space.\n')
        f.write('The closed loop indicates a limit cycle, which is characteristic of sustained oscillations.\n')
        f.write('The trajectory encircles the steady state point (p_i, q_i), which is an unstable equilibrium in this model.\n')
        f.write(f'\n![Phase Portrait](./{filename})\n')
    
    # Close the plot to free memory
    plt.close()
    
    # 3. Talandic Energy plot
    plt.figure(figsize=(10, 6))
    plt.plot(t, G)
    
    plt.xlabel('Time')
    plt.ylabel('Talandic Energy G')
    plt.title('Conservation of Talandic Energy (Eq. 15)')
    plt.grid(True)
    
    # Calculate statistics to assess conservation
    G_mean = np.mean(G)
    G_std = np.std(G)
    G_rel_std = G_std / abs(G_mean) * 100  # Relative standard deviation in percentage
    
    plt.figtext(0.5, 0.01, f'Mean G: {G_mean:.6f}, Std: {G_std:.6f}, Rel. Std: {G_rel_std:.6f}%', 
                ha='center', fontsize=10, bbox={'facecolor':'white', 'alpha':0.5, 'pad':5})
    
    # Save the plot
    filename = 'talandic_energy.png'
    filepath = os.path.join(run_folder, filename)
    plt.savefig(filepath)
    
    # Create markdown documentation
    md_filename = 'talandic_energy.md'
    md_filepath = os.path.join(run_folder, md_filename)
    with open(md_filepath, 'w') as f:
        f.write(f'# Talandic Energy Conservation for Goodwin Oscillator (Eq. 14 & 15)\n\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        f.write('## Parameters:\n')
        f.write(f'- α_i: {alpha_i}\n')
        f.write(f'- β_i: {beta_i}\n')
        f.write(f'- a_i: {a_i}\n')
        f.write(f'- b_i: {b_i}\n')
        f.write(f'- k_i: {k_i}\n')
        f.write(f'- A_i: {A_i}\n\n')
        f.write('## Talandic Energy Formula (Eq. 15):\n')
        f.write('G(X, Y) = (α_i/2)*X^2 - β_i*X + b_i*Y + (a_i/k_i)*log(A_i + k_i*Y)\n\n')
        f.write('## Conservation Statistics:\n')
        f.write(f'- Mean G: {G_mean:.6f}\n')
        f.write(f'- Standard Deviation: {G_std:.6f}\n')
        f.write(f'- Relative Std: {G_rel_std:.6f}%\n\n')
        f.write('## Description:\n')
        f.write('This plot shows the Talandic Energy (G) over time, which should be conserved for the idealized model.\n')
        f.write('The degree to which G remains constant is a measure of the accuracy of the numerical integration.\n')
        f.write('A small relative standard deviation confirms the conservative nature of the system.\n')
        f.write(f'\n![Talandic Energy](./{filename})\n')
    
    # Close the plot to free memory
    plt.close()

def parameter_study(parameter='k_i', range_values=None, base_params=None, t_span=(0, 500), 
                    t_points=1000, steady_state_factor=None):
    """
    Perform a parameter study by varying one parameter and observing its effect
    
    Parameters:
    - parameter: name of parameter to vary ('alpha_i', 'beta_i', 'a_i', 'b_i', 'k_i', 'A_i')
    - range_values: list of values for the parameter
    - base_params: base parameter values [alpha_i, beta_i, a_i, b_i, k_i, A_i]
    - t_span, t_points: simulation parameters
    - steady_state_factor: [X_factor, Y_factor] to set initial state relative to steady state
    """
    # Default base parameters if none provided
    if base_params is None:
        # As suggested in the plan
        alpha_i = 0.2   # protein/mRNA/min
        beta_i = 0.02   # 1/min (~50 min protein half-life)
        a_i = 5.0       # mRNA/min (max synthesis rate)
        b_i = 0.1       # mRNA/min (degradation rate)
        k_i = 0.05      # 1/protein_conc
        A_i = 1.0       # dimensionless
        
        base_params = [alpha_i, beta_i, a_i, b_i, k_i, A_i]
    
    # Default range values if none provided
    if range_values is None:
        if parameter == 'k_i':
            range_values = [0.02, 0.05, 0.1]  # Different repression strengths
        elif parameter == 'alpha_i':
            range_values = [0.1, 0.2, 0.4]    # Different protein synthesis rates
        elif parameter == 'beta_i':
            range_values = [0.01, 0.02, 0.04]  # Different protein degradation rates
        elif parameter == 'a_i':
            range_values = [2.5, 5.0, 10.0]    # Different mRNA synthesis rates
        elif parameter == 'b_i':
            range_values = [0.05, 0.1, 0.2]    # Different mRNA degradation rates
        elif parameter == 'A_i':
            range_values = [0.5, 1.0, 2.0]     # Different repression threshold values
    
    # Parameter index in the parameter list
    param_dict = {'alpha_i': 0, 'beta_i': 1, 'a_i': 2, 'b_i': 3, 'k_i': 4, 'A_i': 5}
    param_index = param_dict[parameter]
    
    # 1. Time series plot for varying parameter
    plt.figure(figsize=(12, 8))
    steady_states = []
    
    for value in range_values:
        # Create a copy of base parameters and update the parameter to vary
        params = base_params.copy()
        params[param_index] = value
        
        # Get parameter values for labeling
        alpha_i, beta_i, a_i, b_i, k_i, A_i = params
        
        # Run simulation (use the same initial conditions adjustment from steady state)
        t, states, _, steady_state = simulate_oscillator(
            params=params, t_span=t_span, t_points=t_points,
            steady_state_factor=steady_state_factor
        )
        
        # Store steady state for each parameter value
        steady_states.append(steady_state)
        
        # Plot the mRNA (X) variable
        plt.plot(t, states[0], label=f'{parameter}={value}')
    
    plt.xlabel('Time')
    plt.ylabel('mRNA (X) Concentration')
    plt.title(f'Effect of {parameter} on Oscillation Dynamics')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    filename = f'parameter_study_{parameter}_time_series.png'
    filepath = os.path.join(run_folder, filename)
    plt.savefig(filepath)
    
    # Create markdown documentation
    md_filename = f'parameter_study_{parameter}_time_series.md'
    md_filepath = os.path.join(run_folder, md_filename)
    with open(md_filepath, 'w') as f:
        f.write(f'# Parameter Study: Effect of {parameter} on Goodwin Oscillator\n\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        f.write('## Base Parameters:\n')
        params_names = ['α_i', 'β_i', 'a_i', 'b_i', 'k_i', 'A_i']
        for i, (name, value) in enumerate(zip(params_names, base_params)):
            if i != param_index:  # Don't list the parameter being varied
                f.write(f'- {name}: {value}\n')
        
        f.write(f'\n## Varied Parameter: {parameter}\n')
        f.write(f'Values: {range_values}\n\n')
        
        f.write('## Steady State Values:\n')
        for i, (value, ss) in enumerate(zip(range_values, steady_states)):
            p_i, q_i = ss
            f.write(f'For {parameter}={value}:\n')
            f.write(f'- p_i (X steady state): {p_i:.4f}\n')
            f.write(f'- q_i (Y steady state): {q_i:.4f}\n\n')
        
        f.write('## Description:\n')
        f.write(f'This plot shows how varying the parameter {parameter} affects the dynamics of the Goodwin oscillator.\n')
        f.write('Changes in oscillation amplitude, period, and waveform can be observed as the parameter value changes.\n')
        f.write(f'\n![Parameter Study Time Series](./{filename})\n')
    
    # Close the plot to free memory
    plt.close()
    
    # 2. Phase portraits for varying parameter
    plt.figure(figsize=(12, 10))
    
    for i, value in enumerate(range_values):
        # Create a copy of base parameters and update the parameter to vary
        params = base_params.copy()
        params[param_index] = value
        
        # Run simulation
        t, states, _, steady_state = simulate_oscillator(
            params=params, t_span=t_span, t_points=t_points,
            steady_state_factor=steady_state_factor
        )
        
        p_i, q_i = steady_state
        
        # Create subplot
        ax = plt.subplot(2, 2, i+1)
        ax.plot(states[0], states[1])
        
        # Mark the initial condition
        ax.scatter(states[0][0], states[1][0], color='blue', s=50, label='Initial state')
        
        # Mark the steady state point
        ax.scatter(p_i, q_i, color='red', marker='x', s=50, label='Steady state')
        
        ax.set_xlabel('mRNA (X)')
        ax.set_ylabel('Protein (Y)')
        ax.set_title(f'{parameter}={value}')
        ax.legend()
        ax.grid(True)
    
    plt.tight_layout()
    
    # Save the plot
    filename = f'parameter_study_{parameter}_phase_portraits.png'
    filepath = os.path.join(run_folder, filename)
    plt.savefig(filepath)
    
    # Create markdown documentation
    md_filename = f'parameter_study_{parameter}_phase_portraits.md'
    md_filepath = os.path.join(run_folder, md_filename)
    with open(md_filepath, 'w') as f:
        f.write(f'# Phase Portraits: Effect of {parameter} on Goodwin Oscillator\n\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        f.write('## Base Parameters:\n')
        params_names = ['α_i', 'β_i', 'a_i', 'b_i', 'k_i', 'A_i']
        for i, (name, value) in enumerate(zip(params_names, base_params)):
            if i != param_index:  # Don't list the parameter being varied
                f.write(f'- {name}: {value}\n')
        
        f.write(f'\n## Varied Parameter: {parameter}\n')
        f.write(f'Values: {range_values}\n\n')
        
        f.write('## Description:\n')
        f.write(f'These phase portraits show how varying the parameter {parameter} affects the shape and size of the limit cycle.\n')
        f.write('Changes in the dynamics of the system are reflected in the geometry of the trajectories in phase space.\n')
        f.write('Note how the steady state point (marked with x) changes with the parameter value.\n')
        f.write(f'\n![Phase Portraits](./{filename})\n')
    
    # Close the plot to free memory
    plt.close()

def create_readme():
    """Create a README file for the run folder summarizing the content"""
    readme_path = os.path.join(run_folder, "README.md")
    with open(readme_path, 'w') as f:
        f.write(f'# Goodwin Oscillator Simulation Results (Eq. 14 & 15)\n\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        f.write('## Contents\n\n')
        f.write('### Basic Simulation\n')
        f.write('- [Time Series](./time_series.md)\n')
        f.write('- [Phase Portrait](./phase_portrait.md)\n')
        f.write('- [Talandic Energy Conservation](./talandic_energy.md)\n\n')
        f.write('### Parameter Studies\n')
        f.write('- Repression Strength (k_i)\n')
        f.write('  - [Time Series Effect](./parameter_study_k_i_time_series.md)\n')
        f.write('  - [Phase Portraits](./parameter_study_k_i_phase_portraits.md)\n\n')

def save_sample_outputs(t, states, G, steady_state, params):
    """
    Save sample output plots to the sample_outputs directory for README display
    
    Parameters:
    - t: time points
    - states: array of state variables at each time point
    - G: Talandic Energy at each time point
    - steady_state: (p_i, q_i) steady state values
    - params: [alpha_i, beta_i, a_i, b_i, k_i, A_i] parameter list
    """
    X_sol, Y_sol = states
    p_i, q_i = steady_state
    
    # 1. Time Series plot
    plt.figure(figsize=(10, 6))
    plt.plot(t, X_sol, label='mRNA (X)')
    plt.plot(t, Y_sol, label='Protein (Y)')
    
    # Add horizontal lines for steady state values
    plt.axhline(y=p_i, color='r', linestyle='--', alpha=0.5, label=f'X steady state ({p_i:.2f})')
    plt.axhline(y=q_i, color='g', linestyle='--', alpha=0.5, label=f'Y steady state ({q_i:.2f})')
    
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.title('Goodwin Oscillator Time Series (Eq. 14)')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    plt.savefig('sample_outputs/time_series.png')
    plt.close()
    
    # 2. Phase Portrait
    plt.figure(figsize=(8, 8))
    plt.plot(X_sol, Y_sol)
    
    # Mark the initial condition
    plt.scatter(X_sol[0], Y_sol[0], color='blue', s=100, label='Initial state')
    
    # Mark the steady state point
    plt.scatter(p_i, q_i, color='red', marker='x', s=100, label='Steady state (p_i, q_i)')
    
    plt.xlabel('mRNA (X)')
    plt.ylabel('Protein (Y)')
    plt.title('Phase Portrait of Goodwin Oscillator')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    plt.savefig('sample_outputs/phase_portrait.png')
    plt.close()
    
    # 3. Talandic Energy plot
    plt.figure(figsize=(10, 6))
    plt.plot(t, G)
    
    plt.xlabel('Time')
    plt.ylabel('Talandic Energy G')
    plt.title('Conservation of Talandic Energy (Eq. 15)')
    plt.grid(True)
    
    # Calculate statistics to assess conservation
    G_mean = np.mean(G)
    G_std = np.std(G)
    G_rel_std = G_std / abs(G_mean) * 100  # Relative standard deviation in percentage
    
    plt.figtext(0.5, 0.01, f'Mean G: {G_mean:.6f}, Std: {G_std:.6f}, Rel. Std: {G_rel_std:.6f}%', 
                ha='center', fontsize=10, bbox={'facecolor':'white', 'alpha':0.5, 'pad':5})
    
    # Save the plot
    plt.savefig('sample_outputs/talandic_energy.png')
    plt.close()

if __name__ == "__main__":
    # Define parameters as suggested in the plan
    alpha_i = 0.2   # protein/mRNA/min
    beta_i = 0.02   # 1/min (~50 min protein half-life)
    a_i = 5.0       # mRNA/min (max synthesis rate)
    b_i = 0.1       # mRNA/min (degradation rate)
    k_i = 0.05      # 1/protein_conc
    A_i = 1.0       # dimensionless
    
    params = [alpha_i, beta_i, a_i, b_i, k_i, A_i]
    
    # 1. Basic simulation with default parameters
    # Start with X at 50% above steady state and Y at 20% below
    steady_state_factor = [1.5, 0.8]
    t, states, G, steady_state = simulate_oscillator(
        params=params, t_span=(0, 500), t_points=1000,
        steady_state_factor=steady_state_factor
    )
    
    # Plot results
    plot_results(t, states, G, steady_state, params)
    
    # Save sample outputs for README display
    save_sample_outputs(t, states, G, steady_state, params)
    
    # 2. Parameter study for repression strength k_i
    parameter_study(
        parameter='k_i', 
        range_values=[0.02, 0.05, 0.1],
        base_params=params,
        steady_state_factor=steady_state_factor
    )
    
    # Create README for the run folder
    create_readme()
    
    print(f"Simulation completed. All plots have been saved to the '{run_folder}' directory.") 
import tomli
from typing import Tuple
import _CavSim


def read_config(filename: str) -> Tuple[_CavSim.CavitySetup, float]:
    """
    Reads a TOML configuration file and returns a CavitySetup object with solver tolerance.
    
    TOML structure:
    [domain_dimensions]
    length = <float>
    height = <float>
    
    [physical_properties]
    density = <float>
    viscosity = <float>
    lid_velocity = <float>
    
    [numerical_properties]
    nx = <int>
    ny = <int>
    solver_tolerance = <float>
    
    Args:
        filename: Path to the TOML configuration file
        
    Returns:
        Tuple of (CavitySetup object, solver_tolerance)
    """
    try:
        with open(filename, 'rb') as f:
            config = tomli.load(f)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        raise
    except tomli.TOMLDecodeError as e:
        print(f"Error: Failed to parse TOML file '{filename}': {e}")
        raise
    
    # Extract domain dimensions
    domain = config.get('domain_dimensions', {})
    length = domain.get('length')
    height = domain.get('height')
    
    if length is None or height is None:
        raise KeyError("Missing 'length' or 'height' in domain_dimensions section")
    
    # Extract physical properties
    physics = config.get('physical_properties', {})
    density = physics.get('density')
    viscosity = physics.get('viscosity')
    lid_velocity = physics.get('lid_velocity')
    
    if density is None or viscosity is None or lid_velocity is None:
        raise KeyError("Missing 'density', 'viscosity', or 'lid_velocity' in physical_properties section")
    
    # Extract numerical properties
    numerical = config.get('numerical_properties', {})
    nx = numerical.get('nx')
    ny = numerical.get('ny')
    solver_tolerance = numerical.get('solver_tolerance')
    
    if nx is None or ny is None or solver_tolerance is None:
        raise KeyError("Missing 'nx', 'ny', or 'solver_tolerance' in numerical_properties section")
    
    # Create CavitySetup object
    cavity_setup = _CavSim.CavitySetup(
        length=length,
        height=height,
        n_x=nx,
        n_y=ny,
        rho=density,
        mu=viscosity,
        U_lid=lid_velocity
    )
    
    return cavity_setup, solver_tolerance


if __name__ == "__main__":
    # Example usage
    cavity_setup, solver_tolerance = read_config("data/cavsim_config.toml")
    print("Configuration loaded:")
    print(f"CavitySetup: length={cavity_setup.L}, height={cavity_setup.H}, nx={cavity_setup.n_x}, ny={cavity_setup.n_y}")
    print(f"Physical: rho={cavity_setup.rho}, mu={cavity_setup.mu}, U={cavity_setup.U}")
    print(f"Grid: dx={cavity_setup.dx}, dy={cavity_setup.dy}")
    print(f"Solver tolerance: {solver_tolerance}")

import pytest
import sys
import os

from IO import read_config


class TestReadConfig:
    """Test suite for the read_config function."""
    
    @pytest.fixture
    def config_file(self):
        """Path to the test configuration file."""
        return os.path.join(os.path.dirname(__file__), '../../../data/cavsim_config.toml')
    
    def test_read_config_returns_tuple(self, config_file):
        """Test that read_config returns a tuple."""
        result = read_config(config_file)
        assert isinstance(result, tuple), "read_config should return a tuple"
        assert len(result) == 2, "read_config should return a tuple of 2 elements"
    
    def test_read_config_returns_cavity_setup_and_float(self, config_file):
        """Test that read_config returns CavitySetup object and float."""
        cavity_setup, solver_tolerance = read_config(config_file)
        
        # Check CavitySetup object exists and has expected attributes
        assert hasattr(cavity_setup, 'L'), "CavitySetup should have 'L' attribute"
        assert hasattr(cavity_setup, 'H'), "CavitySetup should have 'H' attribute"
        assert hasattr(cavity_setup, 'n_x'), "CavitySetup should have 'n_x' attribute"
        assert hasattr(cavity_setup, 'n_y'), "CavitySetup should have 'n_y' attribute"
        assert hasattr(cavity_setup, 'dx'), "CavitySetup should have 'dx' attribute"
        assert hasattr(cavity_setup, 'dy'), "CavitySetup should have 'dy' attribute"
        assert hasattr(cavity_setup, 'rho'), "CavitySetup should have 'rho' attribute"
        assert hasattr(cavity_setup, 'mu'), "CavitySetup should have 'mu' attribute"
        assert hasattr(cavity_setup, 'U'), "CavitySetup should have 'U' attribute"
        
        # Check solver_tolerance is a float
        assert isinstance(solver_tolerance, float), "solver_tolerance should be a float"
    
    def test_read_config_values(self, config_file):
        """Test that read_config extracts correct values from TOML."""
        cavity_setup, solver_tolerance = read_config(config_file)
        
        # Check domain dimensions
        assert cavity_setup.L == 1.0, "Domain length should be 1.0"
        assert cavity_setup.H == 1.0, "Domain height should be 1.0"
        
        # Check grid divisions
        assert cavity_setup.n_x == 10, "Grid should have 10 divisions in x"
        assert cavity_setup.n_y == 10, "Grid should have 10 divisions in y"
        
        # Check physical properties
        assert cavity_setup.rho == 1.0, "Density should be 1.0"
        assert cavity_setup.mu == 0.01, "Viscosity should be 0.01"
        assert cavity_setup.U == 200.0, "Lid velocity should be 200.0"
        
        # Check solver tolerance
        assert solver_tolerance == 1e-4, "Solver tolerance should be 1e-4"
    
    def test_read_config_grid_spacing(self, config_file):
        """Test that read_config calculates correct grid spacing."""
        cavity_setup, _ = read_config(config_file)
        
        # dx should be L / n_x
        expected_dx = 1.0 / 10
        assert abs(cavity_setup.dx - expected_dx) < 1e-10, f"dx should be {expected_dx}"
        
        # dy should be H / n_y
        expected_dy = 1.0 / 10
        assert abs(cavity_setup.dy - expected_dy) < 1e-10, f"dy should be {expected_dy}"
    
    def test_read_config_file_not_found(self):
        """Test that read_config raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            read_config("nonexistent_file.toml")
    
    def test_read_config_missing_section(self, tmp_path):
        """Test that read_config raises KeyError for missing sections."""
        # Create temporary TOML file with missing sections
        incomplete_toml = tmp_path / "incomplete.toml"
        incomplete_toml.write_text("[domain_dimensions]\nlength = 1.0\n")
        
        with pytest.raises(KeyError):
            read_config(str(incomplete_toml))


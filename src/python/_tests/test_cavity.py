import pytest


class TestCavitySetup:
	"""Test suite for the CavitySetup class."""
	
	@pytest.fixture
	def create_cavity_setup(self):
		"""Create a CavitySetup instance with test parameters."""
		from _CavSim import CavitySetup
		return CavitySetup(
			length=200.0,  # m
			height=400.0,  # m
			n_x=10,
			n_y=25,
			rho=900.0,  # kg/m3
			mu=1.002e-3,  # Pa.s
			U_lid=10.0,  # m/s
		)
	
	def test_cavity_setup_attributes(self, create_cavity_setup):
		"""Test that CavitySetup object has all expected attributes."""
		cav_setup = create_cavity_setup
		
		# Check all attributes exist
		assert hasattr(cav_setup, 'L'), "CavitySetup should have 'L' attribute"
		assert hasattr(cav_setup, 'H'), "CavitySetup should have 'H' attribute"
		assert hasattr(cav_setup, 'n_x'), "CavitySetup should have 'n_x' attribute"
		assert hasattr(cav_setup, 'n_y'), "CavitySetup should have 'n_y' attribute"
		assert hasattr(cav_setup, 'dx'), "CavitySetup should have 'dx' attribute"
		assert hasattr(cav_setup, 'dy'), "CavitySetup should have 'dy' attribute"
		assert hasattr(cav_setup, 'rho'), "CavitySetup should have 'rho' attribute"
		assert hasattr(cav_setup, 'mu'), "CavitySetup should have 'mu' attribute"
		assert hasattr(cav_setup, 'U'), "CavitySetup should have 'U' attribute"
	
	def test_cavity_setup_grid_spacing(self, create_cavity_setup):
		"""Test that CavitySetup calculates correct grid spacing."""
		cav_setup = create_cavity_setup
		
		# dx should be L / n_x
		expected_dx = 200 / 10
		assert cav_setup.dx == expected_dx, f"dx should be {expected_dx}"
		
		# dy should be H / n_y
		expected_dy = 400 / 25
		assert cav_setup.dy == expected_dy, f"dy should be {expected_dy}"
	
	def test_cavity_setup_domain_dimensions(self, create_cavity_setup):
		"""Test that CavitySetup stores correct domain dimensions."""
		cav_setup = create_cavity_setup
		
		assert abs(cav_setup.L - 200.0) < 1e-10, "Domain length should be 200.0"
		assert abs(cav_setup.H - 400.0) < 1e-10, "Domain height should be 400.0"
	
	def test_cavity_setup_grid_divisions(self, create_cavity_setup):
		"""Test that CavitySetup stores correct grid divisions."""
		cav_setup = create_cavity_setup
		
		assert cav_setup.n_x == 10, "Grid should have 10 divisions in x"
		assert cav_setup.n_y == 25, "Grid should have 25 divisions in y"
	
	def test_cavity_setup_physical_properties(self, create_cavity_setup):
		"""Test that CavitySetup stores correct physical properties."""
		cav_setup = create_cavity_setup
		
		assert abs(cav_setup.rho - 900.0) < 1e-10, "Density should be 900.0"
		assert abs(cav_setup.mu - 0.001002) < 1e-10, "Viscosity should be 0.001002"
		assert abs(cav_setup.U - 10.0) < 1e-10, "Lid velocity should be 10.0"

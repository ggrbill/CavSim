import pytest


@pytest.fixture
def create_cavity_setup():
    from _CavSim import CavitySetup
    return CavitySetup(
      length=200, # m  
      height=400, # m 
      n_x=10,
      n_y=25,
      rho=900, # kg/m3
      mu=1.002e-3, # Pa.s
      U_lid=10, # m/s
    )


def test_cavity_setup(create_cavity_setup):
    cav_setup = create_cavity_setup

    assert cav_setup.dx == 20
    assert cav_setup.dy == 16

    assert cav_setup.L == 200
    assert cav_setup.H == 400
    assert cav_setup.n_x == 10
    assert cav_setup.n_y == 25
    assert cav_setup.rho == 900
    assert cav_setup.mu == 0.001002
    assert cav_setup.U == 10

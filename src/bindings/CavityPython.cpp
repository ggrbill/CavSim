#include <pybind11/pybind11.h>

#include "Cavity.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_CavSim, m) {
    {
        using T = CavitySetup;  
        py::class_<T, std::shared_ptr<T> > c(m, "CavitySetup");
        c.def(py::init<
            double, // length,
            double, // height,
            int, // n_x, 
            int, // n_y,  
            double, // rho,   
            double, // mu,
            double // U_lid   
        >(),
        py::arg("length"),
        py::arg("height"),
        py::arg("n_x"),
        py::arg("n_y"),
        py::arg("rho"),
        py::arg("mu"),
        py::arg("U_lid"));
    }
}


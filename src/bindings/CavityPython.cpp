#include <pybind11/pybind11.h>

#include "Cavity.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_CavSim, m) {
    m.doc() = "CavSim Library";
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

        c.def_readonly("L", &T::L);
        c.def_readonly("H", &T::H);
        c.def_readonly("dx", &T::dx);
        c.def_readonly("dy", &T::dy);
        c.def_readonly("n_x", &T::n_x);
        c.def_readonly("n_y", &T::n_y);

        c.def_readonly("rho", &T::rho);
        c.def_readonly("mu", &T::mu);
        c.def_readonly("U", &T::U);
    }
}


#ifndef IO_H
#define IO_H

#include <memory>

#include "Cavity.hpp"

/*!
    Input data reader.
*/
std::shared_ptr<CavitySetup> read_input_data(std::string filename);

#endif
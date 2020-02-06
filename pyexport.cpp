#include "pyrunions.h"
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/exception_translator.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(librunIons)
{
    register_exception_translator<PyFieldException>
    (
        &PyFieldException::pyFieldExceptTrans
    );

    class_<PySimFields>("Field")
            .def("addStaticField", &PySimFields::addStaticField, args("files", "volts", "scale_mm"))
            .def("addSinField", &PySimFields::addSinField, args("files", "volts", "freq_MHz", "phase_rad", "scale_mm"))
            .def("potential", &PySimFields::potential, args("x", "y", "z", "t"))
            .def("field", &PySimFields::field, args("x", "y", "z", "t"));

}

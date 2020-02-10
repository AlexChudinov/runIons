#ifndef PYTHONFIELDS_H
#define PYTHONFIELDS_H

#include "timefield.h"
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

/**
 * @brief The PythonFieldException class handles exceptions that were rised inside PythonFieldsClass
 */
class PyFieldException : public std::exception
{
    std::string mMsg;
public:
    PyFieldException(const std::string& msg);
    virtual const char* what() const throw();

    /**
     * @brief pyFieldExceptTrans translates python exception from c++ to python
     * @param e initial exception
     */
    static void pyFieldExceptTrans(const PyFieldException& e);
};

/**
 * @brief The PythonFields class interface for fields to native python
 * datatypes
 */
class PySimFields
{
    friend class PyRunIon;
    friend class PyRunIons;
public:
    using PyList = boost::python::list;
    using PyTuple = boost::python::tuple;

    /**
     * @brief PythonFields constructs empty field
     */
    PySimFields();

    /**
     * @brief addStaticField adds static field
     * @param files files to sum up
     * @param volts voltage values to multiply up
     */
    void addStaticField(const PyList &files, const PyList &volts, double gu_mm);

    /**
     * @brief addSinField adds sinuspidal field
     * @param files same as above
     * @param volts same as above
     * @param freq_MHz field frequency
     * @param phase_rad field phase
     */
    void addSinField
    (
        const PyList &files,
        const PyList &volts,
        double freq_MHz,
        double phase_rad,
        double gu_mm
    );

    /**
     * @brief potential estimates field potential at given coordinate and time
     * @param x
     * @param y
     * @param z
     * @param t
     * @return
     */
    double potential(double x, double y, double z, double t = 0.) const;

    /**
     * @brief field estimates field strength at given coordinate and time
     * @param x
     * @param y
     * @param z
     * @param t
     * @return
     */
    PyTuple field(double x, double y, double z, double t = 0.) const;

private:
    std::shared_ptr<CompoundField> mField;

    /**
     * @brief loadFromFiles creates base field from files
     * @param files
     * @param volts
     * @return
     */
    std::unique_ptr<Field> loadFromFiles(const PyList &files, const PyList &volts);
};

#endif // PYTHONFIELDS_H

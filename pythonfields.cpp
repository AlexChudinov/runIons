#include "pythonfields.h"
#include "simionfield.h"
#include <boost/python/extract.hpp>

PyFieldException::PyFieldException(const std::string &msg)
    :
      mMsg(msg)
{

}

const char *PyFieldException::what() const throw()
{
    return mMsg.c_str();
}

void PyFieldException::pyFieldExceptTrans(const PyFieldException &e)
{
    PyErr_SetString(PyExc_Exception, e.what());
}

PySimFields::PySimFields()
    :
      mField(new CompoundField())
{

}

void PySimFields::addStaticField(const PyList& files, const PyList& volts, double gu_mm)
{
    mField->addField(new StaticField(loadFromFiles(files, volts).release()), gu_mm);
}

void PySimFields::addSinField
(
    const PySimFields::PyList &files,
    const PySimFields::PyList &volts,
    double freq_MHz,
    double phase_rad,
    double gu_mm
)
{
    std::unique_ptr<CompoundField::TimeFunction>
            timeFun(new SinFunction(freq_MHz, phase_rad));
    mField->addField
    (
        new DynamicField
        (
            loadFromFiles(files, volts).release(),
            timeFun.release()
        ),
        gu_mm
    );
}

double PySimFields::potential(double x, double y, double z, double t) const
{
    CompoundField::MayBeDouble mayBePot
            = mField->potential(CompoundField::Vector3d(x, y, z), t);
    if(mayBePot) return mayBePot;
    else return 0.0;
}

PySimFields::PyTuple PySimFields::field(double x, double y, double z, double t) const
{
    CompoundField::MayBeVector3d mayBeVec
            = mField->field(CompoundField::Vector3d(x, y, z), t);
    if(mayBeVec) return boost::python::make_tuple
    (
        mayBeVec.val().x(),
        mayBeVec.val().y(),
        mayBeVec.val().z()
    );
    else return boost::python::make_tuple(0.0, 0.0, 0.0);
}

std::unique_ptr<Field> PySimFields::loadFromFiles
(
    const PySimFields::PyList &files,
    const PySimFields::PyList &volts
)
{
    try{
        std::unique_ptr<Field> field(new SimionField);
        if(boost::python::len(files) != boost::python::len(volts))
            throw PyFieldException("addStaticField: files and volts have different lengths");
        for(int i = 0; i < boost::python::len(files); ++i)
        {
            char const * path = boost::python::extract<char const *>(files[i]);
            const double voltage = boost::python::extract<const double>(volts[i]);
            field->add(path, voltage);
        }
        return field;
    }
    catch(const std::string& e)
    {
        throw PyFieldException("addStaticField: simion code says:" + e);
    }

}

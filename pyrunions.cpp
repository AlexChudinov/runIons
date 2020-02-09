#include "pyrunions.h"
#include <boost/python/iterator.hpp>

PyRunIons::PyObserver::PyObserver(PyRunIons::PyList &list, int stepsPerSample)
    :
      mList(list),
      mStepsPerSample(stepsPerSample)
{
}

void PyRunIons::PyObserver::write(const DefaultTrackIon::State &state, double time_us)
{
    static int cnts;
    if(cnts++ >= mStepsPerSample)
    {
        PyList list;
        list.append(state[0]);
        list.append(state[1]);
        list.append(state[2]);
        list.append(state[3]);
        list.append(state[4]);
        list.append(state[5]);
        list.append(time_us);
        mList.append(list);
        cnts = 0;
    }
}

void PyRunIons::PyObserver::startWrite(const DefaultTrackIon::State &state, double time_us)
{
    PyList list;
    list.append(state[0]);
    list.append(state[1]);
    list.append(state[2]);
    list.append(state[3]);
    list.append(state[4]);
    list.append(state[5]);
    list.append(time_us);
    mList.append(list);
}

void PyRunIons::PyObserver::finalWrite(const DefaultTrackIon::State &state, double time_us)
{
    startWrite(state, time_us);
}


PyRunIons::PyStopCondition::PyStopCondition(double maxTime_us)
    :
      mMaxTime_us(maxTime_us)
{

}

bool PyRunIons::PyStopCondition::stop(const DefaultTrackIon::State &, double time_us) const
{
    return time_us >= mMaxTime_us;
}

PyRunIons::PyRunIons(PySimFields field)
    :
      mField(field),
      mTrackIon(new DefaultTrackIon)
{
}

PyRunIons::PyList PyRunIons::integrators()
{
    PyList res;
    for(const std::string& s : DefaultTrackIon::IntegratorSelector::avaliable())
        res.append(s);
    return res;
}

void PyRunIons::setIntegrator(const std::string &name, double mz_amu)
{
    if(integrators().count(name) == 0)
        throw PyFieldException("setIntegrator: no integrator with name: " + name);
    mTrackIon->setIntegrator(name, mz_amu, mField.mField);
}

void PyRunIons::initState(double x, double y, double z, double vx, double vy, double vz)
{
    mTrackIon->setInitState({x, y, z, vx, vy, vz});
}

void PyRunIons::setTime_us(double time_us)
{
    mTrackIon->setTime(time_us);
}

void PyRunIons::setTimeStep_us(double h_us)
{
    mTrackIon->setTimeStep(h_us);
}

void PyRunIons::run(PyRunIons::PyList &result, double stopTime, int stepsPerSample)
{
    mTrackIon->setObserver(new PyObserver(result, stepsPerSample));
    mTrackIon->setStopCond(new PyStopCondition(stopTime));
    mTrackIon->run();
}

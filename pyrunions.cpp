#include "pyrunions.h"
#include <boost/python/iterator.hpp>

PyRunIon::PyObserver::PyObserver(PyRunIon::PyList &list, int stepsPerSample)
    :
      mList(list),
      mStepsPerSample(stepsPerSample)
{
}

void PyRunIon::PyObserver::write(const DefaultTrackIon::State &state, double time_us)
{
    static int cnts;
    if(cnts++ >= mStepsPerSample)
    {
        startWrite(state, time_us);
        cnts = 0;
    }
}

void PyRunIon::PyObserver::startWrite(const DefaultTrackIon::State &state, double time_us)
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

void PyRunIon::PyObserver::finalWrite(const DefaultTrackIon::State &state, double time_us)
{
    startWrite(state, time_us);
}


PyRunIon::PyStopCondition::PyStopCondition(double maxTime_us)
    :
      mMaxTime_us(maxTime_us)
{

}

bool PyRunIon::PyStopCondition::stop(const DefaultTrackIon::State &, double time_us) const
{
    return time_us >= mMaxTime_us;
}

PyRunIon::PyRunIon(PySimFields field)
    :
      mField(field),
      mTrackIon(new DefaultTrackIon)
{
}

PyRunIon::PyList PyRunIon::integrators()
{
    PyList res;
    for(const std::string& s : DefaultTrackIon::IntegratorSelector::avaliable())
        res.append(s);
    return res;
}

void PyRunIon::setIntegrator(const std::string &name, double mz_amu)
{
    if(integrators().count(name) == 0)
        throw PyFieldException("setIntegrator: no integrator with name: " + name);
    mTrackIon->setIntegrator(name, mz_amu, mField.mField);
}

void PyRunIon::initState(double x, double y, double z, double vx, double vy, double vz)
{
    mTrackIon->setInitState({x, y, z, vx, vy, vz});
}

void PyRunIon::setTime_us(double time_us)
{
    mTrackIon->setTime(time_us);
}

void PyRunIon::setTimeStep_us(double h_us)
{
    mTrackIon->setTimeStep(h_us);
}

void PyRunIon::run(PyRunIon::PyList &result, double stopTime, int stepsPerSample)
{
    mTrackIon->setObserver(new PyObserver(result, stepsPerSample));
    mTrackIon->setStopCond(new PyStopCondition(stopTime));
    mTrackIon->run();
}

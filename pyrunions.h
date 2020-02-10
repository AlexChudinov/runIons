#ifndef PYRUNIONS_H
#define PYRUNIONS_H

#include "pythonfields.h"
#include "trackion.h"

class PyRunIon
{
    PySimFields mField;

    std::shared_ptr<DefaultTrackIon> mTrackIon;

public:
    using PyList = PySimFields::PyList;

    class PyObserver : public DefaultTrackIon::Observer
    {
        PyList & mList;
        int mStepsPerSample;
    public:
        PyObserver(PyList & list, int stepsPerSample);

        void write(const DefaultTrackIon::State &state, double time_us);

        void startWrite(const DefaultTrackIon::State &state, double time_us);

        void finalWrite(const DefaultTrackIon::State &state, double time_us);
    };

    class PyStopCondition : public DefaultTrackIon::TrackStop
    {
        const double mMaxTime_us;
    public:
        PyStopCondition(double maxTime_us);

        bool stop(const DefaultTrackIon::State &, double time_us) const;
    };

    PyRunIon(PySimFields field);

    /**
     * @brief integrators obtains list of available integrators
     * @return
     */
    static PyList integrators();

    void setIntegrator(const std::string& name, double mz_amu);

    void initState(double x, double y, double z, double vx, double vy, double vz);

    void setTime_us(double time_us);

    void setTimeStep_us(double h_us);

    void run(PyList& result, double stopTime, int stepsPerSample);
};

class PyRunIons
{
    std::unique_ptr<TrackIonBunch> mTracker;
public:
    PyRunIons(PySimFields field,
              const std::string& integratorName,
              int nIons,
              double x0_mm = 0.0,
              double y0_mm = 0.0,
              double z0_mm = 0.0,
              double xnorm = 0.0,
              double ynorm = 0.0,
              double znorm = 1.0,
              double tobMax_us = 1.0);

private:
};

#endif // PYRUNIONS_H

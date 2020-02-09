#include "trackion.h"

DefaultTrackIon::IntegratorSelector::MapStrFun
    DefaultTrackIon::IntegratorSelector::sIntegrators;

const DefaultTrackIon::IntegratorSelector::StringList
    &DefaultTrackIon::IntegratorSelector::avaliable()
{
    static std::vector<std::string> res;
    if(res.size() != sIntegrators.size())
    {
        res.reserve(sIntegrators.size());
        for(const auto& e : sIntegrators)
            res.push_back(e.first);
    }
    return res;
}

std::unique_ptr<PhaseStateIntegrator> DefaultTrackIon::IntegratorSelector::select
(
  const std::string &name,
  PhaseStateIntegrator::Eqns *eqns
)
{
    MapStrFun::const_iterator it = sIntegrators.find(name);
    if(it != sIntegrators.end())
        return it->second(eqns);
    else
        return std::unique_ptr<PhaseStateIntegrator>();
}

void DefaultTrackIon::IntegratorSelector::addIntegrator
(
    const std::string &name,
    DefaultTrackIon::IntegratorSelector::IntegratorCreator cons
)
{
    MapStrFun::const_iterator it = sIntegrators.find(name);
    if(it == sIntegrators.end())
    {
        sIntegrators.insert({name, cons});
    }
}

void DefaultTrackIon::IntegratorSelector::addDefaultIntegrators()
{
    addIntegrator
    (
        "modified_midpoint",
        [](PhaseStateIntegrator::Eqns * eqns)->std::unique_ptr<PhaseStateIntegrator>
    {
        return std::unique_ptr<PhaseStateIntegrator>(new ModifiedMidpoint(eqns));
    }
    );

    addIntegrator
    (
        "runge_kutta4",
        [](PhaseStateIntegrator::Eqns * eqns)->std::unique_ptr<PhaseStateIntegrator>
    {
        return std::unique_ptr<PhaseStateIntegrator>(new RK4Integrator(eqns));
    }
    );
}

DefaultTrackIon::DefaultTrackIon()
{
    IntegratorSelector::addDefaultIntegrators();
}

void DefaultTrackIon::notify()
{
    mIsRunning = false;
}

void DefaultTrackIon::run()
{
    if(mIntegrator)
    {
        mIsRunning = true;
        State S1 = mInitState, S0 = S1;
        double h = mH_us;
        double t0 = mTime_us;
        if(mObs) mObs->startWrite(S0, t0); //write first step
        if(mStopCond && mStopCond->stop(mInitState, t0))
            return; //early return if stop condition was achieved before integration
        while(h / mH_us >= 1e-10)
        {
            S1 = mIntegrator->doStep(S0, t0, h);
            if(!mIsRunning || (mStopCond && mStopCond->stop(S1, t0 + h)))
            {//step was unsuccessfull
                mIsRunning = true;
                h *= .5;
            }
            else //if succeded
            {
                S0 = S1;
                t0 += h;
                if(mObs) mObs->write(S0, t0);
            }
        }
        mObs->finalWrite(S0, t0);
    }
}

void DefaultTrackIon::setInitState(const DefaultTrackIon::State &state)
{
    mInitState = state;
}

void DefaultTrackIon::setTimeStep(double h_us)
{
    mH_us = h_us;
}

void DefaultTrackIon::setTime(double time_us)
{
    mTime_us = time_us;
}

void DefaultTrackIon::setIntegrator
(
    const std::string &name,
    double mz_amu,
    const std::shared_ptr<CompoundField> &field
)
{
    constexpr double gAmu_kg = 1.66053892173E-27;
    constexpr double gElemCharge_C = 1.60217663410E-19;
    constexpr double gMzFactor = gAmu_kg / gElemCharge_C;
    double factor = 1e-6/(mz_amu * gMzFactor);
    PhaseStateIntegrator::Eqns * eqns
            = new SimpleEqns(field, factor, this);
    mIntegrator.reset
    (
        IntegratorSelector::select
        (
            name,
            eqns
        ).release()
    );
}

void DefaultTrackIon::setStopCond(DefaultTrackIon::TrackStop *stopCond)
{
    mStopCond.reset(stopCond);
}

void DefaultTrackIon::setObserver(DefaultTrackIon::Observer *obs)
{
    mObs.reset(obs);
}

FileObserver::FileObserver(const std::string &fileName, int stepsPerSample)
    :
      mStepsPerSample(stepsPerSample)
{
    mFile.open(fileName, std::ostream::out);
}

void FileObserver::write(const DefaultTrackIon::State &state, double time_us)
{
    static int cnt;
    if(cnt++ >= mStepsPerSample)
    {
        startWrite(state, time_us);
        cnt = 0;
    }
}

void FileObserver::startWrite(const DefaultTrackIon::State &state, double time_us)
{
    mFile << time_us << "\t" << state[0] << "\t" << state[1]
          << "\t" << state[2] << "\t" << state[3] << "\t"
          << state[4] << "\t" << state[5] << "\n";
}

void FileObserver::finalWrite(const DefaultTrackIon::State &state, double time_us)
{
    startWrite(state, time_us);
}

VectorObserver::VectorObserver
(
    std::vector<DefaultTrackIon::State> &states,
    std::vector<double> &times,
    int stepsPerSample
)
    :
      mStates(states),
      mTimes(times),
      mStepsPerSample(stepsPerSample)
{

}

void VectorObserver::write(const DefaultTrackIon::State &state, double time_us)
{
    static int cnt;
    if(cnt++ >= mStepsPerSample)
    {
        startWrite(state, time_us);
        cnt = 0;
    }
}

void VectorObserver::startWrite(const DefaultTrackIon::State &state, double time_us)
{
    mStates.push_back(state);
    mTimes.push_back(time_us);
}

void VectorObserver::finalWrite(const DefaultTrackIon::State &state, double time_us)
{
    startWrite(state, time_us);
}

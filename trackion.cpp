#include "trackion.h"
#include <QtConcurrent>
#include <ctime>

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
    constexpr double gMzFactor = amu_kg / e_C;
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


TrackIonBunch::TrackIonBunch
(
    TrackIonBunch::InitStateGenerator *source,
    TrackIonBunch::InitTimeGenerator *time
)
    :
      mSource(source),
      mTime(time)
{
}

void TrackIonBunch::run()
{
    for(auto t : mTrackers)
    {
        t->setInitState(mSource->generate());
        t->setTime(mTime->generate());
    }
    QtConcurrent::map(mTrackers, [this](std::shared_ptr<DefaultTrackIon> t)
    {
        t->run();
    }).waitForFinished();
}

void TrackIonBunch::addTracker(DefaultTrackIon * tracker)
{
    mTrackers.push_back(std::shared_ptr<DefaultTrackIon>(tracker));
}

RoundSpotTemp::RoundSpotTemp
(
    double mass_amu,
    double Temp_K,
    double r0_mm,
    const RoundSpotTemp::Vector3d &pos_mm,
    const RoundSpotTemp::Vector3d &norm
)
    :
      mGen(std::time(nullptr)),
      mUniDist(0, r0_mm),
      mNormDist
      (
        0.0,
        std::sqrt(DefaultTrackIon::kB_J * Temp_K / mass_amu / DefaultTrackIon::amu_kg)/1000
      ),
      mPos_mm(pos_mm),
      mEn(r0_mm * makeOrthogonal(norm)),
      mEt(r0_mm * mEn.cross(norm).normalized())
{
}

TrackIonBunch::InitStateGenerator::State RoundSpotTemp::generate()
{
    const double phi = 2. * M_PI * mUniDist(mGen);
    const double r = std::sqrt(mUniDist(mGen));
    const Vector3d pos = (mEn * std::cos(phi) + mEt * std::sin(phi)) * r;
    return {pos.x(), pos.y(), pos.z(), mNormDist(mGen), mNormDist(mGen), mNormDist(mGen)};
}

RoundSpotTemp::Vector3d RoundSpotTemp::makeOrthogonal(const RoundSpotTemp::Vector3d &v)
{
    Vector3d rn;
    if(v.x() != 0)
    {
        rn.x() = -(v.y() + v.z())/v.x();
        rn.y() = 1;
        rn.z() = 1;
    }
    else if(v.y() != 0)
    {
        rn.x() = 1;
        rn.y() = -(v.x() + v.z())/v.y();
        rn.z() = 1;
    }
    else if(v.z() != 0)
    {
        rn.x() = 1;
        rn.y() = 1;
        rn.z() = -(v.x() + v.y())/v.z();
    }
    else
    {
        throw std::string("RoundSpotTemp::makeOrthogonal : zero division error");
    }
    return rn.normalized();
}


UniformStartTime::UniformStartTime(double maxTime_us)
    :
      mGen(std::time(nullptr)),
      mUniDist(0, maxTime_us)
{

}

double UniformStartTime::generate()
{
    return mUniDist(mGen);
}

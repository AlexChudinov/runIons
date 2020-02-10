#ifndef TRACKION_H
#define TRACKION_H

#include "timefield.h"
#include "phasestateintegrator.h"
#include <map>
#include <fstream>
#include <QRunnable>
#include <functional>

/**
 * @brief The TrackIon class interface to track ions trough the fields
 */
class DefaultTrackIon
        :
        public QRunnable,
        public PhaseStateIntegrator::Eqns::Notifier
{
public:
    using State = PhaseStateIntegrator::State;

    /**
     * Ion Constants
     */
    static constexpr double amu_kg = 1.66053892173E-27;
    static constexpr double kB_J = 1.380649e-23;
    static constexpr double e_C = 1.60217663410E-19;


    /**
     * @brief The InegratorSelector class selects integrator using string
     * lict
     */
    class IntegratorSelector
    {
    public:
        using StringList = std::vector<std::string>;
        using IntegratorCreator
            = std::function
            <
                std::unique_ptr<PhaseStateIntegrator>
                (
                    PhaseStateIntegrator::Eqns*
                )
            >;
        using MapStrFun = std::map<std::string, IntegratorCreator>;

        /**
         * @brief avaliable is getter for integrators list
         * @return list of integrators
         */
        static const StringList& avaliable();

        /**
         * @brief select creates integrator by name
         * @param name string integrator name
         * @return new integrator instance
         */
        static std::unique_ptr<PhaseStateIntegrator>
            select(const std::string& name, PhaseStateIntegrator::Eqns * eqns);

        static void addIntegrator
        (
            const std::string& name,
            IntegratorCreator cons
        );

        /**
         * @brief addDefaultIntegrators adds integrators by default
         */
        static void addDefaultIntegrators();

    private:
        /**
         * @brief sIntegrators availuable integrators
         */
        static MapStrFun sIntegrators;
    };

    /**
     * @brief The Observer class is interface to observe state
     * at each integrator step
     */
    class Observer
    {
    public:
        virtual void write(const State& state, double time_us) = 0;
        virtual void startWrite(const State& state, double time_us) = 0;
        virtual void finalWrite(const State& state, double time_us) = 0;
    };

    /**
     * @brief The TrackStop class is interface to control track duration
     */
    class TrackStop
    {
    public:
        /**
         * @brief stop checks current state
         * @param state is current state
         * @return true if track should be interrapted
         */
        virtual bool stop(const State& state, double time_us) const = 0;
    };

    DefaultTrackIon();

    virtual ~DefaultTrackIon()=default;

    /**
     * @brief notify interface to be notified from eqns
     */
    void notify();

    /**
     * @brief run runs ion using given state as initial condition
     * @param s0 init state
     */
    void run();

    //Setters:
    void setInitState(const State& state);
    void setTimeStep(double h_us);
    void setTime(double time_us);
    virtual void setIntegrator
    (
        const std::string& name,
        double mz_amu,
        const std::shared_ptr<CompoundField> &field
    );
    void setStopCond(TrackStop * stopCond);
    void setObserver(Observer * obs);

protected:
    /**
     * @brief mIsRunning flag that track can run further
     */
    bool mIsRunning;

    State mInitState;

    double mH_us;

    double mTime_us;

    std::unique_ptr<PhaseStateIntegrator> mIntegrator;

    std::unique_ptr<TrackStop> mStopCond;

    std::unique_ptr<Observer> mObs;
};

/**
 * @brief The FileObserver class implements observer which saves
 * data to file
 */
class FileObserver : public DefaultTrackIon::Observer
{
public:
    FileObserver(const std::string& fileName, int stepsPerSample = 1);

    void write(const DefaultTrackIon::State &state, double time_us);

    void startWrite(const DefaultTrackIon::State &state, double time_us);

    void finalWrite(const DefaultTrackIon::State &state, double time_us);
private:
    /**
     * @brief mFile file to save data
     */
    std::ofstream mFile;
    /**
     * @brief mStepsPerSample number of inegration time steps
     * between samples
     */
    int mStepsPerSample;
};

/**
 * @brief The FileObserver class implements observer which saves
 * data to file
 */
class VectorObserver : public DefaultTrackIon::Observer
{
public:
    VectorObserver
    (
        std::vector<DefaultTrackIon::State>& states,
        std::vector<double>& times,
        int stepsPerSample = 1
    );

    void write(const DefaultTrackIon::State &state, double time_us);

    void startWrite(const DefaultTrackIon::State &state, double time_us);

    void finalWrite(const DefaultTrackIon::State &state, double time_us);
private:
    std::vector<DefaultTrackIon::State>& mStates;
    std::vector<double>& mTimes;
    int mStepsPerSample;
};

class TrackIonBunch : public QRunnable
{
public:
    /**
     * @brief The InitStateGenerator class interface to generate init phase state
     * conditions
     */
    class InitStateGenerator
    {
    public:
        using State = PhaseStateIntegrator::State;

        virtual ~InitStateGenerator()=default;

        /**
         * @brief generate init phase state of an ion
         * @return return ion phase state
         */
        virtual State generate() = 0;
    };

    /**
     * @brief The InitTimeGenerator class interface to generate init start time
     */
    class InitTimeGenerator
    {
    public:
        virtual ~InitTimeGenerator()=default;

        virtual double generate() = 0;
    };

    TrackIonBunch(InitStateGenerator * source, InitTimeGenerator * time);

    virtual ~TrackIonBunch() = default;

    void run();

    void addTracker(DefaultTrackIon *tracker);
private:
    std::vector<std::shared_ptr<DefaultTrackIon>> mTrackers;
    std::unique_ptr<InitStateGenerator> mSource;
    std::unique_ptr<InitTimeGenerator> mTime;
};

/**
 * @brief The RoundSpotTemp class implements round spot particle source with
 * initial temperature
 */
class RoundSpotTemp : public TrackIonBunch::InitStateGenerator
{
    std::mt19937_64 mGen;
    std::uniform_real_distribution<> mUniDist;
    std::normal_distribution<> mNormDist;

public:
    using Vector3d = Eigen::Vector3d;
    RoundSpotTemp
    (
        double mass_amu,
        double Temp_K,
        double r0_mm,
        const Vector3d& pos_mm,
        const Vector3d& norm
    );

    State generate();

private:
    Vector3d makeOrthogonal(const Vector3d& v);

    const Vector3d mPos_mm;
    const Vector3d mEn;
    const Vector3d mEt;
};

class UniformStartTime : public TrackIonBunch::InitTimeGenerator
{
    std::mt19937_64 mGen;
    std::uniform_real_distribution<> mUniDist;
public:
    UniformStartTime(double maxTime_us);

    double generate();
};

#endif // TRACKION_H

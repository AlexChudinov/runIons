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
    virtual ~TrackIonBunch() = default;

    void run();
private:
    std::vector<std::shared_ptr<DefaultTrackIon>> mTrackers;
    std::unique_ptr<InitStateGenerator> mSource;
}

#endif // TRACKION_H

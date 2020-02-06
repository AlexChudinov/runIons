#ifndef PHASESTATEINTEGRATOR_H
#define PHASESTATEINTEGRATOR_H

#include <memory>
#include <array>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/modified_midpoint.hpp>
#include "maybe.h"

class CompoundField;
/**
 * @brief The PhaseStateIntegrator class integrates phase
 * state of a particle (x, y, z, vx, vy, vz)
 */
class PhaseStateIntegrator
{
public:
    using State = std::array<double, 6>;

    virtual ~PhaseStateIntegrator()=default;

    /**
     * @brief The Eqns class represents equation that will be
     * integrated
     */
    class Eqns
    {
    public:
        using State = PhaseStateIntegrator::State;

        /**
         * @brief The Notifier class use this interface if field
         * value can fail
         */
        class Notifier
        {
        public:
            virtual void notify() = 0;
        };

        virtual ~Eqns()=default;

        Eqns
        (
            const std::shared_ptr<CompoundField>& field,
            double factor,
            Notifier * notifier = nullptr
        );

        /**
         * @brief system system of differential equations that needs
         * to be solved
         * @param x system state vector
         * @param dx derivatives of state
         * @param t time
         */
        virtual void system
        (
            const State& x,
            State& dx,
            const double t_us
        ) const = 0;

    protected:
        std::weak_ptr<CompoundField> mField;

        /**
         * @brief mFactor is a coefficient between field strength and
         * acceleration of a given particle
         */
        const double mFactor;

        Notifier * mNotifier;
    };

    PhaseStateIntegrator(Eqns * eqns);

    virtual State doStep
    (
        const State& x0,
        double t_us,
        double h_us
    ) const = 0;

protected:
    std::unique_ptr<Eqns> mEqns;
};

/**
 * @brief The SimpleEqns class simplest case of motion in vacuum
 */
class SimpleEqns : public PhaseStateIntegrator::Eqns
{
public:

    SimpleEqns
    (
        const std::shared_ptr<CompoundField>& field,
        double factor,
        Notifier * notifier = nullptr
    );

    void system(const State& x0, State& dx, const double t_us) const;
};

class RK4Integrator : public PhaseStateIntegrator
{
    mutable boost::numeric::odeint::runge_kutta4<State> mStepper;

public:

    RK4Integrator(Eqns * eqns);

    State doStep(const State &x0, double t_us, double h_us) const;
};

class ModifiedMidpoint : public PhaseStateIntegrator
{
    mutable boost::numeric::odeint::modified_midpoint<State> mStepper;

public:

    ModifiedMidpoint(Eqns * eqns);

    State doStep(const State &x0, double t_us, double h_us) const;
};

#endif // PHASESTATEINTEGRATOR_H

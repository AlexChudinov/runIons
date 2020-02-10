#include "phasestateintegrator.h"
#include "timefield.h"

PhaseStateIntegrator::PhaseStateIntegrator(PhaseStateIntegrator::Eqns *eqns)
    :
      mEqns(eqns)
{

}

PhaseStateIntegrator::Eqns::Eqns
(
    const std::shared_ptr<CompoundField> &field,
    double factor,
    PhaseStateIntegrator::Eqns::Notifier *notifier
)
    :
      mField(field),
      mFactor(factor),
      mNotifier(notifier)
{

}

SimpleEqns::SimpleEqns
(
    const std::shared_ptr<CompoundField> &field,
    double factor,
    Notifier * notifier
)
    :
      PhaseStateIntegrator::Eqns(field, factor, notifier)
{

}

void SimpleEqns::system
(
    const State &x0,
    State &dx,
    const double t_us
) const
{
    std::shared_ptr<CompoundField> p;
    Field::MayBeVector3d a;
    if
    (
        (p = mField.lock())
        && (a = p->field(TimeField::Vector3d(x0[0], x0[1], x0[2]), t_us))
    )
    {
        a.val(a.val() * mFactor);
        dx = {x0[3], x0[4], x0[5], a.val().x(), a.val().y(), a.val().z()};
    }
    else
    {
        dx = {x0[3], x0[4], x0[5], 0., 0., 0.};
        if(mNotifier) mNotifier->notify();
    }
}

#define SIMPLE_INTEGRATOR_DOSTEP\
    State x1;\
    mStepper.do_step\
    (\
        [this](const State& x0, State& dx, const double t)\
    {\
        mEqns->system(x0, dx, t);\
    },\
        x0,\
        t_us,\
        x1,\
        h_us\
    );\
    return x1;
#define SIMPLE_INTEGRATOR(name)\
    name::name(PhaseStateIntegrator::Eqns *eqns)\
        :\
          PhaseStateIntegrator(eqns){}\
    name::State name::doStep(const State &x0, double t_us, double h_us) const\
    {\
        SIMPLE_INTEGRATOR_DOSTEP\
    }

SIMPLE_INTEGRATOR(RK4Integrator)

SIMPLE_INTEGRATOR(ModifiedMidpoint)


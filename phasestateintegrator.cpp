#include "phasestateintegrator.h"
#include "timefield.h"
#include <cassert>
#include <ctime>

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


const double RoundSpotTemp::kB = 1.380649e-23;
const double RoundSpotTemp::amuKg = 1.66053892173E-27;
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
      mNormDist(0.0, std::sqrt(kB * Temp_K / mass_amu / amuKg)/1000),
      mPos_mm(pos_mm),
      mEn(r0_mm * makeOrthogonal(norm)),
      mEt(r0_mm * mEn.cross(norm).normalized())
{
}

InitStateGenerator::State RoundSpotTemp::generate()
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

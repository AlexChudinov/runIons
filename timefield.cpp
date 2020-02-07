#include "timefield.h"

StaticField::StaticField(Field *field)
    :
      mField(field)
{
}

TimeField::MayBeVector3d StaticField::field(const Vector3d &pos_gu, double) const
{
    return mField->field(pos_gu);
}

TimeField::MayBeDouble StaticField::potential(const Vector3d &pos_gu, double) const
{
    return mField->potential(pos_gu);
}


DynamicField::DynamicField(Field *field, TimeField::TimeFunction *fun)
    :
      mField(field),
      mFun(fun)
{

}

TimeField::MayBeVector3d DynamicField::field(const Vector3d &pos_gu, double t_us) const
{
    MayBeVector3d t = mField->field(pos_gu);
    return t.isValid() ? MayBeVector3d(t.val() * mFun->calc(t_us)) : t;
}

TimeField::MayBeDouble DynamicField::potential(const Vector3d &pos_gu, double t_us) const
{
    MayBeDouble t = mField->potential(pos_gu);
    return t.isValid() ? MayBeDouble(t.val() * mFun->calc(t_us)) : t;
}


TimeField::MayBeVector3d CompoundField::field(const TimeField::Vector3d &pos_mm, double t_us) const
{
    MayBeVector3d res(Vector3d(0.0, 0.0, 0.0));
    res.isValid(false);
    for(size_t i = 0; i < mScaling.size(); ++i)
    {
        const TimeField::Vector3d pos_gu = pos_mm / mScaling[i];
        MayBeVector3d t = mFields[i]->field(pos_gu, t_us);
        if(t) res.val(res.val() + t.val() / mScaling[i]);
        res.isValid(res.isValid() | t.isValid());
    }
    return res;
}

TimeField::MayBeDouble CompoundField::potential
(
    const TimeField::Vector3d &pos_mm,
    double t_us
) const
{
    MayBeDouble res(0.0);
    res.isValid(false);
    for(size_t i = 0; i < mScaling.size(); ++i)
    {
        const TimeField::Vector3d pos_gu = pos_mm / mScaling[i];
        MayBeDouble t = mFields[i]->potential(pos_gu, t_us);
        if(t) res.val(res.val() + t.val());
        res.isValid(res.isValid() | t.isValid());
    }
    return res;
}

void CompoundField::addField(TimeField *field, double gu_mm)
{
    mFields.push_back(std::shared_ptr<TimeField>(field));
    mScaling.push_back(gu_mm);
}

void CompoundField::clear()
{
    mFields.clear();
    mScaling.clear();
}


SinFunction::SinFunction(double freq_MHz, double phase_rad)
    :
      mOmega_us(2. * M_PI * freq_MHz),
      mPhase_rad(phase_rad)
{

}

double SinFunction::calc(double t_us) const
{
    return std::sin(mOmega_us * t_us + mPhase_rad);
}

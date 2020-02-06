#include "simionfield.h"
#include <QVector>
#include <QtConcurrent>

//global vector to lazy parallel indexing
QVector<int> xs;
//iterate over potential array in parallel
#define ITERATE_OVER_PAR(pa, OPERATION) \
    if(xs.size() != pa.nx())\
    {\
        xs.resize(pa.nx());\
        std::iota(xs.begin(), xs.end(), 0);\
    }\
    QtConcurrent::map(xs, [&](QVector<int>::const_reference x)\
    {\
        for(int y = 0; y < pa.ny(); ++y)\
        {\
            for(int z = 0; z < pa.nz(); ++z)\
            {\
                OPERATION;\
            }\
        }\
    }).waitForFinished();\

//sequenced iteration
#define ITERATE_OVER(pa, OPERATION)\
    for(int x = 0; x < pa.nx(); ++x)\
    {\
        for(int y = 0; y < pa.ny(); ++y)\
        {\
            for(int z = 0; z < pa.nz(); ++z)\
            {\
                OPERATION;\
            }\
        }\
    }\

SimionField::SimionField()
{
}

void SimionField::add(const char *fileName, double voltage)
{
    if(!isValid())
    {
        mPotentialArray.load(fileName);
        init(voltage);
    }
    else
    {
        simion::PA pa;
        pa.load(fileName);
        const double voltage0 = electVoltage(pa);
        ITERATE_OVER_PAR
        (
            mPotentialArray,
            mPotentialArray.potential
            (
                x,
                y,
                z,
                mPotentialArray.potential(x, y, z)
                        + voltage/voltage0 * pa.potential(x, y, z)
            )
        )
    }
}

Field::MayBeDouble SimionField::potential(const Field::Vector3d &pos_gu) const
{
    if(isValid() && checkPos(pos_gu))
    {
        return MayBeDouble(mPotentialArray.potential(pos_gu.x(), pos_gu.y(), pos_gu.z()));
    }
    else
    {
        return MayBeDouble();
    }
}

Field::MayBeVector3d SimionField::field(const Field::Vector3d &pos_gu) const
{
    if(isValid() && checkPos(pos_gu))
    {
        simion::Vector3R v = mPotentialArray.field(pos_gu.x(), pos_gu.y(), pos_gu.z());
        return MayBeVector3d(Vector3d(v.x(), v.y(), v.z()));
    }
    else
    {
        return MayBeVector3d();
    }
}

bool SimionField::isValid() const
{
    return mPotentialArray.nx() > 3 && mPotentialArray.ny() > 3;
}

double SimionField::electVoltage(const simion::PA &pa)
{
    double result = 0.0;
    ITERATE_OVER
    (
        pa,
        if(pa.electrode(x, y, z) && pa.potential(x, y, z) != 0.0)
        {
            result = pa.potential(x, y, z);
            goto end;
        }
    )
    end:
    return result;
}

void SimionField::init(double voltage)
{
    const double voltage0 = electVoltage(mPotentialArray);
    ITERATE_OVER_PAR
    (
        mPotentialArray,
        mPotentialArray.potential
        (
            x, y, z,
            voltage / voltage0 * mPotentialArray.potential(x, y, z)
        );
    )
}

bool SimionField::checkPos(const Field::Vector3d pos_gu) const
{
    return
            mPotentialArray.inside(pos_gu.x(), pos_gu.y(), pos_gu.z())
            && !mPotentialArray.electrode(pos_gu.x(), pos_gu.y(), pos_gu.z());
}

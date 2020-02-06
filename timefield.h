#ifndef FIELDS_H
#define FIELDS_H

#include "field.h"
#include <vector>
#include <memory>

/**
 * @brief The TimeField class interafce for field time dependence calculation
 */
class TimeField
{
public:
    using Vector3d = Field::Vector3d;
    using MayBeDouble = Field::MayBeDouble;
    using MayBeVector3d = Field::MayBeVector3d;

    virtual ~TimeField(){;}

    /**
     * @brief The TimeFunction class incapsulates every field dependence on time
     */
    class TimeFunction
    {
    public:
        virtual ~TimeFunction(){;}

        /**
         * @brief run calculates dependent on time multiplier for a given field
         * @param t_us time (us)
         * @return field multiplier
         */
        virtual double calc(double t_us) const = 0;
    };

    /**
     * @brief field returns field at given space point and time
     * @param pos position in space
     * @param time current time
     * @return may be field strength if pos is inside the mesh
     */
    virtual MayBeVector3d field(const Vector3d& pos_gu, double time) const = 0;
    /**
     * @brief potential returns potential at given space point and time
     * @param pos position in space
     * @param time current time
     * @return may be potential value if pos is inside the mesh
     */
    virtual MayBeDouble potential(const Vector3d& pos_gu, double time) const = 0;
};

/**
 * @brief The StaticField class implements constant field
 */
class StaticField : public TimeField
{
    std::unique_ptr<Field> mField;
public:

    StaticField(Field* field);

    virtual MayBeVector3d field(const Vector3d& pos_gu, double) const;

    virtual MayBeDouble potential(const Vector3d& pos_gu, double) const;
};

/**
 * @brief The DynamicField class implements field changing in time
 */
class DynamicField : public TimeField
{
    std::unique_ptr<Field> mField;

    std::unique_ptr<TimeFunction> mFun;
public:

    DynamicField(Field* field, TimeFunction* fun);

    virtual MayBeVector3d field(const Vector3d& pos_gu, double t_us) const;

    virtual MayBeDouble potential(const Vector3d& pos_gu, double t_us) const;
};

/**
 * @brief The CompoundField class calculates the sum of fields
 */
class CompoundField : public TimeField
{
    std::vector<std::shared_ptr<TimeField>> mFields;
    /**
     * @brief mScaling is field scale units represented by grid unit in mm
     * so [mm]/mScaling = [gu]
     */
    std::vector<double> mScaling;
public:
    virtual MayBeVector3d field(const Vector3d& pos_mm, double t_us) const;

    virtual MayBeDouble potential(const Vector3d& pos_mm, double t_us) const;

    void addField(TimeField * field, double gu_mm = 1);

    void clear();
};

/**
 * @brief The SinFunction class implements sinusoidal field time function
 */
class SinFunction : public CompoundField::TimeFunction
{
    const double mOmega_us;
    const double mPhase_rad;
public:
    SinFunction(double freq_us, double phase_rad);

    virtual double calc(double t_us) const;
};
#endif // FIELDS_H

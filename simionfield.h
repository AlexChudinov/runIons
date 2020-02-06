#ifndef SIMIONFIELD_H
#define SIMIONFIELD_H

#include "field.h"
#include "pa.h"

class SimionField : public Field
{
public:
    SimionField();

    virtual void add(const char * fileName, double voltage = 1.0);

    virtual MayBeDouble potential(const Vector3d& pos_gu) const;

    virtual MayBeVector3d field(const Vector3d &pos_gu) const;

    virtual bool isValid() const;
private:

    /**
     * @brief electVoltage calculates voltage on simion electrode from file
     * @param pa electrode downloaded from file
     * @return voltage value
     */
    static double electVoltage(const simion::PA &pa);

    /**
     * @brief init initialises potential values with given voltage
     * @param voltage voltage value
     */
    void init(double voltage);

    /**
     * @brief checkPos checks that position is inside the mesh
     * @param pos_gu position to check in grid units
     * @return true if position is inside, false - otherwise
     */
    bool checkPos(const Vector3d pos_gu) const;

    simion::PA mPotentialArray;
};

#endif // SIMIONFIELD_H

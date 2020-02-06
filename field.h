#ifndef FIELD_H
#define FIELD_H

#include "maybe.h"
#include <eigen3/Eigen/Dense>

/**
 * @brief The Field class is base field interface
 */
class Field
{
public:
    using Vector3d = Eigen::Vector3d;
    using MayBeDouble = MayBe<double>;
    using MayBeVector3d = MayBe<Vector3d>;

    Field() = default;
    virtual ~Field() = default;

    /**
     * @brief add adds potentials from the file to field
     * @param fileName name of the file
     * @param voltage adding field will be multiplied by this value
     */
    virtual void add(const char * fileName, double voltage = 1.0) = 0;

    /**
     * @brief potential returns potential value at a given position
     * @param pos position in space
     * @return maybe potential value if position is inside the mesh
     */
    virtual MayBeDouble potential(const Vector3d& pos_mm) const = 0;

    /**
     * @brief field returns field at a given point
     * @param pospoint in space
     * @return maybe field if position is inside the mesh
     */
    virtual MayBeVector3d field(const Vector3d& pos_mm) const = 0;

    /**
     * @brief isValid checks the field validity
     * @return true if field is valid
     */
    virtual bool isValid() const = 0;
};

#endif // FIELD_H

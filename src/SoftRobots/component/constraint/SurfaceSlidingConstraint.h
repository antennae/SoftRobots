#pragma once

#include <sofa/core/behavior/ConstraintResolution.h>

#include <sofa/core/objectmodel/Link.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>

#include <SoftRobots/component/behavior/SoftRobotsConstraint.h>
#include <SoftRobots/component/initSoftRobots.h>

#include <sofa/core/visual/VisualParams.h>


namespace softrobots::constraint {

using sofa::core::ConstraintParams;
using sofa::core::behavior::ConstraintResolution;
using sofa::core::objectmodel::Data;
using sofa::core::topology::BaseMeshTopology;
using sofa::core::visual::VisualParams;
using sofa::helper::ReadAccessor;
using sofa::linearalgebra::BaseVector;
using sofa::type::vector;
using softrobots::behavior::SoftRobotsConstraint;
// using sofa::core::behavior::MechanicalState;


class BilateralConstraintResolution : public ConstraintResolution
{
public:
    BilateralConstraintResolution(SReal* initF=nullptr) 
        : ConstraintResolution(1)
        , _f(initF) {}
    void resolution(int line, SReal** w, SReal* d, SReal* force, SReal *dfree) override
    {
        SOFA_UNUSED(dfree);
        force[line] -= d[line] / w[line][line];
    }

    void init(int line, SReal** /*w*/, SReal* force) override
    {
        if(_f) { force[line] = *_f; }
    }

    void initForce(int line, SReal* force) override
    {
        if(_f) { force[line] = *_f; }
    }

    void store(int line, SReal* force, bool /*convergence*/) override
    {
        if(_f) *_f = force[line];
    }

protected:
    SReal* _f;
};


// // Custom constraint resolution class (add to header file)
// class SurfaceSlidingConstraintResolution : public sofa::core::behavior::ConstraintResolution
// {
// public:
//     SurfaceSlidingConstraintResolution() : ConstraintResolution(1) {}
    
//     void resolution(int line, SReal** w, SReal* d, SReal* force, SReal* dfree) override
//     {
//         // Bilateral constraint: force can be positive or negative
//         // This keeps the point exactly on the surface
//         force[line] = -dfree[line] / w[line][line];
//     }
// };




template <class DataTypes>
class SurfaceSlidingConstraint : public SoftRobotsConstraint<DataTypes> {
public:
  SOFA_CLASS(SOFA_TEMPLATE(SurfaceSlidingConstraint, DataTypes),
             SOFA_TEMPLATE(SoftRobotsConstraint, DataTypes));

  typedef typename DataTypes::VecCoord VecCoord;
  typedef typename DataTypes::VecDeriv VecDeriv;
  typedef typename DataTypes::Coord Coord;
  typedef typename DataTypes::Deriv Deriv;
  typedef typename DataTypes::MatrixDeriv MatrixDeriv;
  typedef typename Coord::value_type Real;
  typedef typename sofa::core::behavior::MechanicalState<DataTypes>
      MechanicalState;

  typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
  typedef Data<VecCoord> DataVecCoord;
  typedef Data<VecDeriv> DataVecDeriv;
  typedef Data<MatrixDeriv> DataMatrixDeriv;
  typedef sofa::type::vector<unsigned int> SetIndexArray;

  typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
  typedef sofa::core::topology::BaseMeshTopology::Quad Quad;
  typedef sofa::core::topology::BaseMeshTopology::Edge Edge;

  // Constructor
  SurfaceSlidingConstraint(
      sofa::core::behavior::MechanicalState<DataTypes> *mm = nullptr);
  ~SurfaceSlidingConstraint() override;

  ////////////////////////// Inherited from BaseObject ////////////////////
  void init() override;
  void reinit() override;
  void bwdInit() override;
  void reset() override;
  void draw(const VisualParams *vparams) override;

  // Inherited from SoftRobotsConstraint
  void buildConstraintMatrix(const sofa::core::ConstraintParams *cParams,
                             DataMatrixDeriv &c, unsigned int &cIndex,
                             const DataVecCoord &x) override;
  void getConstraintViolation(const sofa::core::ConstraintParams *cParams,
                              BaseVector *resV, const BaseVector *Jdx) override;
  void getConstraintResolution(const ConstraintParams *,
                                std::vector<ConstraintResolution *> &,
                                unsigned int &) override;
  ////////////////////////// Inherited from BaseConstraint ////////////////
  void storeLambda(const ConstraintParams *cParams,
                   sofa::core::MultiVecDerivId res,
                   const BaseVector *lambda) override;

protected:
  Data<sofa::type::vector<unsigned int>> d_pointIndex;
  Data<sofa::type::vector<Triangle>> d_triangles;

  sofa::type::vector<Edge> m_edges;
  Data<Real> m_force; // Store constraint forces
  // Distance from the point to the surface
  Data<sofa::type::vector<Real>> m_distance;

  sofa::SingleLink<SurfaceSlidingConstraint<DataTypes>, 
          sofa::core::behavior::MechanicalState<DataTypes>, 
          sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK> d_surfaceState;
  MechanicalState *m_pointState;

  ////////////////////////// Inherited attributes ////////////////////////////
  using SoftRobotsConstraint<DataTypes>::m_state;
  using SoftRobotsConstraint<DataTypes>::getContext;
  using SoftRobotsConstraint<DataTypes>::m_nbLines;
  using SoftRobotsConstraint<DataTypes>::d_constraintIndex;
  using SoftRobotsConstraint<DataTypes>::d_componentState;
  ////////////////////////////////////////////////////////////////////////////

private:
  //   void setUpData();
  void internalInit();

  void computeEdges();

  // Surface computation
  void findClosestPointOnSurface(const Coord &P, Coord &closestPoint,
                                 Deriv &normal, Real &distance);
  Coord getTriangleClosestPoint(const Coord &P, const Triangle &tri,
                                const VecCoord &positions);
  Deriv getTriangleNormal(const Triangle &tri, const VecCoord &positions);

  // Cached results
  Coord m_closestPoint;
  Deriv m_surfaceNormal;
};

} // namespace softrobots::constraint
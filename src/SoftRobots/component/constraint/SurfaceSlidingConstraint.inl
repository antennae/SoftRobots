#pragma once

#include <SoftRobots/component/constraint/SurfaceSlidingConstraint.h>

namespace softrobots::constraint {

using sofa::core::objectmodel::ComponentState;
using sofa::helper::WriteAccessor;
using sofa::type::Vec3d;
using sofa::type::vector;

template <class DataTypes>
SurfaceSlidingConstraint<DataTypes>::SurfaceSlidingConstraint(
    MechanicalState *object)
    : SoftRobotsConstraint<DataTypes>(object),
      d_pointIndex(initData(&d_pointIndex, "pointIndex",
                            "Index of the point on the surface to apply the "
                            "sliding constraint.")),
      d_surfaceState(
          initLink("surfaceState",
                   "Path to the mechanical state of the surface on which the "
                   "sliding constraint is applied.")),
      d_triangles(initData(
          &d_triangles, "triangles",
          "List of triangles on which the sliding constraint is applied. \n"
          "If no list is given, the component will \n"
          "fill the list with the context topology.")),
      m_force(initData(&m_force, "force",
                       "Constraint forces applied to the surface points.")) {
  m_force.setReadOnly(true);
}

template <class DataTypes>
SurfaceSlidingConstraint<DataTypes>::~SurfaceSlidingConstraint() {}

template <class DataTypes> void SurfaceSlidingConstraint<DataTypes>::init() {
  d_componentState.setValue(ComponentState::Valid);
  SoftRobotsConstraint<DataTypes>::init();

  internalInit();

  m_pointState =
      dynamic_cast<MechanicalState *>(getContext()->getMechanicalState());

  d_componentState.setValue(ComponentState::Valid);
}

template <class DataTypes> void SurfaceSlidingConstraint<DataTypes>::bwdInit() {
  if (d_componentState.getValue() != ComponentState::Valid)
    return;

  // Optional: Any backward initialization logic
  // Usually empty for constraints
}

template <class DataTypes> void SurfaceSlidingConstraint<DataTypes>::reinit() {
  if (d_componentState.getValue() != ComponentState::Valid)
    return;

  SoftRobotsConstraint<DataTypes>::reinit();
  internalInit();
}

template <class DataTypes> void SurfaceSlidingConstraint<DataTypes>::reset() {
  if (d_componentState.getValue() != ComponentState::Valid)
    return;
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::buildConstraintMatrix(
    const ConstraintParams *cParams, DataMatrixDeriv &cMatrix,
    unsigned int &cIndex, const DataVecCoord &x) {
  if (d_componentState.getValue() != ComponentState::Valid)
    return;

  SOFA_UNUSED(cParams);

  d_constraintIndex.setValue(cIndex);
  const auto &constraintIndex =
      sofa::helper::getReadAccessor(d_constraintIndex);

  WriteAccessor<Data<sofa::type::vector<Real>>> stored_distance =
      sofa::helper::getWriteAccessor(m_distance);

  stored_distance.clear();

  // Build constraint matrix: single row with normal direction
  MatrixDeriv &matrix = *cMatrix.beginEdit();
  MatrixDerivRowIterator rowIterator = matrix.writeLine(constraintIndex);

  const VecCoord &positions = x.getValue();
  std::vector<unsigned int> pointIndices = d_pointIndex.getValue();

  for (unsigned int i = 0; i < pointIndices.size(); ++i) {
    // Get the position of the sliding point
    const Coord &pointPos = positions[pointIndices[i]];

    // Find closest point on surface and compute normal
    Coord closestPoint;
    // normal direction to the surface
    Deriv normal;
    Real distance;
    findClosestPointOnSurface(pointPos, closestPoint, normal, distance);
    stored_distance.push_back(distance);

    // Slide point onto the surface
    rowIterator.addCol(pointIndices[i], normal);
    cIndex++;
  }

  cMatrix.endEdit();
  m_nbLines = cIndex - constraintIndex;
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::getConstraintViolation(
    const ConstraintParams *cParams, BaseVector *resV, const BaseVector *Jdx) {
  if (d_componentState.getValue() != ComponentState::Valid)
    return;

  SOFA_UNUSED(cParams);

  for (unsigned int i = 0; i < d_pointIndex.getValue().size(); ++i) {
    // Get the index of the point
    unsigned int pointIndex = d_pointIndex.getValue()[i];

    // Get the distance from the point to the surface
    Real distance = m_distance.getValue()[pointIndex];

    Real dfree = Jdx->element(0) + distance;
    // Set the constraint violation in the result vector
    resV->set(d_constraintIndex.getValue() + i, dfree);
  }
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::getConstraintResolution(
    const ConstraintParams *cParams,
    std::vector<ConstraintResolution *> &resolutions, unsigned int &cIndex) {
  if (!this->isComponentStateValid()) {
    return;
  }

  SOFA_UNUSED(cParams);

  // For each sliding point, add a bilateral constraint resolution
  for (unsigned int i = 0; i < d_pointIndex.getValue().size(); ++i) {
    // Use BilateralConstraintResolution for bilateral constraints
    // This ensures the point stays exactly on the surface (no inequality)
    resolutions[cIndex++] = new BilateralConstraintResolution();
  }
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::findClosestPointOnSurface(
    const Coord &P, Coord &closestPoint, Deriv &normal, Real &distance) {
  ReadAccessor<Data<VecCoord>> positions = d_surfaceState->readPositions();
  ReadAccessor<Data<vector<Triangle>>> triangles = d_triangles;

  Real minDist = std::numeric_limits<Real>::max();
  closestPoint = P;
  normal = Deriv(0, 0, 1); // default normal
  distance = 0;

  // Check all triangles
  for (const Triangle &tri : triangles) {
    Coord triClosest = getTriangleClosestPoint(P, tri, positions.ref());
    Real dist = (P - triClosest).norm();

    if (dist < minDist) {
      minDist = dist;
      closestPoint = triClosest;
      normal = getTriangleNormal(tri, positions.ref());

      // Determine if point is above or below surface
      Deriv toPoint = P - triClosest;
      distance = toPoint.norm();
      if (dot(toPoint, normal) < 0)
        distance = -distance;
    }
  }
}

template <class DataTypes>
typename SurfaceSlidingConstraint<DataTypes>::Coord
SurfaceSlidingConstraint<DataTypes>::getTriangleClosestPoint(
    const Coord &P, const Triangle &tri, const VecCoord &positions) {
  const Coord &A = positions[tri[0]];
  const Coord &B = positions[tri[1]];
  const Coord &C = positions[tri[2]];

  // Simple projection to triangle plane - can be improved with barycentric
  // coordinates
  Deriv AB = B - A;
  Deriv AC = C - A;
  Deriv normal = cross(AB, AC);
  normal.normalize();

  Deriv AP = P - A;
  Real projDist = dot(AP, normal);

  return P - normal * projDist;
}

template <class DataTypes>
typename SurfaceSlidingConstraint<DataTypes>::Deriv
SurfaceSlidingConstraint<DataTypes>::getTriangleNormal(
    const Triangle &tri, const VecCoord &positions) {
  const Coord &A = positions[tri[0]];
  const Coord &B = positions[tri[1]];
  const Coord &C = positions[tri[2]];

  Deriv AB = B - A;
  Deriv AC = C - A;
  Deriv normal = cross(AB, AC);
  normal.normalize();

  return normal;
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::storeLambda(
    const ConstraintParams *cParams, sofa::core::MultiVecDerivId res,
    const BaseVector *lambda) {
  if (d_componentState.getValue() != ComponentState::Valid)
    return;

  SOFA_UNUSED(cParams);
  SOFA_UNUSED(res);

  for (unsigned int i = 0; i < d_pointIndex.getValue().size(); ++i) {
    // Get the index of the point
    unsigned int pointIndex = d_pointIndex.getValue()[i];

    // Store the lambda value for this point
    m_force.setValue(lambda->element(d_constraintIndex.getValue() + i));
  }
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::internalInit() {
  if (m_state == nullptr) {
    msg_error() << "There is no mechanical state associated with this node. "
                   "The object is then deactivated. "
                   "To remove this error message, fix your scene possibly by "
                   "adding a MechanicalObject.";
    return;
  }

  /// Check that the triangles datafield contain something, otherwise
  /// get context topology
  if (d_triangles.getValue().size() == 0) {
    msg_info() << "No triangles given. Get context topology.";
    BaseMeshTopology *topology = getContext()->getMeshTopology();

    if (topology == nullptr) {
      msg_error() << "There is no topology state associated with this node. "
                     "To remove this error message, fix your scene possibly by "
                     "adding a Topology in the parent node or by giving a list "
                     "of triangles"
                     "indices as nodes parameters .";
      return;
    }

    d_triangles.setValue(topology->getTriangles());
    m_edges = topology->getEdges();
  } else {
    computeEdges();
  }

  /// Check that the triangles datafield does not contains indices that would
  /// crash the component.
  ReadAccessor<Data<VecCoord>> positions = d_surfaceState->readPositions();
  int numTris = d_triangles.getValue().size();
  auto triangles = d_triangles.getValue();
  for (int i = 0; i < numTris; i++) {
    for (int j = 0; j < 3; j++) {
      if (triangles[i][j] >= positions.size())
        msg_error() << "triangles[" << i << "][" << j << "]=" <<
        triangles[i][j]
                    << ". is too large regarding mechanicalState size of("
                    << positions.size() << ")";
    }
  }

  // check that the pointIndex is valid
  if (d_pointIndex.getValue().size() >= positions.size())
    msg_error() << "pointIndex=" << d_pointIndex.getValue()
                << " is too large regarding mechanicalState size of("
                << positions.size() << ")";
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::computeEdges() {
  ReadAccessor<Data<vector<Triangle>>> triList = d_triangles;

  std::map<Edge, unsigned int> edgeMap;
  unsigned int edgeIndex;
  m_edges.clear();

  for (Triangle t : triList) {
    std::map<Edge, unsigned int>::iterator ite;
    Edge e;
    for (unsigned int j = 0; j < 3; ++j) {
      unsigned int v1 = t[(j + 1) % 3];
      unsigned int v2 = t[(j + 2) % 3];
      // sort vertices in lexicographics order
      if (v1 < v2)
        e = Edge(v1, v2);
      else
        e = Edge(v2, v1);
      ite = edgeMap.find(e);
      if (ite == edgeMap.end()) {
        // edge not in edgeMap so create a new one
        edgeIndex = (unsigned int)m_edges.size();
        edgeMap[e] = edgeIndex;
        m_edges.push_back(e);
      }
    }
  }
}

template <class DataTypes>
void SurfaceSlidingConstraint<DataTypes>::draw(const VisualParams *vparams) {
  // const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();
  // Draw Triangles
  ReadAccessor<Data<VecCoord>> positions = d_surfaceState->readPositions();
  ReadAccessor<Data<vector<Triangle>>> triangles = d_triangles;

  std::vector<sofa::type::Vec3> pos;
  pos.reserve(triangles.size() * 3u);
  for (unsigned int i = 0; i < triangles.size(); i++) {
    const Triangle &c = triangles[i];
    const Coord &p0 = positions[c[0]];
    const Coord &p1 = positions[c[1]];
    const Coord &p2 = positions[c[2]];
    pos.emplace_back(p0[0], p0[1], p0[2]);
    pos.emplace_back(p1[0], p1[1], p1[2]);
    pos.emplace_back(p2[0], p2[1], p2[2]);
  }
  vparams->drawTool()->drawTriangles(
      pos, sofa::type::RGBAColor(0.4f, 1.0f, 0.3f, 1.0f));
}

} // namespace softrobots::constraint
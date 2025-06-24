
#define SOFTROBOTS_SURFACEDSLIDINGCONSTRAINT_CPP

#include <sofa/defaulttype/VecTypes.h>

#include <SoftRobots/component/constraint/SurfaceSlidingConstraint.inl>

#include <sofa/core/ObjectFactory.h>

namespace softrobots::constraint {

using namespace sofa::defaulttype;

void registerSurfaceSlidingConstraint(sofa::core::ObjectFactory *factory) {

  factory->registerObjects(sofa::core::ObjectRegistrationData(
                               "This component constrains a model by applying "
                               "sliding constraints on surfaces.")
                               .add<SurfaceSlidingConstraint<Vec3Types>>(true));
}

template class SOFA_SOFTROBOTS_API SurfaceSlidingConstraint<Vec3Types>;

} // namespace softrobots::constraint

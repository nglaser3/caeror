#pragma once

#include <cmath>

#include <RAJA/RAJA.hpp>
#include <camp/tuple.hpp>

#include "caeror/data.h"
#include "caeror/raja_layouts.h"

namespace caeror {
/**
 * @brief Enum for sense of a surface
 */
enum class SenseResult { Negative, On, Positive };

/**
 * @brief Intersection result container
 */
struct IntersectionResult {
  double distance_;
  bool valid_;
};

/**
 * @brief Enum of allowed axes
 */
enum Axis : int { X, Y, Z };

/**
 * @brief Enum of all allowed surface types
 */
enum class SurfaceType {
  Plane,
  AxisAlignedPlane,
  AxisAlignedCylinder,
  Sphere,
};

template <SurfaceType Type> struct SurfaceTraits;

struct SurfaceBase {
  /// @brief Unique ID associated with this surface
  const SurfaceID id_;
};

/**
 * @brief Container for surface type and underlying data to represent the
 * surface
 */
template <SurfaceType Type, typename... DataTypes>
struct Surface : public SurfaceBase {

  /// @brief Type of the surface
  static constexpr SurfaceType surface_type_ = Type;

  /// @brief Data requirements
  static constexpr size_t data_size_ = sizeof...(DataTypes) + 1;

  /// @brief Data used to define the surface
  camp::tuple<DataTypes...> data_;
};

template <SurfaceType Type>
RAJA_HOST_DEVICE SenseResult Sense(const Particle &p, const double *data);

template <SurfaceType Type>
RAJA_HOST_DEVICE IntersectionResult Intersection(const Particle &p,
                                                 const double *data);

template <SurfaceType Type>
RAJA_HOST_DEVICE bool AlignedNormal(const Particle& p, const SenseResult sense, const double *data);

// ============================================================================
// Generic Plane
// ============================================================================
/// @brief Generic plane surface
using Plane = Surface<SurfaceType::Plane, double, double, double, double>;

template <> struct SurfaceTraits<SurfaceType::Plane> {
  using Type = Plane;
};

/// @brief Helper enum for accessing generic surface data
enum class PlaneData { A, B, C, D };

template <>
RAJA_HOST_DEVICE SenseResult Sense<SurfaceType::Plane>(const Particle &p,
                                                       const double *data) {
  const double a = data[static_cast<int>(PlaneData::A)];
  const double b = data[static_cast<int>(PlaneData::B)];
  const double c = data[static_cast<int>(PlaneData::C)];
  const double d0 = data[static_cast<int>(PlaneData::D)];

  const double value = a * p.x + b * p.y + c * p.z + d0;

  if (value < 0.0)
    return SenseResult::Negative;
  if (value > 0.0)
    return SenseResult::Positive;

  return SenseResult::On;
}

template <>
RAJA_HOST_DEVICE IntersectionResult
Intersection<SurfaceType::Plane>(const Particle &p, const double *data) {
  const double a = data[static_cast<int>(PlaneData::A)];
  const double b = data[static_cast<int>(PlaneData::B)];
  const double c = data[static_cast<int>(PlaneData::C)];
  const double d0 = data[static_cast<int>(PlaneData::D)];

  const double denom = a * p.ux + b * p.uy + c * p.uz;

  if (denom == 0.0)
    return {INFINITY, false};

  const double dist = -(a * p.x + b * p.y + c * p.z + d0) / denom;

  if (dist < 0.0)
    return {INFINITY, false};

  return {dist, true};
}

template <>
RAJA_HOST_DEVICE bool
AlignedNormal<SurfaceType::Plane>(const Particle &p, const SenseResult sense, const double* data) {
  const double a = data[static_cast<int>(PlaneData::A)];
  const double b = data[static_cast<int>(PlaneData::B)];
  const double c = data[static_cast<int>(PlaneData::C)];

  double dot = a*p.ux + b*p.uy + c*p.uz;
  if (dot == 0.0) return false;
  switch (sense) {
    case SenseResult::Positive:
      return dot > 0;
    case SenseResult::Negative:
      return dot < 0;
    default:
      return false;
  }
}

// ============================================================================
// Axis Aligned Plane
// ============================================================================
/// @brief Axis aligned plane surface
using AxisAlignedPlane = Surface<SurfaceType::AxisAlignedPlane, Axis, double>;

template <> struct SurfaceTraits<SurfaceType::AxisAlignedPlane> {
  using Type = AxisAlignedPlane;
};

/// @brief Helper enum for accessing axis aligned plane data
enum class AxisAlignedPlaneData { Axis, Intercept };

template <>
RAJA_HOST_DEVICE SenseResult
Sense<SurfaceType::AxisAlignedPlane>(const Particle &p, const double *data) {
  const Axis axis =
      static_cast<Axis>(data[static_cast<int>(AxisAlignedPlaneData::Axis)]);
  const double intercept =
      data[static_cast<int>(AxisAlignedPlaneData::Intercept)];

  double value = 0.0;

  switch (axis) {
  case Axis::X:
    value = p.x - intercept;
    break;
  case Axis::Y:
    value = p.y - intercept;
    break;
  case Axis::Z:
    value = p.z - intercept;
    break;
  }

  if (value < 0.0)
    return SenseResult::Negative;
  if (value > 0.0)
    return SenseResult::Positive;

  return SenseResult::On;
}

template <>
RAJA_HOST_DEVICE IntersectionResult Intersection<SurfaceType::AxisAlignedPlane>(
    const Particle &p, const double *data) {
  const Axis axis =
      static_cast<Axis>(data[static_cast<int>(AxisAlignedPlaneData::Axis)]);
  const double intercept =
      data[static_cast<int>(AxisAlignedPlaneData::Intercept)];

  double pos{0.0}, dir{0.0};

  switch (axis) {
  case Axis::X:
    pos = p.x;
    dir = p.ux;
    break;
  case Axis::Y:
    pos = p.y;
    dir = p.uy;
    break;
  case Axis::Z:
    pos = p.z;
    dir = p.uz;
    break;
  }

  if (dir == 0.0)
    return {INFINITY, false};

  const double distance = (intercept - pos) / dir;

  if (distance < 0.0)
    return {INFINITY, false};

  return {distance, true};
}

template <>
RAJA_HOST_DEVICE bool
AlignedNormal<SurfaceType::AxisAlignedPlane>(const Particle &p, const SenseResult sense, const double* data) {
  const Axis axis =
      static_cast<Axis>(data[static_cast<int>(AxisAlignedPlaneData::Axis)]);

  double dot;
  switch (axis) {
    case Axis::X:
      dot = p.ux;
      break;
    case Axis::Y:
      dot = p.uy;
      break;
    case Axis::Z:
      dot = p.uz;
      break;
  }

  if (dot == 0.0) return false;
  switch (sense) {
    case SenseResult::Positive:
      return dot > 0;
    case SenseResult::Negative:
      return dot < 0;
    default:
      return false;
  }
}

// ============================================================================
// Axis Aligned Cylinder
// ============================================================================
/// @brief Axis aligned cylinder surface
using AxisAlignedCylinder =
    Surface<SurfaceType::AxisAlignedCylinder, Axis, double, double, double>;

template <> struct SurfaceTraits<SurfaceType::AxisAlignedCylinder> {
  using Type = AxisAlignedCylinder;
};

/// @brief Helper enum for accessing axis aligned cylinder data, centers are
/// stored (x, y), (x, z), or (y, z)
enum class AxisAlignedCylinderData { Axis, Radius, Center1, Center2 };

template <>
RAJA_HOST_DEVICE SenseResult
Sense<SurfaceType::AxisAlignedCylinder>(const Particle &p, const double *data) {
  const auto axis =
      static_cast<Axis>(data[static_cast<int>(AxisAlignedCylinderData::Axis)]);
  const auto radius = data[static_cast<int>(AxisAlignedCylinderData::Radius)];
  const auto center1 = data[static_cast<int>(AxisAlignedCylinderData::Center1)];
  const auto center2 = data[static_cast<int>(AxisAlignedCylinderData::Center2)];

  double diff1{0.0}, diff2{0.0};

  switch (axis) {
  case Axis::X:
    diff1 = center1 - p.y;
    diff2 = center2 - p.z;
    break;
  case Axis::Y:
    diff1 = center1 - p.x;
    diff2 = center2 - p.z;
    break;
  case Axis::Z:
    diff1 = center1 - p.x;
    diff2 = center2 - p.y;
    break;
  }

  double value = diff1 * diff1 + diff2 * diff2 - radius * radius;

  if (value < 0.0)
    return SenseResult::Negative;
  if (value > 0.0)
    return SenseResult::Positive;
  return SenseResult::On;
}

template <>
RAJA_HOST_DEVICE IntersectionResult
Intersection<SurfaceType::AxisAlignedCylinder>(const Particle &p,
                                               const double *data) {
  const auto axis =
      static_cast<Axis>(data[static_cast<int>(AxisAlignedCylinderData::Axis)]);
  const auto radius = data[static_cast<int>(AxisAlignedCylinderData::Radius)];
  const auto center1 = data[static_cast<int>(AxisAlignedCylinderData::Center1)];
  const auto center2 = data[static_cast<int>(AxisAlignedCylinderData::Center2)];

  double diff1{0.0}, diff2{0.0};
  double dir1{0.0}, dir2{0.0};

  switch (axis) {
  case Axis::X:
    diff1 = p.y - center1;
    diff2 = p.z - center2;
    dir1 = p.uy;
    dir2 = p.uz;
    break;
  case Axis::Y:
    diff1 = p.x - center1;
    diff2 = p.z - center2;
    dir1 = p.ux;
    dir2 = p.uz;
    break;
  case Axis::Z:
    diff1 = p.x - center1;
    diff2 = p.y - center2;
    dir1 = p.ux;
    dir2 = p.uy;
    break;
  }

  const double a = dir1 * dir1 + dir2 * dir2;
  if (a == 0.0)
    return {0.0, false};

  const double b = 2.0 * (diff1 * dir1 + diff2 * dir2);
  const double c = diff1 * diff1 + diff2 * diff2 - radius * radius;
  const double disc = b * b - 4.0 * a * c;

  if (disc < 0.0)
    return {INFINITY, false};

  const double root1 = (-b - sqrt(disc)) / (2.0 * a);
  if (root1 >= 0.0)
    return {root1, true};
  const double root2 = (-b + sqrt(disc)) / (2.0 * a);
  if (root2 >= 0.0)
    return {root2, true};

  return {INFINITY, false};
}

template <>
RAJA_HOST_DEVICE bool
AlignedNormal<SurfaceType::AxisAlignedCylinder>(const Particle &p, const SenseResult sense, const double* data) {
  const Axis axis =
      static_cast<Axis>(data[static_cast<int>(AxisAlignedCylinderData::Axis)]);
  const auto center1 = data[static_cast<int>(AxisAlignedCylinderData::Center1)];
  const auto center2 = data[static_cast<int>(AxisAlignedCylinderData::Center2)];

  double dot;
  switch (axis) {
    case Axis::X:
      dot = (p.y - center1) * p.uy + (p.z - center2) * p.uz;
      break;
    case Axis::Y:
      dot = (p.x - center1) * p.ux + (p.z - center2) * p.uz;
      break;
    case Axis::Z:
      dot = (p.x - center1) * p.ux + (p.y - center2) * p.uy;
      break;
  }

  if (dot == 0.0) return false;
  switch (sense) {
    case SenseResult::Positive:
      return dot > 0;
    case SenseResult::Negative:
      return dot < 0;
    default:
      return false;
  }
}

// ============================================================================
// Sphere
// ============================================================================
/// @brief Sphere surface
using Sphere = Surface<SurfaceType::Sphere, double, double, double, double>;

template <> struct SurfaceTraits<SurfaceType::Sphere> {
  using Type = Sphere;
};

/// @brief Helper enum for accessing sphere data
enum class SphereData { Radius, XCenter, YCenter, ZCenter };

template <>
RAJA_HOST_DEVICE SenseResult Sense<SurfaceType::Sphere>(const Particle &p,
                                                        const double *data) {
  const double radius = data[static_cast<int>(SphereData::Radius)];
  const double xc = data[static_cast<int>(SphereData::XCenter)];
  const double yc = data[static_cast<int>(SphereData::YCenter)];
  const double zc = data[static_cast<int>(SphereData::ZCenter)];

  double x = xc - p.x;
  double y = yc - p.y;
  double z = zc - p.z;

  const double value = x * x + y * y + z * z - radius * radius;

  if (value < 0.0)
    return SenseResult::Negative;
  if (value > 0.0)
    return SenseResult::Positive;
  return SenseResult::On;
}

template <>
RAJA_HOST_DEVICE IntersectionResult
Intersection<SurfaceType::Sphere>(const Particle &p, const double *data) {
  const double radius = data[static_cast<int>(SphereData::Radius)];
  const double xc = data[static_cast<int>(SphereData::XCenter)];
  const double yc = data[static_cast<int>(SphereData::YCenter)];
  const double zc = data[static_cast<int>(SphereData::ZCenter)];

  const double x = p.x - xc;
  const double y = p.y - yc;
  const double z = p.z - zc;

  const double b = x * p.ux + y * p.uy + z * p.uz;
  const double c = x * x + y * y + z * z - radius * radius;

  double disc = b * b - c;
  if (disc < 0.0)
    return {INFINITY, false};

  const double root1 = -b - sqrt(disc);
  if (root1 >= 0.0)
    return {root1, true};
  const double root2 = -b + sqrt(disc);
  if (root2 >= 0.0)
    return {root2, true};
  return {INFINITY, false};
}

template<>
RAJA_HOST_DEVICE bool
AlignedNormal<SurfaceType::Sphere>(const Particle &p, const SenseResult sense, const double* data) {
  const double xc = data[static_cast<int>(SphereData::XCenter)];
  const double yc = data[static_cast<int>(SphereData::YCenter)];
  const double zc = data[static_cast<int>(SphereData::ZCenter)];

  double dot = (p.x - xc) * p.ux + (p.y - yc) * p.uy + (p.z - zc) * p.uz;

  if (dot == 0.0) return false;
  switch (sense) {
    case SenseResult::Positive:
      return dot > 0;
    case SenseResult::Negative:
      return dot < 0;
    default:
      return false;
  }
}

// ============================================================================
// Bookkeeping
// ============================================================================
template <SurfaceType Type>
RAJA_HOST_DEVICE const auto &CastSurface(const SurfaceBase &s) {
  using FullSurfType = typename SurfaceTraits<Type>::Type;
  return static_cast<const FullSurfType &>(s);
}

template <SurfaceType Type> RAJA_HOST_DEVICE size_t SurfaceDataReq() {
  using FullSurfType = typename SurfaceTraits<Type>::Type;
  return FullSurfType::data_size_;
}

RAJA_HOST_DEVICE
SenseResult Sense(SurfaceType stype, const Particle &p, const double *data) {
  switch (stype) {
  case SurfaceType::Plane:
    return Sense<SurfaceType::Plane>(p, data);
  case SurfaceType::AxisAlignedPlane:
    return Sense<SurfaceType::AxisAlignedPlane>(p, data);
  case SurfaceType::AxisAlignedCylinder:
    return Sense<SurfaceType::AxisAlignedCylinder>(p, data);
  case SurfaceType::Sphere:
    return Sense<SurfaceType::Sphere>(p, data);
  }
}

RAJA_HOST_DEVICE
IntersectionResult Intersection(SurfaceType stype, const Particle &p,
                                const double *data) {
  switch (stype) {
  case SurfaceType::Plane:
    return Intersection<SurfaceType::Plane>(p, data);
  case SurfaceType::AxisAlignedPlane:
    return Intersection<SurfaceType::AxisAlignedPlane>(p, data);
  case SurfaceType::AxisAlignedCylinder:
    return Intersection<SurfaceType::AxisAlignedCylinder>(p, data);
  case SurfaceType::Sphere:
    return Intersection<SurfaceType::Sphere>(p, data);
  }
}

RAJA_HOST_DEVICE bool AlignedNormal(SurfaceType stype, const Particle &p, const SenseResult sense, const double* data) {
  switch (stype) {
  case SurfaceType::Plane:
    return AlignedNormal<SurfaceType::Plane>(p, sense, data);
  case SurfaceType::AxisAlignedPlane:
    return AlignedNormal<SurfaceType::AxisAlignedPlane>(p, sense, data);
  case SurfaceType::AxisAlignedCylinder:
    return AlignedNormal<SurfaceType::AxisAlignedCylinder>(p, sense, data);
  case SurfaceType::Sphere:
    return AlignedNormal<SurfaceType::Sphere>(p, sense, data);
  } 
}

} // namespace caeror
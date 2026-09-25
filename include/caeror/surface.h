#pragma once

#include <cmath>

#include <camp/tuple.hpp>
#include <RAJA/RAJA.hpp>

#include "caeror/data.h"

namespace caeror
{
/**
 * @brief Enum for sense of a surface
 */
enum class SenseResult {
  Negative,
  On,
  Positive
};

/**
 * @brief Intersection result container
 */
struct IntersectionResult{
  double distance_;
  bool valid_;
};

/**
 * @brief Enum of allowed axes
 */
enum Axis {
  X,
  Y, 
  Z
};

/**
 * @brief Enum of all allowed surface types
 */
enum class SurfaceType {
  Plane,
  AxisAlignedPlane,
  AxisAlignedCylinder,
  Sphere,
};

/**
 * @brief Container for surface type and underlying data to represent the surface
 */
template<SurfaceType Type, typename... DataTypes>
struct Surface{
  /// @brief Type of the surface
  static constexpr SurfaceType surface_type_ = Type;

  /// @brief Data used to define the surface
  camp::tuple<DataTypes...> data_;

  template <auto Index>
  RAJA_HOST_DEVICE
  decltype(auto) GetValue() const {
    return camp::get<static_cast<size_t>(Index)>(data_);
  }

  RAJA_HOST_DEVICE
  SenseResult Sense(const Point& p) const; 

  RAJA_HOST_DEVICE
  IntersectionResult Intersection(const Point& p, const Direction& d) const;
};

// ============================================================================
// Generic Plane
// ============================================================================
/// @brief Generic plane surface
using Plane = Surface<SurfaceType::Plane, double, double, double, double>;
/// @brief Helper enum for accessing generic surface data
enum class PlaneData {
  A,
  B, 
  C,
  D
};

template<>
RAJA_HOST_DEVICE
SenseResult Plane::Sense(const Point& p) const {
  const double a = GetValue<PlaneData::A>();
  const double b = GetValue<PlaneData::B>();
  const double c = GetValue<PlaneData::C>();
  const double d0 = GetValue<PlaneData::D>();

  const double value = a * p.x + b * p.y + c * p.z + d0;

  if (value < 0.0)
    return SenseResult::Negative;
  if (value > 0.0)
    return SenseResult::Positive;

  return SenseResult::On;
}

template<>
RAJA_HOST_DEVICE
IntersectionResult Plane::Intersection(const Point& p, const Direction& d) const {
  const double a = GetValue<PlaneData::A>();
  const double b = GetValue<PlaneData::B>();
  const double c = GetValue<PlaneData::C>();
  const double d0 = GetValue<PlaneData::D>();

  const double denom = a * d.x + b * d.y + c * d.z;

  if (denom == 0.0)
    return {0.0, false};

  const double dist = -(a * p.x + b * p.y + c * p.z + d0) / denom;

  if (dist < 0.0)
    return {dist, false};

  return {dist, true};
}

// ============================================================================
// Axis Aligned Plane
// ============================================================================
/// @brief Axis aligned plane surface
using AxisAlignedPlane = Surface<SurfaceType::AxisAlignedPlane, Axis, double>;
/// @brief Helper enum for accessing axis aligned plane data
enum class AxisAlignedPlaneData {
  Axis,
  Intercept
};

template<>
RAJA_HOST_DEVICE
SenseResult AxisAlignedPlane::Sense(const Point& p) const {
  const Axis axis = GetValue<AxisAlignedPlaneData::Axis>();
  const double intercept = GetValue<AxisAlignedPlaneData::Intercept>();

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
  }

  if (value < 0.0)
      return SenseResult::Negative;
  if (value > 0.0)
      return SenseResult::Positive;

  return SenseResult::On;
}

template<>
RAJA_HOST_DEVICE
IntersectionResult AxisAlignedPlane::Intersection(const Point& p, const Direction& d) const {
  const Axis axis = GetValue<AxisAlignedPlaneData::Axis>();
  const double intercept = GetValue<AxisAlignedPlaneData::Intercept>();

  double pos{0.0}, dir{0.0};

  switch (axis) {
    case Axis::X:
      pos = p.x;
      dir = d.x;
      break;
    case Axis::Y:
      pos = p.y;
      dir = d.y;
      break;
    case Axis::Z:
      pos = p.z;
      dir = d.z;
      break;
  }

  if (dir == 0.0)
    return {0.0, false};

  const double distance = (intercept - pos) / dir;

  if (distance < 0.0)
    return {distance, false};

  return {distance, true};
}

// ============================================================================
// Axis Aligned Cylinder
// ============================================================================
/// @brief Axis aligned cylinder surface
using AxisAlignedCylinder = Surface<SurfaceType::AxisAlignedCylinder, Axis, double, double, double>;
/// @brief Helper enum for accessing axis aligned cylinder data, centers are stored (x, y), (x, z), or (y, z)
enum class AxisAlignedCylinderData {
  Axis,
  Radius,
  Center1,
  Center2
};

template <>
RAJA_HOST_DEVICE
SenseResult AxisAlignedCylinder::Sense(const Point& p) const {
  const auto axis = GetValue<AxisAlignedCylinderData::Axis>();
  const auto radius = GetValue<AxisAlignedCylinderData::Radius>();
  const auto center1 = GetValue<AxisAlignedCylinderData::Center1>();
  const auto center2 = GetValue<AxisAlignedCylinderData::Center2>();

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
      diff2 = center2 - p.z;
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
RAJA_HOST_DEVICE
IntersectionResult AxisAlignedCylinder::Intersection(const Point& p, const Direction& d) const {
  const auto axis = GetValue<AxisAlignedCylinderData::Axis>();
  const auto radius = GetValue<AxisAlignedCylinderData::Radius>();
  const auto center1 = GetValue<AxisAlignedCylinderData::Center1>();
  const auto center2 = GetValue<AxisAlignedCylinderData::Center2>();

  double diff1{0.0}, diff2{0.0};
  double dir1{0.0}, dir2{0.0};

  switch (axis) {
    case Axis::X:
      diff1 = p.y - center1;
      diff2 = p.z - center2;
      dir1 = d.y;
      dir2 = d.z;
      break;
    case Axis::Y:
      diff1 = p.x - center1;
      diff2 = p.z - center2;
      dir1 = d.x;
      dir2 = d.z;
      break;
    case Axis::Z:
      diff1 = p.x - center1;
      diff2 = p.y - center2;
      dir1 = d.x;
      dir2 = d.y;
      break;
  }

  const double a = dir1 * dir1 + dir2 *dir2;
  if (a == 0.0)
    return {0.0, false};
  const double b = 2.0 * (diff1 * dir1 + diff2 * dir2);
  const double c = diff1 * diff1 + diff2 * diff2 - radius * radius;
  const double disc = b * b - 4.0 * a * c;

  if (disc < 0.0)
    return {0.0, false};

  const double root1 = (-b - sqrt(disc)) / (2.0 * a);
  if (root1 >= 0.0)
    return {root1, true};
  
  const double root2 = (-b + sqrt(disc)) / (2.0 * a);
  if (root2 >= 0.0)
    return {root2, true};

  return {0.0, false};
}

// ============================================================================
// Sphere
// ============================================================================
/// @brief Sphere surface
using Sphere = Surface<SurfaceType::Sphere, double, double, double, double>;
/// @brief Helper enum for accessing sphere data
enum class SphereData {
  Radius,
  XCenter,
  YCenter,
  ZCenter
};

template <>
RAJA_HOST_DEVICE
SenseResult Sphere::Sense(const Point& p) const {
  const double radius = GetValue<SphereData::Radius>();
  const double xc = GetValue<SphereData::XCenter>();
  const double yc = GetValue<SphereData::YCenter>();
  const double zc = GetValue<SphereData::ZCenter>();

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
RAJA_HOST_DEVICE
IntersectionResult Sphere::Intersection(const Point& p, const Direction& d) const {
  const double radius = GetValue<SphereData::Radius>();
  const double xc = GetValue<SphereData::XCenter>();
  const double yc = GetValue<SphereData::YCenter>();
  const double zc = GetValue<SphereData::ZCenter>();

  double x = p.x - xc;
  double y = p.y - yc;
  double z = p.z - zc;

  const double b = x * d.x + y * d.y + z * d.z;
  const double c = x * x + y * y + z * z - radius * radius;

  double disc = b * b - c;
  if (disc < 0.0) 
    return {0.0, false};
  
  const double root1 = -b - sqrt(disc);
  if (root1 >= 0.0)
    return {root1, true};
  const double root2 = -b + sqrt(disc);
  if (root2 >= 0.0)
    return {root2, true};
  return {0.0, false};
}
} // namespace caeror
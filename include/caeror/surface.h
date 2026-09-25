#pragma once
#include <camp/tuple.hpp>
#include <RAJA/RAJA.hpp>

#include "caeror/data.h"

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
using AxisAlignedCylinder = Surface<SurfaceType::AxisAlignedCylinder, Axis, double, double>;
/// @brief Helper enum for accessing axis aligned cylinder data, centers are stored (x, y), (x, z), or (y, z)
enum class AxisAlignedCylinderData {
  Axis,
  Center1,
  Center2
};

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
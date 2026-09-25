#pragma once
#include <camp/tuple.hpp>

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
};

/// @brief Generic plane surface
using Plane = Surface<SurfaceType::Plane, double, double, double, double>;
/// @brief Helper enum for accessing generic surface data
enum class PlaneData {
  A,
  B, 
  C,
  D
};

/// @brief Axis aligned plane surface
using AxisAlignedPlane = Surface<SurfaceType::AxisAlignedPlane, Axis, double>;
/// @brief Helper enum for accessing axis aligned plane data
enum class AxisAlignedPlaneData {
  Axis,
  Intercept
};

/// @brief Axis aligned cylinder surface
using AxisAlignedCylinder = Surface<SurfaceType::AxisAlignedCylinder, Axis, double, double>;
/// @brief Helper enum for accessing axis aligned cylinder data, centers are stored (x, y), (x, z), or (y, z)
enum class AxisAlignedCylinderData {
  Axis,
  Center1,
  Center2
};

/// @brief Sphere surface
using Sphere = Surface<SurfaceType::Sphere, double, double, double, double>;
/// @brief Helper enum for accessing sphere data
enum class SphereData {
  Radius,
  XCenter,
  YCenter,
  ZCenter
};
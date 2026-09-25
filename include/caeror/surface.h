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
/// @brief Axis aligned plane surface
using AxisAlignedPlane = Surface<SurfaceType::AxisAlignedPlane, Axis, double>;
/// @brief Axis aligned cylinder surface
using AxisAlignedCylinder = Surface<SurfaceType::AxisAlignedCylinder, Axis, double, double>;
/// @brief Sphere surface
using Sphere = Surface<SurfaceType::Sphere, double, double, double, double>;
#pragma once

#include <camp/tuple.hpp>
#include <RAJA/RAJA.hpp>

#include "caeror/geometry_constructor.h"
#include "caeror/raja_layouts.h"

namespace caeror {

  class CaerorGeometry {
    public:
      CaerorGeometry(const GeometryConstructor& input) {
        GenerateSurfacesData(input);
      }

      ~CaerorGeometry() {
        delete[] surfaces_data_;
        delete[] surface_data_offsets_;

        delete[] rpn_logic_data_;
        delete[]cell_logic_offsets_;
      }

      void GenerateSurfacesData(const GeometryConstructor& input){
        const auto& surface_types = input.GetSurfaceTypes();
        size_t surf_data_size = 0;
        surface_data_offsets_ = new size_t[surface_types.size() + 1];
        surface_data_offsets_[0] = 0;
        for (size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
          auto stype = surface_types[surf_id];
          size_t data_size;
          switch (stype) {
            case SurfaceType::Plane:
              data_size = SurfaceDataReq<SurfaceType::Plane>();
              break;
            case SurfaceType::AxisAlignedPlane:
              data_size = SurfaceDataReq<SurfaceType::AxisAlignedPlane>();
              break;
            case SurfaceType::AxisAlignedCylinder:
              data_size = SurfaceDataReq<SurfaceType::AxisAlignedCylinder>();
              break;
            case SurfaceType::Sphere:
              data_size = SurfaceDataReq<SurfaceType::Sphere>();
              break;
          }
          surface_data_offsets_[surf_id + 1] = surface_data_offsets_[surf_id] + data_size;
          surf_data_size += data_size;
        }

        surfaces_data_ = new double[surf_data_size];
        const auto& surfaces = input.GetSurfaces();
        for(size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
          auto stype = surface_types[surf_id];
          auto offset = surface_data_offsets_[surf_id];
          surfaces_data_[offset] = static_cast<double>(stype);

          const auto& surf = *surfaces[surf_id];
          switch (stype) {
            case SurfaceType::Plane: {
              FillSurfaceData<SurfaceType::Plane, PlaneData::A>(surf, offset); 
              FillSurfaceData<SurfaceType::Plane, PlaneData::B>(surf, offset); 
              FillSurfaceData<SurfaceType::Plane, PlaneData::C>(surf, offset); 
              FillSurfaceData<SurfaceType::Plane, PlaneData::D>(surf, offset); 
              break;
            }
            case SurfaceType::AxisAlignedPlane: {
              FillSurfaceData<SurfaceType::AxisAlignedPlane, AxisAlignedPlaneData::Axis>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedPlane, AxisAlignedPlaneData::Intercept>(surf, offset); 
              break;
            }
            case SurfaceType::AxisAlignedCylinder: {
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Axis>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Radius>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Center1>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Center2>(surf, offset); 
              break;
            }
            case SurfaceType::Sphere: {
              FillSurfaceData<SurfaceType::Sphere, SphereData::Radius>(surf, offset); 
              FillSurfaceData<SurfaceType::Sphere, SphereData::XCenter>(surf, offset); 
              FillSurfaceData<SurfaceType::Sphere, SphereData::YCenter>(surf, offset); 
              FillSurfaceData<SurfaceType::Sphere, SphereData::ZCenter>(surf, offset); 
              break;
            }
          }
        }
      }

      template <SurfaceType Type, auto Index>
      void FillSurfaceData(const SurfaceBase& s, size_t offset) {
        const auto& surface = CastSurface<Type>(s);
        constexpr size_t index = static_cast<size_t>(Index);
        surfaces_data_[offset + 1 + index] = static_cast<double>(camp::get<index>(surface.data_));
      }

    private:
      /// @brief 1d array for surface data and type (0th index)
      double* surfaces_data_;
      /// @brief 1d array for surface data offsets
      size_t* surface_data_offsets_;

      /// @brief 1d array for region logic
      CaerorTokenType* rpn_logic_data_;
      /// @brief 1d array for cell offsets in rpn logic
      size_t* cell_logic_offsets_;
  };
} // namespace caeror

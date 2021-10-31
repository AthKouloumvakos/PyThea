find $BUILD_DIR -type f -name "libvtkAcceleratorsVTKm-9.0.so" -exec rm -f {} \;
find $BUILD_DIR -type f -name "libvtkm_filter-9.0.a" -exec rm -f {} \;
find $BUILD_DIR -type f -name "libvtkm_cont-9.0.a" -exec rm -f {} \;
find $BUILD_DIR -type f -name "all.py" -exec sed -i -e 's/from .vtkAcceleratorsVTKm/#from .vtkAcceleratorsVTKm/g' {} \;

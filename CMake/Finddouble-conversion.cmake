find_path(double-conversion_INCLUDE_DIR
  NAMES double-conversion.h
  PATH_SUFFIXES double-conversion
  DOC "double-conversion include directory"
)

find_library(double-conversion_LIBRARY_RELEASE
  NAMES double-conversion
  PATH_SUFFIXES lib
)

find_library(double-conversion_LIBRARY_DEBUG
  NAMES double-conversion
  PATH_SUFFIXES debug/lib
)

include(SelectLibraryConfigurations)
select_library_configurations(double-conversion)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(double-conversion
  REQUIRED_VARS double-conversion_LIBRARY double-conversion_INCLUDE_DIR
)

if (double-conversion_FOUND)
  set(double-conversion_INCLUDE_DIRS "${double-conversion_INCLUDE_DIR}")
  set(double-conversion_LIBRARIES "${double-conversion_LIBRARY}")

  if (NOT TARGET double-conversion::double-conversion AND double-conversion_LIBRARY)
    add_library(double-conversion::double-conversion UNKNOWN IMPORTED)

    set_target_properties(double-conversion::double-conversion PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${double-conversion_INCLUDE_DIR}"
    )

    if (double-conversion_LIBRARY_RELEASE)
      set_property(TARGET double-conversion::double-conversion PROPERTY
        IMPORTED_LOCATION_RELEASE "${double-conversion_LIBRARY_RELEASE}")
    endif()

    if (double-conversion_LIBRARY_DEBUG)
      set_property(TARGET double-conversion::double-conversion PROPERTY
        IMPORTED_LOCATION_DEBUG "${double-conversion_LIBRARY_DEBUG}")
    endif()
  endif()
endif()
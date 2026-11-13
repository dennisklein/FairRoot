# Included by project() of a package that names the shared Geant4 targets
# (Geant4::G4clhep, ...), when Geant4 is built static only (T2b and later).
# Its config exports Geant4::G4clhep-static, ...; this hook adds the plain
# names as aliases.
find_package(Geant4 CONFIG REQUIRED GLOBAL)
get_property(_g4_targets DIRECTORY PROPERTY IMPORTED_TARGETS)
foreach(_t IN LISTS _g4_targets)
  if(_t MATCHES "^Geant4::(.+)-static$" AND NOT TARGET "Geant4::${CMAKE_MATCH_1}")
    add_library("Geant4::${CMAKE_MATCH_1}" ALIAS "${_t}")
  endif()
endforeach()

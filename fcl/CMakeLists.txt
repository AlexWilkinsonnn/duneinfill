# Install fcl files in /job subdirectory.

install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

add_subdirectory(nov2019bad)
add_subdirectory(nov2019badnoisy)
add_subdirectory(protoduneana)

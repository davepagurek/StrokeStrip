# Cornucopia library CMakeLists.txt

FILE(GLOB Cornucopia_CPP "*.cpp")
FILE(GLOB Cornucopia_H "*.h")
LIST(APPEND Cornucopia_Sources ${Cornucopia_CPP} ${Cornucopia_H})

ADD_LIBRARY( Cornucopia STATIC ${Cornucopia_Sources} )

# Eigen
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(Cornucopia Eigen3::Eigen)

INSTALL( TARGETS Cornucopia ARCHIVE DESTINATION lib )

INSTALL( DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.h"  )

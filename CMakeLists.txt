cmake_minimum_required(VERSION 3.10...3.31)
project(VascuSynth CXX)

# ITK is located via ITK_DIR or CMAKE_PREFIX_PATH (e.g. an activated conda env
# that provides libitk-devel); no machine-specific path is hard-coded here.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(VascuSynth VascuSynth.cpp SupplyMap.cpp OxygenationMap.cpp NodeTable.cpp VascularTree.cpp TreeDrawer.cpp Export.cpp)
target_link_libraries(VascuSynth ${ITK_LIBRARIES})

install(TARGETS VascuSynth RUNTIME DESTINATION bin)

# CMakeList.txt: проект CMake для shanks-transformation; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
set (CMAKE_CXX_STANDARD 17)

add_executable (shanks_transformation "main.cpp" "series.h" "shanks_transformation.h" "epsilon_algorithm.h" "test_framework.h" "test_functions.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET shanks_transformation PROPERTY CXX_STANDARD 20)
endif()

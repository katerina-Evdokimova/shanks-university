# CMakeList.txt: ������ CMake ��� shanks-transformation; �������� �������� ��� � �����������,
# ������� ����� ������ ��� ����������� �������.
#

# �������� �������� � ����������� ���� ����� �������.
add_executable (shanks_transformation "main.cpp" "series.h" "shanks_transformation.h" "epsilon_algorithm.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET shanks_transformation PROPERTY CXX_STANDARD 20)
endif()

# TODO: �������� ����� � ������� �������, ���� ��� ����������.

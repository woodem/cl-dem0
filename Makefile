# uncomment one of the following lines to enable/disable VTK build
vtkOpts = -DCLDEM_VTK -I/usr/include/vtk-5.6 -I/usr/include/vtk-5.8 -lvtkCommon -lvtkIO
# vtkOpts = 

default:
		ccache g++ -ggdb2 -O2 -std=gnu++0x $(vtkOpts) -I/usr/include/eigen3 -shared -fPIC -o _clDem.so cl/Simulation.cpp cl/Collider.cpp cl/wrapper.cpp -lboost_python -lboost_serialization -lboost_iostreams -lOpenCL -L/usr/lib64 `pkg-config python --cflags` -Wno-deprecated

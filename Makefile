default:
		g++ -std=gnu++0x -ggdb2 -lOpenCL test-chain.cc -o test-chain
		g++ -std=gnu++0x -ggdb2 -lOpenCL test-scene.cc -o test-scene
		g++ -std=gnu++0x -DMINIDEM_VTK -I/usr/include/vtk-5.6 -lvtkCommon -lvtkIO -shared -fPIC -o _miniDem.so _miniDem.cc -lboost_python -lOpenCL -L/usr/lib64 `pkg-config python --cflags` -Wno-deprecated

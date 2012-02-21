#include"Simulation.hpp"


#ifdef CLDEM_VTK
	#include<vtkSmartPointer.h>
	#include<vtkLine.h>
	#include<vtkPoints.h>
	#include<vtkPointData.h>
	#include<vtkFloatArray.h>
	#include<vtkIntArray.h>
	#include<vtkCellArray.h>
	#include<vtkCellData.h>
	#include<vtkZLibDataCompressor.h>
	#include<vtkUnstructuredGrid.h>
	#include<vtkXMLUnstructuredGridWriter.h>
	#include<vtkPolyData.h>
	#include<vtkXMLPolyDataWriter.h>
#endif

namespace clDem{
	void Simulation::writeBufs(bool wait){
		bufs.parSize=par.size(); bufs.conSize=con.size();
		clumps.resize(1); // not yet working
		bboxes.resize(par.size()*6);
		bufs.scene=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Scene),NULL);
		queue.enqueueWriteBuffer(bufs.scene,CL_FALSE,0,sizeof(Scene),&scene);
		bufs.par=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Particle)*par.size(),NULL);
		queue.enqueueWriteBuffer(bufs.par,CL_FALSE,0,sizeof(Particle)*par.size(),&(par[0]));
		bufs.con=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Contact )*con.size(),NULL);
		queue.enqueueWriteBuffer(bufs.con,CL_FALSE,0,sizeof(Contact )*con.size(),&(con[0]));
		// read-only
		bufs.clumps=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(par_id_t)*clumps.size(),NULL);
		queue.enqueueWriteBuffer(bufs.clumps,CL_FALSE,0,sizeof(par_id_t)*clumps.size(),&(clumps[0]));
		// this one could be not written at all
		bufs.bboxes=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(float)*6*par.size(),NULL);
		queue.enqueueWriteBuffer(bufs.bboxes,CL_FALSE,0,sizeof(float)*6*par.size(),&(bboxes[0]));
		if(wait) queue.finish();
	}

	void Simulation::readBufs(bool wait){
		par.resize(bufs.parSize); con.resize(bufs.conSize);
		queue.enqueueReadBuffer(bufs.scene,CL_FALSE,0,sizeof(Scene),&scene);
		queue.enqueueReadBuffer(bufs.par,CL_FALSE,0,sizeof(Particle)*par.size(),&(par[0]));
		queue.enqueueReadBuffer(bufs.con,CL_FALSE,0,sizeof(Contact )*con.size(),&(con[0]));
		// don't read clumps back
		//queue.enqueueReadBuffer(bufs.clumps,CL_FALSE,0,sizeof(par_id_t)*clumps.size(),&(clumps[0]));
		queue.enqueueReadBuffer(bufs.bboxes,CL_FALSE,0,sizeof(float)*6*par.size(),&(bboxes[0]));
		if(wait) queue.finish();
	};

	cl::Kernel Simulation::makeKernel(const char* name){
		cl::Kernel k(program,name);
		k.setArg(0,bufs.scene);
		k.setArg(1,bufs.par);
		k.setArg(2,bufs.con);
		k.setArg(3,bufs.clumps);
		k.setArg(4,bufs.bboxes);
		return k;
	}

	void Simulation::initCl(int pNum, int dNum, const string& opts){
		std::vector<cl::Platform> platforms;
		std::vector<cl::Device> devices;
		cl::Platform::get(&platforms);
		if(pNum<0){
			std::cerr<<"==================================="<<std::endl;
			for(int i=0; i<platforms.size(); i++){
				std::cerr<<i<<". platform: "<<platforms[i].getInfo<CL_PLATFORM_NAME>()<<std::endl;
				platforms[i].getDevices(CL_DEVICE_TYPE_ALL,&devices);
				for(int j=0; j<devices.size(); j++){
					std::cerr<<"\t"<<j<<". device: "<<devices[j].getInfo<CL_DEVICE_NAME>()<<std::endl;
				}
			}
			std::cerr<<"==================================="<<std::endl;
		}
		cl::Platform::get(&platforms);
		if(pNum>=(int)platforms.size()){ std::cerr<<"Only "<<platforms.size()<<" platforms available, taking 0th platform."<<std::endl; pNum=0; }
		if(pNum<0) pNum=0;
		platform=platforms[pNum];
		platforms[pNum].getDevices(CL_DEVICE_TYPE_ALL,&devices);
		if(dNum>=(int)devices.size()){ std::cerr<<"Only "<<devices.size()<<" devices available, taking 0th device."<<std::endl; dNum=0; }
		if(dNum<0) dNum=0;
		context=cl::Context(devices);
		device=devices[dNum];
		std::cerr<<"** OpenCL ready: platform \""<<platform.getInfo<CL_PLATFORM_NAME>()<<"\", device \""<<device.getInfo<CL_DEVICE_NAME>()<<"\"."<<std::endl;
		queue=cl::CommandQueue(context,device);
		// compile source
		const char* src="#include\"cl/kernels.cl\"\n\n\0";
		cl::Program::Sources source(1,std::make_pair(src,std::string(src).size()));
		program=cl::Program(context,source);
		try{
			string opts2(opts+" -I.");
			std::cerr<<"** compile with otions: "<<opts2<<endl;
			program.build(std::vector<cl::Device>({device}),opts2.c_str(),NULL,NULL);
		}catch(cl::Error& e){
			std::cerr<<"Error building source. Build log follows."<<std::endl;
			std::cerr<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)<<std::endl;
			throw std::runtime_error("Error compiling OpenCL code.");
		}
		//auto log=program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		//if(!log.empty()) std::cerr<<log<<std::endl;
		std::cerr<<"** Program compiled.\n";
	};
	void Simulation::run(int nSteps){
		if(par.empty() || con.empty()) throw std::runtime_error("There must be some particles (now "+lexical_cast<string>(par.size())+") and contacts (now "+lexical_cast<string>(con.size())+")");
		// create buffers, enqueue copies to the device
		writeBufs();
		// create kernels, set their args
		cl::Kernel stepK=makeKernel("nextTimestep_1"), integratorK=makeKernel("integrator_P"), bboxK=makeKernel("updateBboxes_P"), contactK=makeKernel("contCompute_C"), forcesK=makeKernel("forcesToParticles_C");
		for(int i=0; i<nSteps; i++){
			queue.enqueueTask(stepK);
			queue.enqueueNDRangeKernel(integratorK,cl::NDRange(0),cl::NDRange(par.size()),cl::NDRange(1));
			queue.enqueueNDRangeKernel(bboxK      ,cl::NDRange(0),cl::NDRange(par.size()),cl::NDRange(1));
			queue.enqueueNDRangeKernel(contactK   ,cl::NDRange(0),cl::NDRange(con.size()),cl::NDRange(1));
			queue.enqueueNDRangeKernel(forcesK    ,cl::NDRange(0),cl::NDRange(con.size()),cl::NDRange(1));
		}
		readBufs();
	}
	Real Simulation::pWaveDt(){
		Real ret=INFINITY;
		for(size_t i=0; i<par.size(); i++){
			const Particle& p(par[i]);
			int matId=par_matId_get(&p);
			if(matId<0 || matId>=SCENE_MAT_NUM) throw std::runtime_error("#"+lexical_cast<string>(i)+" has matId our of range 0.."+lexical_cast<string>(SCENE_MAT_NUM)+".");
			if(par_shapeT_get(&p)!=Shape_Sphere) continue;
			const Material& m=scene.materials[matId];
			if(mat_matT_get(&m)==0) throw std::runtime_error("#"+lexical_cast<string>(i)+" references void matId "+lexical_cast<string>(matId));
			Real density, young;
			switch(mat_matT_get(&m)){
				case Mat_ElastMat: density=m.mat.elast.density; young=m.mat.elast.young; break;
				default: throw std::runtime_error("Material "+lexical_cast<string>(matId)+" contains unhandled matT "+lexical_cast<string>(mat_matT_get(&m)));
			}
			ret=std::min(ret,p.shape.sphere.radius/sqrt(young/density));
		}
		return ret;
	}

	#ifdef CLDEM_VTK
		py::list Simulation::saveVtk(string prefix, bool compress, bool ascii){
			py::list savedFiles;
			string suffix="."+lexical_cast<string>(scene.step);
			vtkSmartPointer<vtkDataCompressor> compressor;
			if(compress) compressor=vtkSmartPointer<vtkZLibDataCompressor>::New();
			#define _NEW_ARR(type,name,dim) vtkSmartPointer<type> name=vtkSmartPointer<type>::New(); name->SetNumberOfComponents(dim); name->SetName(#name)
			#define _NEW_FLOAT_ARR(name,dim) _NEW_ARR(vtkFloatArray,name,dim)
			#define _NEW_INT_ARR(name,dim) _NEW_ARR(vtkFloatArray,name,dim)
			#define _INSERT_V3(name,v3) { float _foo[3]={(float)v3[0],(float)v3[1],(float)v3[2]}; name->InsertNextTupleValue(_foo); }
			/* spheres */
			{
				auto pos=vtkSmartPointer<vtkPoints>::New();
				auto cells=vtkSmartPointer<vtkCellArray>::New();
				_NEW_INT_ARR(id,1);
				_NEW_INT_ARR(matId,1);
				_NEW_FLOAT_ARR(radius,1);
				_NEW_FLOAT_ARR(mass,1);
				_NEW_FLOAT_ARR(inertia,3);
				_NEW_FLOAT_ARR(vel,3);
				_NEW_FLOAT_ARR(angVel,3);
				_NEW_FLOAT_ARR(force,3);
				_NEW_FLOAT_ARR(torque,3);
				for(size_t i=0; i<par.size(); i++){
					const Particle& p(par[i]);
					if(par_shapeT_get(&p)!=Shape_Sphere) continue;
					vtkIdType seqId=pos->InsertNextPoint(p.pos[0],p.pos[1],p.pos[2]);
					cells->InsertNextCell(1,&seqId);
					id->InsertNextValue(i);
					matId->InsertNextValue(par_matId_get(&p));
					radius->InsertNextValue(p.shape.sphere.radius);
					mass->InsertNextValue(p.mass);
					_INSERT_V3(inertia,p.inertia);
					_INSERT_V3(vel,p.vel);
					_INSERT_V3(angVel,p.angVel);
					_INSERT_V3(force,p.force);
					_INSERT_V3(torque,p.torque);
				}
				auto grid=vtkSmartPointer<vtkUnstructuredGrid>::New();
				grid->SetPoints(pos);
				grid->SetCells(VTK_VERTEX,cells);
				grid->GetPointData()->AddArray(id);
				grid->GetPointData()->AddArray(matId);
				grid->GetPointData()->AddArray(radius);
				grid->GetPointData()->AddArray(mass);
				grid->GetPointData()->AddArray(inertia);
				grid->GetPointData()->AddArray(vel);
				grid->GetPointData()->AddArray(angVel);
				grid->GetPointData()->AddArray(force);
				grid->GetPointData()->AddArray(torque);
				auto writer=vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
				if(compress) writer->SetCompressor(compressor);
				if(ascii) writer->SetDataModeToAscii();
				writer->SetFileName((prefix+".spheres"+suffix+".vtu").c_str());
				writer->SetInput(grid);
				writer->Write();
				savedFiles.append(string(writer->GetFileName()));
			}
			/* contacts */
			{
				auto parPos=vtkSmartPointer<vtkPoints>::New();
				auto cells=vtkSmartPointer<vtkCellArray>::New();
				_NEW_FLOAT_ARR(fN,1);
				_NEW_FLOAT_ARR(kN,1);
				_NEW_FLOAT_ARR(uN,1);
				//_NEW_FLOAT_ARR(absFt,1);
				//_NEW_FLOAT_ARR(kT,1);
				// particle positions are point between which contact lines are spanned
				for(const Particle& p: par) parPos->InsertNextPoint(p.pos[0],p.pos[1],p.pos[2]);
				for(const Contact& c: con){
					if(con_geomT_get(&c)==0) continue; // potential/deleted contact
					auto line=vtkSmartPointer<vtkLine>::New();
					line->GetPointIds()->SetId(0,c.ids.x);
					line->GetPointIds()->SetId(1,c.ids.y);
					cells->InsertNextCell(line);
					fN->InsertNextValue(c.force[0]); // in local coords
					// RTTI
					if(con_geomT_get(&c)==Geom_L1Geom){ uN->InsertNextValue(c.geom.l1g.uN); }
					else if(con_geomT_get(&c)==Geom_L6Geom){ uN->InsertNextValue(c.geom.l6g.uN); }
					else throw std::runtime_error("Unknown geomT=="+lexical_cast<string>(con_geomT_get(&c)));
					//
					if(con_physT_get(&c)==Phys_NormPhys){ kN->InsertNextValue(c.phys.normPhys.kN); }
					else throw std::runtime_error("Unknown physT=="+lexical_cast<string>(con_physT_get(&c)));
				}
				auto poly=vtkSmartPointer<vtkPolyData>::New();
				poly->SetPoints(parPos);
				poly->SetLines(cells);
				poly->GetCellData()->AddArray(fN);
				poly->GetCellData()->AddArray(kN);
				poly->GetCellData()->AddArray(uN);
				auto writer=vtkSmartPointer<vtkXMLPolyDataWriter>::New();
				if(compress) writer->SetCompressor(compressor);
				if(ascii) writer->SetDataModeToAscii();
				writer->SetFileName((prefix+".con"+suffix+".vtp").c_str());
				writer->SetInput(poly);
				writer->Write();
				savedFiles.append(string(writer->GetFileName()));
			}
			return savedFiles;
		}
	#endif

};

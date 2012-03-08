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
	template<typename T>
	void ensureArrayNonempty(vector<T>& arr, int arrIx, Simulation* s, const T& noItem){
		if(arr.empty()){ arr.push_back(noItem); s->scene.arrAlloc[arrIx]=1; s->scene.arrSize[arrIx]=0; }
	};

	// operates on a copy of Scene
	size_t Simulation::arrNewSize(const Scene& SS, int arrIx){
		size_t currSize=SS.arrAlloc[arrIx], reqSize=SS.arrSize[arrIx];
		if(currSize>=reqSize){
			cerr<<"?? arrIx="<<arrIx<<", currSize="<<currSize<<", reqSize="<<reqSize<<endl;
			return currSize;
		}
		size_t ret=max(currSize,max((size_t)(1.5*reqSize),reqSize+100));
		#if 1
			cerr<<"** ARR_";
			switch(arrIx){
				case ARR_CON: cerr<<"CON"; break;
				case ARR_CONFREE: cerr<<"CONFREE"; break;
				case ARR_POT: cerr<<"POT"; break;
				case ARR_POTFREE: cerr<<"CONFREE"; break;
				case ARR_CJOURNAL: cerr<<"CJOURNAL"; break;
				default: throw std::logic_error("Simulation.arrNewSize received invalid arrIx parameter "+lexical_cast<string>(arrIx));
			}
			cerr<<": "<<currSize<<"→ "<<ret<<" (≥"<<reqSize<<")."<<endl;
		#endif
		return ret;
	}
	void Simulation::writeBufs(bool wait){
		// make sure arrays are not empty
		clumps.resize(1); // not yet working
		bboxes.resize(par.size()*6,NAN);
		if(par.empty()) throw std::runtime_error("There must be some particles (now "+lexical_cast<string>(par.size())+").");
		cl_long2 no2={-1,-1};
		ensureArrayNonempty(con,ARR_CON,this,Contact());
		ensureArrayNonempty(conFree,ARR_CONFREE,this,-1);
		ensureArrayNonempty(pot,ARR_POT,this,no2);
		ensureArrayNonempty(potFree,ARR_POTFREE,this,-1);
		ensureArrayNonempty(cJournal,ARR_CJOURNAL,this,CJournalItem());
		if(scene.step==-1 && scene.rollback==0){
			// allow pre-created contacts in the very first step
			long i;
			// don't count invalid items at the end
			for(i=con.size()-1; i>=0 && con[i].ids.s0<0; i--);
			scene.arrSize[ARR_CON]=i+1;
			for(i=pot.size()-1; i>=0 && pot[i].s0<0; i--);
			scene.arrSize[ARR_POT]=i+1;
		}

		/* write actually allocated buffer sizes to scene */
		scene.arrAlloc[ARR_CON]=con.size();
		scene.arrAlloc[ARR_CONFREE]=conFree.size();
		scene.arrAlloc[ARR_POT]=pot.size();
		scene.arrAlloc[ARR_POTFREE]=potFree.size();
		scene.arrAlloc[ARR_CJOURNAL]=cJournal.size();

		// write scene
		sceneBuf=writeBuf(scene);

		// write arrays
		writeVecBuf(par,bufSize[_par]);
		// quasi-dynamic arrays
		writeVecBuf(con,bufSize[_con]);
		writeVecBuf(conFree,bufSize[_conFree]);
		writeVecBuf(pot,bufSize[_pot]);
		writeVecBuf(potFree,bufSize[_potFree]);
		writeVecBuf(cJournal,bufSize[_cJournal]);
		// fixed-size arrays
		writeVecBuf(clumps,bufSize[_clumps]);
		writeVecBuf(bboxes,bufSize[_bboxes]);
		if(wait) queue.finish();
	}

	void Simulation::readBufs(bool wait){
		readBuf(sceneBuf,scene);
		// read arrays
		readVecBuf(par,bufSize[_par]);
		readVecBuf(con,bufSize[_con]);
		readVecBuf(conFree,bufSize[_conFree]);
		readVecBuf(pot,bufSize[_pot]);
		readVecBuf(potFree,bufSize[_potFree]);
		readVecBuf(cJournal,bufSize[_cJournal]);
		//readVecBuf(clumps,bufSize[_clumps]); // not changed in the simulation
		readVecBuf(bboxes,bufSize[_bboxes]);
		if(wait) queue.finish();
	};

	cl::Kernel Simulation::makeKernel(const char* name){
		cl::Kernel k(program,name);
		k.setArg(0,sceneBuf);
		k.setArg(1,bufSize[_par].buf);
		k.setArg(2,bufSize[_con].buf);
		k.setArg(3,bufSize[_conFree].buf);
		k.setArg(4,bufSize[_pot].buf);
		k.setArg(5,bufSize[_potFree].buf);
		k.setArg(6,bufSize[_cJournal].buf);
		k.setArg(7,bufSize[_clumps].buf);
		k.setArg(8,bufSize[_bboxes].buf);
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
		const char* src="#include\"cl/kernels.cl.h\"\n\n\0";
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
	//#define LOOP_DBG(a) std::cerr<<a<<std::endl;
	#define LOOP_DBG(a)
	void Simulation::run(int _nSteps){
		// create buffers, enqueue copies to the device
		writeBufs(/*wait*/true);
		long goalStep=scene.step+_nSteps;
		/*
		TODO: create loop which will handle interrupted runs here
		and will not return until the required number of steps
		completes successfully
		*/
		int nSteps(_nSteps);
		int substepStart=0;
		const int rollbacksMax=1000;
		scene.rollback=0;
		bool allOk;
		Scene SS;
		const int maxLoop=-1;
		int loop=0;
		while(true){
			if(maxScheduledSteps>0) nSteps=std::min(nSteps,maxScheduledSteps);
			loop++;
			if(maxLoop>0 && loop>maxLoop) throw std::runtime_error("Too many interrupt loops!?");
			LOOP_DBG("scheduling ["<<nSteps<<","<<substepStart<<"]");
			// create kernels, set their args, enqueue them (non-blocking)
			runKernels(nSteps,substepStart);
			// read the scene object to a temporary object
			readBuf(sceneBuf,SS,/*wait*/true);
			LOOP_DBG("{"<<SS.step<<"}");
			// check interrupts here etc
			if(SS.interrupt.step>=0){
				LOOP_DBG(cerr<<"Interrupt: "<<(SS.interrupt.flags&INT_DESTRUCTIVE?"destructive":"non-destructive")<<", step="<<SS.interrupt.step<<", substep="<<SS.interrupt.substep<<", what="<<SS.interrupt.what);
				// non-destructive interrupts: handle them and queue remaining kernels
				if(!(SS.interrupt.flags & INT_DESTRUCTIVE)){
					/*
					non-destructive interrupt is a new rollback point; therefore:
					1. read _all_ buffers back from the device
					2. operate on *scene* rather than the temporary *SS* (since *scene* is what writeBufs writes)
					*/
					readBufs(/*wait*/true);
					if(SS.interrupt.flags & INT_BBOXES){ cerr<<"[C]"; cpuCollider->run(this); }
					else {
						throw std::runtime_error("Unknown non-destructive interrupt: flags="+lexical_cast<string>(SS.interrupt.flags));
					}
					substepStart=scene.interrupt.substep+1;
					scene.step=scene.interrupt.step;
					nSteps=goalStep-scene.interrupt.step;
					// write scene without interrupts back to the device
					Scene_interrupt_reset(&scene);
					LOOP_DBG("scene.step="<<scene.step);
					//sceneBuf=writeBuf(SS,/*wait*/true); 
					writeBufs(/*wait*/false);
				} else { // destructive: rollback device state to what we sent last time, re-run everything since then
					const float enlargeFactor=2.;
					if(SS.interrupt.flags & INT_ARRAYS){
						/* the interrupt is set just for one, but more may need attention;
						actually the information on which array to reallocate is discarded
						*/
						assert(SS.arrAlloc[ARR_CON]==con.size());
						assert(SS.arrAlloc[ARR_CONFREE]==conFree.size());
						assert(SS.arrAlloc[ARR_POT]==pot.size());
						assert(SS.arrAlloc[ARR_POTFREE]==potFree.size());
						assert(SS.arrAlloc[ARR_CJOURNAL]==cJournal.size());
						if(SS.arrAlloc[ARR_CON]<SS.arrSize[ARR_CON]){
							con.resize(arrNewSize(SS,ARR_CON));
						}
						if(SS.arrAlloc[ARR_CONFREE]<SS.arrSize[ARR_CONFREE]){
							conFree.resize(arrNewSize(SS,ARR_CONFREE),-1);
						}
						if(SS.arrAlloc[ARR_POT]<SS.arrSize[ARR_POT]){
							cl_long2 no2={-1,-1};
							pot.resize(arrNewSize(SS,ARR_POT),no2);
						}
						if(SS.arrAlloc[ARR_POTFREE]<SS.arrSize[ARR_POTFREE]){
							potFree.resize(arrNewSize(SS,ARR_POTFREE),-1);
						}
						if(SS.arrAlloc[ARR_CJOURNAL]<SS.arrSize[ARR_CJOURNAL]){
							cJournal.resize(arrNewSize(SS,ARR_CJOURNAL));
						}
					} else {
						throw std::runtime_error("Unknown destructive interrupt: flags="+lexical_cast<string>(SS.interrupt.flags));
					}
					if(SS.rollback>=rollbacksMax) throw std::runtime_error("Too many rollbacks ("+lexical_cast<string>(SS.rollback)+"), giving up.");
					cerr<<"Rollback no. "<<SS.rollback+1<<" to step "<<scene.step<<" / "<<substepStart<<endl;
					//
					nSteps=goalStep-scene.step;
					scene.rollback=SS.rollback+1;
					// substepStart=0; // use substep which was run the last time --
					writeBufs(/*wait*/false);
				}
			} else {
				if(SS.step==goalStep) break;
				if(SS.step>goalStep) abort();
				nSteps=goalStep-SS.step;
				substepStart=0;
			}
		};
		/* end interrupt loop */
		readBufs(/*wait*/true);
	}
	/* enqueue kernels for nSteps; the first step can start at step substepStart,
	which is 0 by default (beginning of the whole timestep).
	nSteps can be 0, in which case only remaining substeps will be added only
	*/
	void Simulation::runKernels(int nSteps, int substepStart){
		for(int step=0; step<nSteps+(substepStart>0?1:0); step++){
			int substepLast=-1;
			for(int substep=(step==0?substepStart:0); clDemKernels[substep].name; substep++){
				const KernelInfo& ki=clDemKernels[substep];
				if(substepLast>=ki.substep) throw std::runtime_error("Kernel substep numbers are not an increasing sequence (error in kernel.cl.h)");
				substepLast=ki.substep;
				cl::Kernel k=makeKernel(ki.name);
				switch(ki.argsType){
					case KARGS_SINGLE: queue.enqueueTask(k); break;
					case KARGS_PAR: queue.enqueueNDRangeKernel(k,cl::NDRange(0),cl::NDRange(bufSize[_par].size),cl::NDRange(1)); break;
					case KARGS_CON: queue.enqueueNDRangeKernel(k,cl::NDRange(0),cl::NDRange(bufSize[_con].size),cl::NDRange(1)); break;
					case KARGS_POT: queue.enqueueNDRangeKernel(k,cl::NDRange(0),cl::NDRange(bufSize[_pot].size),cl::NDRange(1)); break;
					default: throw std::runtime_error("Invalid KernelInfo.argsType value "+lexical_cast<string>(ki.argsType));
				}
				LOOP_DBG("{"<<ki.name<<"}");
			}
		}
	};
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

py::tuple Simulation::getBbox(par_id_t id){
	if(bboxes.size()<id*6 || id<0) throw std::runtime_error("No bbox defined for particle #"+lexical_cast<string>(id));
	return py::make_tuple(Vector3r(bboxes[id*6+0],bboxes[id*6+1],bboxes[id*6+2]),Vector3r(bboxes[id*6+3],bboxes[id*6+4],bboxes[id*6+5]));
}


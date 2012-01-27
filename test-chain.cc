
#include"scene.cl"
#include<vector>
using std::vector;

const Real r=0.005; // radius
const Real E=1e6;   // young's modulus
const Real rho=1e3; // density

int main(int argc, char* argv[]){
	if(argc<=1){
		cerr<<"Usage: "<<argv[0]<<"N support1 support2 ..."<<endl;
		return 1;
	}
	int N=atoi(argv[1])
	vector<int> supports;
	for(int i=2; i<argc; i++) supports.push_back(atoi(argv[i]));

	struct Scene scene=Scene_new();
	// scene setup
	scene.dt=.5*(r/sqrt(E/rho)); // .5 × p-Wave critical timestep
	scene.gravity={0,0,-10};
	scene.damping=.3;
	scene.materials[0]=ElastMat_new();
	scene.materials[0].young=E;

	// create chain of particles connected with each other
	vector<Particle> par(N);
	vector<Contact> con(N-1);
	for(int i=0; i<N; i++){
		Particle& p(scene->par[i]);
		p=Particle_new();
		p.pos={2*r*i,0,0};
		p.mass=rho*(4/3.)*M_PI*pow(r,3);
		Real inert=p.mass*(2/5.)*pow(r,2);
		p.inertia={inert,inert,inert};
		par_dofs_set(&p,(i==0||supports.find(i)!=supports.end()) ? 0 : par_dofs_all);
		if(i>0){ // first particle has no contact with the previous one
			scene->con[i-1]=Contact_new();
			scene->con[i-1].ids={i-1,i};
		}
	}

	/* setup OpenCL */
	/* enqueue kernels (args: scene, par, con) */
	/* gather results */
}

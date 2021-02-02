#include "std_include.h"

#include "cuda_runtime.h"
#include "cuda_profiler_api.h"


#include "DatabaseNetCDFSPV.h"
#include "Simulation.h"
#include "voronoiQuadraticEnergy.h"
#include "selfPropelledParticleDynamics.h"
#include "brownianParticleDynamics.h"
#include "analysisPackage.h"

/*!
This file compiles to produce an executable that can be used to reproduce the timing information
in the main cellGPU paper. It sets up a simulation that takes control of a voronoi model and a simple
model of active motility
NOTE that in the output, the forces and the positions are not, by default, synchronized! The NcFile
records the force from the last time "computeForces()" was called, and generally the equations of motion will 
move the positions. If you want the forces and the positions to be sync'ed, you should call the
Voronoi model's computeForces() funciton right before saving a state.
*/
int main(int argc, char*argv[])
{
    //...some default parameters
    int numpts = 256; //number of cells
    int USE_GPU = -1; //0 or greater uses a gpu, any negative number runs on the cpu
    int c;
    int tSteps = 10000; //number of time steps to run after initialization
    int initSteps = 1000; //number of initialization steps

    double dt = 0.01; //the time step size
    double p0 = 3.8;  //the preferred perimeter
    double a0 = 1.0;  // the preferred area
    double v0 = 0.1;  // the self-propulsion
    double tau = 0.00000000001;  // the memory of the cells

    int SubIntType = 2;  // parameter change (0 for none, 1 for a0 and 2 for p0)
    double SubMean = p0;  // mean value
    int SubSize = 1024;  // size of substrate grid
    double SubSigma = 0.01;  // std deviation
    int SubRadius = ((sqrt(a0)/PI)/sqrt(numpts))*SubSize;  //radius used for substrate average

    //The defaults can be overridden from the command line
    while((c=getopt(argc,argv,"n:t:g:i:e:p:a:v:s:r:")) != -1)
        switch(c)
        {
            case 'n': numpts = atoi(optarg); break;
            case 't': tSteps = atoi(optarg); break;
            case 'g': USE_GPU = atoi(optarg); break;
            case 'i': initSteps = atoi(optarg); break;
            case 'e': dt = atof(optarg); break;
            case 'p': p0 = atof(optarg); break;
            case 'a': a0 = atof(optarg); break;
            case 'v': v0 = atof(optarg); break;
            case 's': SubSigma = atof(optarg); break;
            case 'r': SubRadius = atoi(optarg); break;
            case '?':
                    if(optopt=='c')
                        std::cerr<<"Option -" << optopt << "requires an argument.\n";
                    else if(isprint(optopt))
                        std::cerr<<"Unknown option '-" << optopt << "'.\n";
                    else
                        std::cerr << "Unknown option character.\n";
                    return 1;
            default:
                       abort();
        };

    SubMean=p0;
    SubSize=1024*(sqrt(numpts)/32);
    clock_t t1,t2; //clocks for timing information
    bool reproducible = false; // if you want random numbers with a more random seed each run, set this to false
    //check to see if we should run on a GPU
    bool initializeGPU = true;
    bool gpu = chooseGPU(USE_GPU);
    if (!gpu) 
        initializeGPU = false;

    //possibly save output in netCDF format
    char dataname[256];
    sprintf(dataname,"./test_voronoi.nc");
    int Nvert = numpts;
    SPVDatabaseNetCDF ncdat(Nvert,dataname,NcFile::Replace);

    cout << "initializing a system of " << numpts << " cells at temperature " << v0 << endl;
    //define an equation of motion object...here for self-propelled cells
//    EOMPtr spp = make_shared<selfPropelledParticleDynamics>(numpts);
    shared_ptr<selfPropelledParticleDynamics> spp = make_shared<selfPropelledParticleDynamics>(numpts);
    //define a voronoi configuration with a quadratic energy functional
    shared_ptr<VoronoiQuadraticEnergy> spv  = make_shared<VoronoiQuadraticEnergy>(numpts,1.0,4.0,reproducible);

    //set the cell preferences to uniformly have A_0 = 1, P_0 = p_0
    spv->setCellPreferencesUniform(1.0,p0);
    //set the cell activity to have D_r = 1. and a given v_0
    spv->setv0Dr(v0,1.0);
    //bd->setT(v0);
    spv->setCellAdaptationTimeUniform(tau);
    spv->setSubstratePreferencesRnd(SubIntType, SubMean, SubSize, SubSigma, SubRadius);
    spv->setCellPreferencesSubstrate(a0, p0, false, SubSize);

    //combine the equation of motion and the cell configuration in a "Simulation"
    SimulationPtr sim = make_shared<Simulation>();
    sim->setConfiguration(spv);
    sim->addUpdater(spp,spv);
    //set the time step size
    sim->setIntegrationTimestep(dt);
    //initialize Hilbert-curve sorting... can be turned off by commenting out this line or seting the argument to a negative number
    //sim->setSortPeriod(initSteps/10);
    //set appropriate CPU and GPU flags
    sim->setCPUOperation(!initializeGPU);
    if (!gpu) 
        sim->setOmpThreads(abs(USE_GPU));
    sim->setReproducible(reproducible);

    /*cout<<"copying"<<endl;
    int b=0;
    char data[1000000];
    sprintf(data, "/home/diogo/CellGPU/SPV/Results/scripts1/Job_89/test_voronoi.nc");
    //sprintf(data, "./test_voronoi_SPV.nc");
    SPVDatabaseNetCDF nc(numpts,data,NcFile::ReadOnly);
    b=nc.GetNumRecs();
    nc.ReadState(spv,b-1,true);
    //exit(911);   */

    //run for a few initialization timesteps
    printf("starting initialization\n");
    ncdat.WriteState(spv); 
    for(int ii = 0; ii < initSteps; ++ii)
        {
//cout<<ii<<endl;
//ncdat.WriteState(spv);
        sim->performTimestep();
//if(ii>10)exit (911);
        };
    spv->computeGeometry();
    printf("Finished with initialization\n");
    cout << "current q = " << spv->reportq() << endl;
    //the reporting of the force should yield a number that is numerically close to zero.
    spv->reportMeanCellForce(false);
    if(!initializeGPU)
        spv->setCPU(false);//turn off globabl-ony mode

    //run for additional timesteps, compute dynamical features, and record timing information
    dynamicalFeatures dynFeat(spv->returnPositions(),spv->Box);
    logSpacedIntegers logInts(0,0.05);
    t1=clock();
    cudaProfilerStart();
    ncdat.WriteState(spv);
    for(int ii = 0; ii < tSteps; ++ii)
        {

        //if(ii%100 ==0)
        if(ii == logInts.nextSave)
            {
            //printf("timestep %i\t\t energy %f \t msd %f \t overlap %f\n",ii,spv->computeEnergy(),dynFeat.computeMSD(spv->returnPositions()),dynFeat.computeOverlapFunction(spv->returnPositions()));
            logInts.update();
            ncdat.WriteState(spv);
            };
        sim->performTimestep();
        };
    ncdat.WriteState(spv);
    cudaProfilerStop();
    t2=clock();
    printf("final state:\t\t energy %f \t msd %f \t overlap %f\n",spv->computeEnergy(),dynFeat.computeMSD(spv->returnPositions()),dynFeat.computeOverlapFunction(spv->returnPositions()));
    double steptime = (t2-t1)/(double)CLOCKS_PER_SEC/tSteps;
    cout << "timestep ~ " << steptime << " per frame; " << endl;
    spv->reportMeanCellForce(false);
    cout << spv->reportq() << endl;

    if(initializeGPU)
        cudaDeviceReset();
    return 0;
};

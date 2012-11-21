#include <lmc/lmc.h>
using namespace lmc;

int L = 50;

int main(int argc, char** argv)
{
    LatticeSquare lattice(L*L);
    RandomNumberGenerator48 rng;
    
    std::vector<double> s0(L*L, -1.0);
    SimulationIsingModel ising(lattice, rng);
    ising.setSpins(s0);
    //Tc ~ 2.3
    ising.setBeta(1.0/2.0);
    ising.setField(1.5);
    
    VisualizerXYZ xyz("traj.xyz");
    
    for (int i=0; i<500; i++) {
        xyz.visualize(ising);
        ising.run(10);
    }
}

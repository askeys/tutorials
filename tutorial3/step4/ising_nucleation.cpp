/*
    Copyright (C) 2012  Aaron S. Keys

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//c++ ising_nucleation.cpp ../step2/plugin/TpsSimulationPluginLibLMC.cpp -llmc -ltps

#include "../step2/plugin/TpsSimulationPluginLibLMC.h"
#include <lmc/lmc.h>
#include <cassert>
#include <sstream>

using namespace lmc;
using namespace std;

int L = 50;

//Define a class to compute the nucleus size:
class TpsOrderParameterNucleusSize : public TpsOrderParameter
{
    public:
        TpsOrderParameterNucleusSize()
        {
        }
        
        ~TpsOrderParameterNucleusSize()
        {
        }

        bool hA(TpsTrajectory& traj)
        {
            //nucleus size at the beginning of the trajectory:
            int n = evaluate(traj, 0);
            
            //cerr << traj.getTimeslice(0).getFilename() << "\t" << traj.getSimulationPlugin().computeHamiltonian() << "\thA: " << n << endl;
                        
            if (n < 5) {
                return true;
            }
            else {
                return false;
            }
        }

        bool hB(TpsTrajectory& traj)
        {   
            //nucleus size at the end of the trajectory:
            int n = evaluate(traj, traj.getNumberOfTimeslices()-1);
            
            //cerr << traj.getTimeslice(traj.getNumberOfTimeslices()-1).getFilename() << "\t" << traj.getSimulationPlugin().computeHamiltonian() << "\thB: " << n << endl;
            
            if (n > 100) {
                return true;
            }
            else {
                return false;
            }
        }

        double evaluate(TpsTrajectory& traj, int i)
        {
            //it takes a bit of work to get the simulation out of the trajectory
            //this is not always ncessary, but in this case our order parameter
            //is complicated and we need access to all of the variables...
            
            //either change the cast using pointers or better yet
            //use the helper function we added to change the cast:
            TpsSimulationPluginLibLMC& simplugin 
                = safeDowncastToLibLMCPlugin(traj.getSimulationPlugin());
            
            simplugin.readRestart(traj.getTimeslice(i).getFilename().c_str());

            //get a reference to the underlying Ising simulation:
            return maxCluster(simplugin.getSimulation());
        }
        
    private:

        void clusterRecursive(int i, std::vector<int>& cluster, 
            std::vector<bool>& counted, lmc::Lattice& lattice, std::vector<double>& spin)
        {
            counted[i] = true;
            cluster.push_back(i);
            std::vector<int>& nbr = lattice.getNbrs(i);

            //all spins start at -1 to start; part of cluster if 1
            for (unsigned int j=0; j<nbr.size(); j++) {
                if(!(counted[nbr[j]])) {
                    if (spin[i] > 0 && spin[nbr[j]] > 0) {
                        clusterRecursive(nbr[j], cluster, counted, lattice, spin);
                    }
                }
            }
        }

        int maxCluster(lmc::Simulation& sim)
        {
            lmc::Lattice& lattice = sim.getLattice();
            int ncell = lattice.getN();
            vector<bool> counted(ncell, false);
            std::vector<double>& spin = sim.getSpins();
            
            int nmax = -1;
			for(int i=0; i<ncell; i++) {
				if(!counted[i]) {
					vector<int> cluster_i;
					clusterRecursive(i, cluster_i, counted, lattice, spin);
					int ni = cluster_i.size();
                    if (ni > nmax) {
                        nmax = ni;
                    }
				}
			}
            return nmax;
        }

};

int NEQUIL=20;
int LTRAJ = 100;
int STEP = 1;
int NTRAJ = 1000;

int main(int argc, char** argv)
{
    LatticeSquare lattice(L*L);
    RandomNumberGenerator48 rng48;
    
    std::vector<double> s0(L*L, -1.0);
    SimulationIsingModel ising(lattice, rng48);
    //ising.setNMovesPerCycle(10);
    ising.setSpins(s0);
    ising.setBeta(1.0/2.0);
    ising.setField(1.5);
    
    TpsSimulationPluginLibLMC sim(ising);
    
	TpsTrajectoryUniformStep traj(sim, LTRAJ, STEP);
	TpsInitializerBruteForce init;
	TpsTrajectoryEnsemble tpe(traj, 0, false);

    TpsOrderParameterNucleusSize op;
	TpsRNG48 rng(1, 2, 3);    
    
	cerr << "Initializing...\n";
    //TpsAlgorithmTPS tps(tpe, rng, op, init);
    TpsAlgorithmTPS tps(tpe, rng, op, init);
	cerr << "Sampling...\n";
	
	TpsTrialMoveShotForward fshot;
	TpsTrialMoveShotBackward bshot;	
    TpsTrialMoveShift shift(0, LTRAJ, STEP);
    
	tps.addTrialMove(fshot, 0.5);
	tps.addTrialMove(bshot, 0.5);
    tps.addTrialMove(shift, 1.0);

    //********new:**************
    
    TpsAnalysisCommitter committer(op, LTRAJ/5);
    committer.setUseShortcut(true);
    TpsAnalysisTransitionStates tse(committer, TpsAnalysisTransitionStates::STEPWISE);
    tse.setAlpha(2.0);
    tse.setNmax(50);
    
    lmc::VisualizerXYZ xyz("transition_states.xyz");
    
	for (int i=0; i<NTRAJ; i++) {
		cout << i << "\n";
		tps.doStep();
        //trajectories are created in order, starting with index 0;
        //we can erase the old ones to save memory
        
        if (i%50 == 0) {
            tse.analyze(tpe.getLastTrajectory());
            std::vector<TpsTimeslice>& ts = tse.getTransitionStates();
            cerr << "found " << ts.size() << " transition states\n";
            for (int j=0; j<ts.size(); j++) {
                sim.readRestart(ts[j].getFilename().c_str());
                xyz.visualize(ising);
            }
            ts.clear();
        }
        tpe.eraseTrajectories(0, tpe.getLastTrajectory().getID()-2);
    }
}

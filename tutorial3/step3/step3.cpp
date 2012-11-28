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
//c++ -O3 step3.cpp TpsSimulationPluginLibLMC.cpp -llmc -ltps

#include "TpsSimulationPluginLibLMC.h"
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
            
            //nucleus must be smaller than N=5 to be in basin A
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
                        
            //nucleus must be larger than N=200 to be in basin B
            if (n > 200) {
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


//Define a class for visualizing trajectories:
class MyTpsAlgorithm : public TpsAlgorithmTPS
{
    public: 
        MyTpsAlgorithm(TpsTrajectoryEnsemble& t, 
            TpsRNG& r, TpsOrderParameter& f, TpsInitializer& i)
        :	TpsAlgorithmTPS(t, r, f, i)
        {
            _count = 0;
        }
        
        //overwrite this function to visualize the path after every 10 steps:
        void callbackOnTrialMoveAccepted()
        {
            if (_count++ % 10 == 0) {
                ostringstream filename;
                filename << "traj_" << _count++ << ".xyz";
                lmc::VisualizerXYZ xyz(filename.str().c_str());

                TpsTrajectory& traj = _trajectory_factory.getLastTrajectory();
                TpsSimulationPluginLibLMC& simplugin 
                    = safeDowncastToLibLMCPlugin(traj.getSimulationPlugin());
                                    
                for (int i=0; i<traj.getNumberOfTimeslices(); i++) {                    
                    simplugin.readRestart(traj.getTimeslice(i).getFilename().c_str());
                    xyz.visualize(simplugin.getSimulation());
                }
            }        
        }
        
    private:
        int _count;

};

int main(int argc, char** argv)
{
    //Initialize simulation engine:
    LatticeSquare lattice(L*L);
    RandomNumberGenerator48 rng48;    
    std::vector<double> s0(L*L, -1.0);
    SimulationIsingModel ising(lattice, rng48);
    ising.setSpins(s0);
    ising.setBeta(1.0/2.0);
    ising.setField(1.5);
    
    //Plug engine into TPS framework:
    TpsSimulationPluginLibLMC sim(ising);
    
    //Declare the type of trajectory: (100 steps total separated by 1 MC cycle) 
    int trajectory_length = 100;
    int step_size = 1;
	TpsTrajectoryUniformStep traj(sim, trajectory_length, step_size);

    //A class to get the first trajectory (just wait until an event happens) 
	TpsInitializerBruteForce init;

    //A class to hold the ensemble of reactive trajectories (TPE) 
	TpsTrajectoryEnsemble tpe(traj, 0, false);

    //A class to compute the order parameter and define different basins
    TpsOrderParameterNucleusSize op;
	
    //Random number generator for TPS
    TpsRNG48 rng(1, 2, 3);    
    
    //Initialize the first trajectory:
	cerr << "Initializing...\n";
    MyTpsAlgorithm tps(tpe, rng, op, init);
	cerr << "Sampling...\n";
	
    //Do shooting moves and shifting moves:
	TpsTrialMoveShotForward fshot;
	TpsTrialMoveShotBackward bshot;	
    TpsTrialMoveShift shift(0, trajectory_length, step_size);
    
	tps.addTrialMove(fshot, 0.5);
	tps.addTrialMove(bshot, 0.5);
    tps.addTrialMove(shift, 1.0);
    
    //********new:**************
    
    //Compute the committer.  The settings are explained in libtps documentation
    TpsAnalysisCommitter committer(op, trajectory_length/5);
    committer.setUseShortcut(true);
    TpsAnalysisTransitionStates tse(committer, TpsAnalysisTransitionStates::STEPWISE);
    tse.setAlpha(2.0);
    tse.setNmax(50);
    
    lmc::VisualizerXYZ xyz("transition_states.xyz");    
    int number_of_trajectories = 1000;
	for (int i=0; i<number_of_trajectories; i++) {
		cout << i << "\n";
		tps.doStep();        
        //every 50 steps, look for transition states:
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
        //trajectories are created in order, starting with index 0;
        //we can erase the old ones to save memory
        tpe.eraseTrajectories(0, tpe.getLastTrajectory().getID()-2);
    }
}

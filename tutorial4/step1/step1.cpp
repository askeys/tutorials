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
//c++ -O3 step1.cpp -llmc
//./a.out
//python mbar-umbrella-sampling.py

#include <lmc/lmc.h>
#include <cassert>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace lmc;
using namespace std;

#define L 50                                                                    ///<Side length of the square lattice 

//A recursive clustering algorithm, for identifying nuclei
void cluster_recursive(int i, vector<int>& cluster, 
    vector<bool>& counted, lmc::Lattice& lattice, vector<double>& spin)
{
    counted[i] = true;
    cluster.push_back(i);
    vector<int>& nbr = lattice.getNbrs(i);

    //all spins start at -1 to start; part of cluster if 1
    for (unsigned int j=0; j<nbr.size(); j++) {
        if(!(counted[nbr[j]])) {
            if (spin[i] > 0 && spin[nbr[j]] > 0) {
                cluster_recursive(nbr[j], cluster, counted, lattice, spin);
            }
        }
    }
}

//Identifies the lattice sites belonging to the largest cluster in the system
vector<int> max_cluster(lmc::Simulation& sim)
{
    lmc::Lattice& lattice = sim.getLattice();
    int ncell = lattice.getN();
    vector<bool> counted(ncell, false);
    vector<double>& spin = sim.getSpins();
    vector<int> nucleus;
    
    int nmax = -1;
    for(int i=0; i<ncell; i++) {
        if(!counted[i]) {
            vector<int> cluster_i;
            cluster_recursive(i, cluster_i, counted, lattice, spin);
            int ni = cluster_i.size();
            if (ni > nmax) {
                nmax = ni;
                nucleus = cluster_i;
            }
        }
    }
    return nucleus;
}

//A class to perform umbrella sampling over a window
class Umbrella
{    
    protected:

        //Umbrella sampling variables:
        
        double _E;                                                              ///< The energy associated with the fictitious umbrella potential
        double _N;                                                              ///< Value of the order parameter (N=nucleus size)
        double _N0;                                                             ///< Ideal value of the order parameter (minimum of the potential)
        double _kappa;                                                          ///< 1/2 * the spring constant for umbrella potential
        vector<double> _x;                                                      ///< Microstate (list of Ising spins)        

        //Data collection variables:
        
        ofstream _ofile;                                                        ///< File to write umbrella data to
        bool _equilibrate;                                                      ///< Do not write data if equilibrating

        //Simulation variables:

        LatticeSquare _lattice;                                                 ///< The lattice the simulation will run on
        RandomNumberGenerator48 _rng;                                           ///< Random number generator for the simulation
        SimulationIsingModel _sim;                                              ///< The simulation

        //run control variables:

        int _step_size;                                                         ///< Number of MC sweeps to take for each umbrella potential move attempt 
        int _steps;                                                             ///< Total number of umbrella potential moves attempted 

    public:

        //Constructor takes a prototype lattice and a simulation algorithm
        Umbrella() : _lattice(L*L), _sim(_lattice, _rng)
        {    
            //set up some initial values (see step 0):
            _sim.setBeta(1.0/2.0);
            _sim.setField(1.5);
            _x = vector<double>(L*L, -1.0);
            _sim.setSpins(_x);
            _E = 100000.0;
            _N = -10;
            _steps = 0;
            _step_size = 5;
            _equilibrate = false;
        }

        //Opens an output file that will contain a time series for N
        void openFile()
        {
            ostringstream os;
            os << "window_N=" << _N0 << ".txt";
            _ofile.open(os.str().c_str(), ios::app);
            if (_ofile.fail()) {
                cerr << "ERROR: Failed to open window " << _N0 << ".\n";
                exit(0);
            }        
        }

        //Closes the output file (we will need this in the later steps)
        void closeFile()
        {
            _ofile.close();
        }
        
        //Set the spring constant and the minimum of the umbrella potential
        void setParam(double kappa, double N0)
        {
            _kappa = kappa;
            _N0 = N0;
        }
        
        //Run nstep umbrella sweeps
        void run(int nstep)
        {
            for (int i=0; i<nstep; i++) {
                
                //for each umbrella step run _step_size steps
                _sim.run(_step_size);
                
                //compute the size of the largest cluster
                vector<int> nucleus = max_cluster(_sim);
                int N_new = nucleus.size();
                
                //compute the umbrella energy
                double E_new = _kappa*(N_new-_N0)*(N_new-_N0);

                //accept or reject according to the Metroplis criterion
                if (E_new <= _E) {
                    _N = N_new;
                    _E = E_new;
                    _x = _sim.getSpins();
                }
                else if (exp(-(E_new-_E)) < drand48()) {                    
                    _sim.setSpins(_x);
                }
                else {
                    _N = N_new;
                    _E = E_new;
                    _x = _sim.getSpins();
                }
            }
            
            //write the size of the nucleus, if not equilibrating
            if (!_equilibrate) {
                _ofile << _steps++ << "\t" << _N << endl;
            }

        }
        
        //Set the value of the equilibrate variable (determines whether or 
        //  not to write data).
        void setEquilibrate(bool e)
        {
            _equilibrate = e;
        }
        
        //Get the current value of the umbrella energy
        double getE()
        {
            return _E;
        }
        
        //Get the current size of the nucleus
        double getN()
        {
            return _N;
        }
        
        //Get the current tartge nucleus size
        double getN0()
        {
            return _N0;
        }
        
        //Get the current value of 1/2 x the spring constant
        double getKappa()
        {
            return _kappa;
        }
        
        //Set the number of MC sweeps per umbrella move
        void setStepSize(int nstep)
        {
            _step_size = nstep;
        }
};


int main(int argc, char** argv)
{
    int nsamples=10000;                                                         ///< Number of samples per window
    int nwindows=8;                                                             ///< Number of umbrella windows
    int N0min=0;                                                                ///< Minumum value of the target nucleus size
    int N0max=24;                                                               ///< Maximum value of the target nucleus size
    int N0step=(N0max-N0min)/nwindows;                                          ///< Steps between values of N0
    double k_spring = 0.4;                                                      ///< Umbrella spring constant

    //Write umbrella information to this file, for use with pymbar
    ofstream os("umbrellas");
    for (int i=0; i<nwindows; i++) {
        double N0i = N0min + i*N0step;        
        os << N0i<< "\t" << k_spring << "\n";
    }
    os.close();
    
    //Use a single Umbrella for each window by sliding the bias potential
    Umbrella window;
    for (int i=0; i<nwindows; i++) {
        //update the bias
        double N0 = N0min + i*N0step;
        
        //equilibrate the window at the new value of the bias
        cerr << "Equilibrating N0 = " << N0 << endl;
        window.setParam(0.5*k_spring, N0);
        window.setEquilibrate(true);
        window.run(100);
        window.setEquilibrate(false);
        window.openFile();    
        
        //sample the window
        cerr << "Sampling N0 = " << N0 << endl;
        for (int ii=0; ii<nsamples; ii++) {
            window.run(1);
        }
        window.closeFile();
    }
}

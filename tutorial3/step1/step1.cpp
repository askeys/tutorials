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
//c++ lmc_energy_test.cpp ../plugin/TpsSimulationPluginLibLMC.cpp -llmc -ltps

#include "TpsSimulationPluginLibLMC.h"
#include <lmc/lmc.h>
#include <cassert>

using namespace lmc;
using namespace std;

int L = 50;

int main(int argc, char** argv)
{
    //Initialize simulation engine:
    LatticeSquare lattice(L*L);
    RandomNumberGenerator48 rng;    
    std::vector<double> s0(L*L, -1.0);
    SimulationIsingModel ising(lattice, rng);
    ising.setSpins(s0);
    ising.setBeta(1.0/2.0);
    ising.setField(1.5);
    
    //Plug engine into TPS framework:
    TpsSimulationPluginLibLMC sim(ising);

    //Make sure the simulation is reading and writing restart data correctly:
    sim.writeRestart("start.dat");
    double ei = sim.computeHamiltonian();
    sim.run(1000);
    double em = sim.computeHamiltonian();
    sim.readRestart("start.dat");
    double ee = sim.computeHamiltonian();    
    //(Energies should be the same before and after restarting):
    assert(em != ei);
    assert(ei == -4*L*L + 1.5*L*L);
    assert(ee == ei);
    
    //Make sure simulation is working properly through the TPS interface:
    sim.run(10000);
    double e = 0.0;
    for (int i=0; i<100; i++) {
        sim.run(100);
        e += sim.computeHamiltonian();
    }
    e/=100;
    //The energy should belower now:
    assert(e < -4*L*L);
}

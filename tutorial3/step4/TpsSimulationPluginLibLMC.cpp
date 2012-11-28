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

#include "TpsSimulationPluginLibLMC.h"
#include <iostream>
#include <fstream>
#include <cassert>

TpsSimulationPluginLibLMC::TpsSimulationPluginLibLMC(lmc::Simulation& sim)
    : _sim(sim)
{	
}

void TpsSimulationPluginLibLMC::run(int nsteps)
{
    _sim.run(nsteps);
	callbackOnRunCommand();
}

int TpsSimulationPluginLibLMC::getDimensions()
{
	return _sim.getLattice().getDimensions();
}

void TpsSimulationPluginLibLMC::writeRestart(const char* filename)
{
	std::string timeslice_name = filename;
    
    std::vector<double> microstate = _sim.getSpins();
	_saved_timeslices[timeslice_name] = microstate;
}

void TpsSimulationPluginLibLMC::readRestart(const char* filename)
{			
	std::string timeslice_name = filename;
    std::vector<double>& microstate = _saved_timeslices[timeslice_name];
    _sim.setSpins(microstate);
}

void TpsSimulationPluginLibLMC::freeRestart(const char* filename)
{
	std::string timeslice_name = filename;
	_saved_timeslices.erase(timeslice_name);
}

void TpsSimulationPluginLibLMC::copyRestart(const char* src_file, const char* dest_file)
{
	std::string name1 = src_file;
	std::string name2 = dest_file;
	_saved_timeslices[name2] = _saved_timeslices[name1];
}

double TpsSimulationPluginLibLMC::computeHamiltonian()
{
	return _sim.computeEnergy();
}

lmc::Simulation& TpsSimulationPluginLibLMC::getSimulation()
{
    return _sim;
}

TpsSimulationPluginLibLMC& safeDowncastToLibLMCPlugin(TpsSimulationPlugin& sim)
{
	try {
		TpsSimulationPluginLibLMC& sim_LibLMC 
			= dynamic_cast<TpsSimulationPluginLibLMC&>(sim);
			return sim_LibLMC;
	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		int cast_to_LibLMC_simulation = 0;
		assert(cast_to_LibLMC_simulation);
	}
}

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
/*
 *  TpsTestRandomWalk.h
 *  tpsapi
 *
 *  Created by askeys on 6/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef TPSSIMULATIONPLUGINLIBLMC_H
#define TPSSIMULATIONPLUGINLIBLMC_H

#include <tps/tps.h>
#include <lmc/lmc.h>

#include <string>
#include <cmath>
#include <vector>
#include <map>

/**
\brief uses LibLMC as the underlying simulation engine
\ingroup extensions
\ingroup sim
*/
class TpsSimulationPluginLibLMC : public TpsSimulationPlugin
{
	public:	

		//re-implemented:
        TpsSimulationPluginLibLMC(lmc::Simulation&);
		void run(int);
		void writeRestart(const char*);
		void readRestart(const char*);
		void freeRestart(const char*);		
		void copyRestart(const char*, const char*);		
		double computeHamiltonian();
        
        //optional:
		//void writeData(const char* filename, int mode=0); 
        void invertVelocities();
        
        //extra:
        int getDimensions();
        lmc::Simulation& getSimulation();
    
	protected:
        lmc::Simulation& _sim;
		std::map< std::string, std::vector<double> > _saved_timeslices;						
};

TpsSimulationPluginLibLMC& safeDowncastToLibLMCPlugin(TpsSimulationPlugin&);

#endif

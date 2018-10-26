#include <string>
#include <iostream>
#include <exception>
#include <vector>
#include "Pythia8/Pythia.h"

using namespace Pythia8;

struct item{
	int id;
	double M;
	double e, px, py, pz;
	double weight;
	std::vector<item> radlist;
	friend std::ostream& operator<<(std::ostream& os, const item& p){
    		os << p.id << " " << p.M << " " << p.e << " " << p.px
                   << " " << p.py << " " << p.pz << " " << p.weight;
		return os;
  	}
};

int main(int argc, char* argv[]){
	if (argc < 3){
		std::cout << "need a pythia setting file and number of event" << std::endl;
		exit(-1);
	}
	std::string filename(argv[1]);
	int Neve = atoi(argv[2]);
	
	// read settings
	Pythia pythia;
	pythia.readFile(filename);
	Event& event = pythia.event;
	Info& info = pythia.info;

	// suppress output
	pythia.readString("SoftQCD:all = off");
	pythia.readString("PromptPhoton:all=off");
  	pythia.readString("WeakSingleBoson:all=off");
  	pythia.readString("WeakDoubleBoson:all=off");

	pythia.readString("Init:showProcesses = off");  
	pythia.readString("Init:showMultipartonInteractions = off");  
	pythia.readString("Init:showChangedSettings = off");  
	pythia.readString("Init:showChangedParticleData = off");  
	pythia.readString("Next:numberCount = 1000");  
	pythia.readString("Next:numberShowInfo = 0");  
	pythia.readString("Next:numberShowProcess = 0");  
	pythia.readString("Next:numberShowEvent = 0"); 

	pythia.init();

	std::vector<item> plist;
	for (int iEvent = 0; iEvent < Neve; ++iEvent) 
    	{
	    	pythia.next();
		double weight = pythia.info.weight();
		for (size_t i = 0; i < event.size(); ++i) 
		{
			auto p = event[i];
			if (p.isFinal() && p.idAbs() == 4 && std::abs(p.eta()) < 3.) 
			{
				item c_entry;
				// final momenta 
				c_entry.M = p.m(); // mass
				c_entry.id = 4; // charm quark
				c_entry.weight = weight;
				c_entry.e = p.e();
				c_entry.px = p.px();
				c_entry.py = p.py();
				c_entry.pz = p.pz();
				c_entry.radlist.clear();
				// trace back to its first production to 
				// find out final state gluon radiation
				while(true){
					auto m1 = p.mother1();
					auto m2 = p.mother2();
					auto d1 = p.daughter1();
					auto d2 = p.daughter2();

					// trace the fermion line back, else break
					if (event[m1].id() == p.id()) 
						p = event[m1];
					else if (event[m2].id() == p.id()) 
						p = event[m2];
					else break;
					
					int gluon_eid = -1;
					if (event[p.daughter1()].idAbs() == 21 
					&& std::abs(event[p.daughter1()].status()) == 51)
						gluon_eid = p.daughter1();
						
					else if (event[p.daughter2()].idAbs() == 21 
					&& std::abs(event[p.daughter2()].status()) == 51) 
						gluon_eid = p.daughter2();

					if (gluon_eid >0) { // a radiated gluon:
						// make it pre-formed gluon
						auto g = event[gluon_eid];
						c_entry.px += g.px();
						c_entry.py += g.py();
						c_entry.pz += g.pz();
						c_entry.e = std::sqrt(c_entry.px*c_entry.px + c_entry.py*c_entry.py + c_entry.pz*c_entry.pz + c_entry.M*c_entry.M);
						item gluon;
						gluon.id = 21; 
						gluon.M = 0.;
						gluon.weight = c_entry.weight;
        	                        	gluon.e = g.e();
                	                	gluon.px = g.px();
                        		        gluon.py = g.py();
                                		gluon.pz = g.pz();
						gluon.radlist.clear();
						c_entry.radlist.push_back(gluon);	
					}
				}
				plist.push_back(c_entry); 
			}
    		}
  	}
	std::ofstream of("PythiaSummary.dat");
	of << "# weight-sum\tsigma-gen\tNevent\tNc+Ncbar" << std::endl;
	of << pythia.info.weightSum() << " " << pythia.info.sigmaGen() << " " 
	   << Neve << " " << plist.size() << std::endl;
        std::cout << pythia.info.weightSum() << " " << pythia.info.sigmaGen() << " "
           << Neve << " " << plist.size() << std::endl;

	of.close();
	std::ofstream os("PythiaInput.dat");
	for (auto & p : plist) {
		os << p << std::endl;
		for (auto & g : p.radlist) {
			os << g << std::endl;
		}
	}
	os.close();
	return 0;
}



#include "BeamMC.hh"


BeamMC::BeamMC(double Ee_, double Eh_) {
    Ee = Ee_;
    Eh = Eh_;
}

BeamMC::~BeamMC() {
}

void BeamMC::SetEvent(const podio::Frame* event) {
	mEvent = event;
}

void BeamMC::GetSpecInfo(std::vector<int> &SpecPBG, std::vector<int> &SpecIndex, std::vector<PxPyPzEVector> &SpecVec, std::vector<int> &OtherPBG, std::vector<int> &OtherIndex, std::vector<PxPyPzEVector> &OtherVec)
{
    SpecPBG = spec_pbg;
    SpecIndex = spec_index;
    SpecVec = spec_vec;
    OtherPBG = other_pbg;
    OtherIndex = other_index;
    OtherVec = other_vec;
    
    return;
}

void BeamMC::GetMCinfo(PxPyPzEVector &mc_e,  PxPyPzEVector &mc_p, int &n_pbg) 
{ 
    edm4hep::MCParticle beam_electron;
    edm4hep::MCParticle beam_nucleon;

    spec_pbg.clear();
    spec_index.clear();
    spec_vec.clear();
    other_pbg.clear();
    other_index.clear();
    other_vec.clear();

    const auto& mcparts = static_cast<const edm4hep::MCParticleCollection&>(*(mEvent->get("MCParticles")));

    // cout << endl << "** event" << endl << endl;

	for(const auto& mcp : mcparts)
    {
        if ( mcp.getGeneratorStatus() == 4 )
        {
            if( mcp.getPDG() == ID_ELECTRON )
            {
                beam_electron = mcp;
            } 
            else
            {
                beam_nucleon = mcp;
                n_pbg = mcp.getPDG();
            }

            // cout << mcp << endl;
        }

        if ( mcp.getGeneratorStatus() == 1 )
        {
            if( mcp.getPDG() == ID_NEUTRON || mcp.getPDG() == ID_PROTON )
            {
                int count_id = 0;
                int collect_id[2] = {0,0};
                for (auto it = mcp.parents_begin(), end = mcp.parents_end(); it != end; ++it) 
                {
                    if ( count_id < 2)
                        collect_id[count_id] = it->getObjectID().index;

                    count_id ++;
                }

                if( collect_id[0] == 3 && collect_id[1] == 4 ) // spectators from He3
                {
                    PxPyPzEVector vec;
                    vec.SetPxPyPzE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getEnergy());
                    spec_vec.push_back(vec);
                    spec_pbg.push_back(mcp.getPDG());
                    spec_index.push_back(mcp.getObjectID().index);
                }
                else
                {
                    PxPyPzEVector vec;
                    vec.SetPxPyPzE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getEnergy());
                    other_vec.push_back(vec);
                    other_pbg.push_back(mcp.getPDG());
                    other_index.push_back(mcp.getObjectID().index);
                }
            }
        }
    }

    if ( !beam_electron.isAvailable() || !beam_nucleon.isAvailable() )
        return;

    mc_e.SetPxPyPzE(beam_electron.getMomentum().x, beam_electron.getMomentum().y, beam_electron.getMomentum().z, beam_electron.getEnergy());
    mc_p.SetPxPyPzE(beam_nucleon.getMomentum().x, beam_nucleon.getMomentum().y, beam_nucleon.getMomentum().z, beam_nucleon.getEnergy());

	return;
}
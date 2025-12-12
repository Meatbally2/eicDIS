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

void BeamMC::GetMCinfo(PxPyPzEVector &mc_e,  PxPyPzEVector &mc_p, int &n_pbg) 
{ 
    edm4hep::MCParticle beam_electron;
    edm4hep::MCParticle beam_nucleon;

    auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

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
        }
    }

    if ( !beam_electron.isAvailable() || !beam_nucleon.isAvailable() )
        return;

    mc_e.SetPxPyPzE(beam_electron.getMomentum().x, beam_electron.getMomentum().y, beam_electron.getMomentum().z, beam_electron.getEnergy());
    mc_p.SetPxPyPzE(beam_nucleon.getMomentum().x, beam_nucleon.getMomentum().y, beam_nucleon.getMomentum().z, beam_nucleon.getEnergy());

	return;
}
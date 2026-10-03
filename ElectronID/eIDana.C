// Cut-independent eID inspection ntuple.
// Run from ElectronID/ with eic-shell, for example:
//   root -b -q 'eIDana.C(10,100,1,0,0,-1,500,"../data/sample/your_10x100_sample.root")'

#include "../GlobalUtil/preLoadLib.hh"
#include <TFile.h>
#include <TList.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TString.h>
#include <cmath>

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/AnaManager.cc"

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

constexpr unsigned int kEIDStudySchemaVersion = 7;
constexpr Long64_t kIsolationReferenceCheckLimit = 100;
const std::vector<double> kEIDStudyIsolationRadii = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

std::uint64_t EIDStudyEventKey(const std::string& source, Long64_t entry) {
    // Canonicalize XRootD aliases to the path below EPIC/RECO so changing
    // mirrors does not change event identity. Keep local paths as supplied.
    const std::size_t reco_marker = source.find("/EPIC/RECO/");
    const std::string identity = reco_marker == std::string::npos ? source : source.substr(reco_marker);

    // FNV-1a over the source identity and local entry. The source string and
    // entry are also persisted, so the composite identity remains auditable.
    std::uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : identity) {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    for (int i = 0; i < 8; ++i) {
        hash ^= static_cast<unsigned char>((static_cast<std::uint64_t>(entry) >> (8 * i)) & 0xff);
        hash *= 1099511628211ULL;
    }
    return hash;
}

double EIDStudyDeltaPhi(double a, double b) {
    double dphi = a - b;
    while (dphi > M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    return dphi;
}

std::string EIDStudyOutputName(int Ee, int Eh, int beam_type,
                               int select_region, int sr, int file0) {
    static const char* sample_names[] = {"eHe3", "ep", "piBG", "beamBG", "ep", "ep"};
    std::string name = Form("tmp/%s_%dx%d", sample_names[beam_type], Ee, Eh);
    if (select_region) {
        const double q2_min = (beam_type == AnaManager::PI_BG) ? 0.0 : std::pow(10.0, sr);
        name += Form("_minQ2=%.0f", q2_min);
    }
    name += "_eid_ana";
    if (file0 >= 0)
        name += Form("_f%d", file0);
    return name + ".root";
}

double EIDStudyQ2Minimum(const std::string& source) {
    const std::string marker = "minQ2=";
    const auto position = source.find(marker);
    if (position == std::string::npos)
        return -1.0;
    try {
        return std::stod(source.substr(position + marker.size()));
    } catch (...) {
        return -1.0;
    }
}

std::vector<edm4hep::MCParticle> EIDStudyTruthElectrons(
        const edm4hep::MCParticleCollection& particles) {
    std::vector<int> beam_indices;
    for (const auto& particle : particles) {
        if (particle.getPDG() == 11 && particle.getGeneratorStatus() == 4)
            beam_indices.push_back(particle.getObjectID().index);
    }
    std::sort(beam_indices.begin(), beam_indices.end());

    struct ClosestCandidates {
        int depth = std::numeric_limits<int>::max();
        std::vector<int> indices;
    };
    std::unordered_map<int, ClosestCandidates> per_beam;
    for (int beam_index : beam_indices)
        per_beam.emplace(beam_index, ClosestCandidates{});

    for (const auto& candidate : particles) {
        if (candidate.getPDG() != 11 || candidate.getGeneratorStatus() != 1)
            continue;

        std::queue<std::pair<edm4hep::MCParticle, int>> frontier;
        std::unordered_set<int> visited;
        visited.insert(candidate.getObjectID().index);
        frontier.emplace(candidate, 0);

        while (!frontier.empty()) {
            const auto current = frontier.front().first;
            const int depth = frontier.front().second;
            frontier.pop();
            const int current_index = current.getObjectID().index;

            const auto beam = per_beam.find(current_index);
            if (current.getPDG() == 11 && current.getGeneratorStatus() == 4 && beam != per_beam.end()) {
                if (depth < beam->second.depth) {
                    beam->second.depth = depth;
                    beam->second.indices.clear();
                    beam->second.indices.push_back(candidate.getObjectID().index);
                } else if (depth == beam->second.depth) {
                    beam->second.indices.push_back(candidate.getObjectID().index);
                }
                break;
            }

            for (auto parent_it = current.parents_begin(); parent_it != current.parents_end(); ++parent_it) {
                const auto parent = *parent_it;
                const int parent_index = parent.getObjectID().index;
                if (parent_index >= 0 && visited.insert(parent_index).second)
                    frontier.emplace(parent, depth + 1);
            }
        }
    }

    std::unordered_map<int, edm4hep::MCParticle> by_index;
    for (const auto& particle : particles)
        by_index.emplace(particle.getObjectID().index, particle);

    std::vector<edm4hep::MCParticle> selected;
    for (int beam_index : beam_indices) {
        for (int candidate_index : per_beam[beam_index].indices) {
            const auto found = by_index.find(candidate_index);
            if (found != by_index.end())
                selected.push_back(found->second);
        }
    }
    return selected;
}

void EIDStudyKinematics(int Ee, int Eh, const edm4hep::MCParticle& electron,
                        double& xB, double& Q2, double& W2, double& y, double& nu) {
    const auto p = electron.getMomentum();
    const double incoming_e = std::hypot(Ee, MASS_ELECTRON);
    const double outgoing_e = std::hypot(std::hypot(p.x, p.y, p.z), MASS_ELECTRON);
    const double target_x = Eh * std::sin(CROSSING_ANGLE);
    const double target_z = Eh * std::cos(CROSSING_ANGLE);
    const double target_e = std::hypot(Eh, MASS_PROTON);
    const double q_x = -p.x;
    const double q_y = -p.y;
    const double q_z = -Ee - p.z;
    const double q_e = incoming_e - outgoing_e;
    const double q_dot_target = q_x * target_x + q_z * target_z - q_e * target_e;
    const double incoming_dot_target = -Ee * target_z - incoming_e * target_e;
    Q2 = -(q_e * q_e - q_x * q_x - q_y * q_y - q_z * q_z);
    nu = q_dot_target / MASS_PROTON;
    xB = (nu != 0.0) ? Q2 / (2.0 * MASS_PROTON * nu) : -999.0;
    y = (incoming_dot_target != 0.0) ? q_dot_target / incoming_dot_target : -999.0;
    W2 = MASS_PROTON * MASS_PROTON + 2.0 * MASS_PROTON * nu - Q2;
}

struct EIDStudyEventRow {
    std::string source_file;
    Long64_t source_entry = -1;
    Long64_t source_file_entries = -1;
    ULong64_t event_key = 0;
    double source_q2_min = -1.0;
    int n_candidates = 0;
    int n_truth_electrons = 0;
    int event_header_valid = 0;
    double event_weight = 0.0;
    std::vector<double> event_weights;
    int truth_e_valid = 0;
    double truth_px = -999.0, truth_py = -999.0, truth_pz = -999.0, truth_energy = -999.0;
    double truth_xB = -999.0, truth_Q2 = -999.0, truth_W2 = -999.0, truth_y = -999.0, truth_nu = -999.0;
};

struct EIDStudyCandidateRow {
    ULong64_t event_key = 0;
    int candidate_index = -1;
    int reco_object_index = -1;
    int reco_pdg = 0;
    double charge = 0.0;
    int truth_match_valid = 0;
    int truth_pdg = -999;
    int is_truth_scattered_electron = 0;
    double truth_px = -999.0, truth_py = -999.0, truth_pz = -999.0, truth_energy = -999.0;
    int n_tracks = 0;
    int n_clusters = 0;
    double px = 0.0, py = 0.0, pz = 0.0, energy = 0.0;
    double cluster_energy_sum = 0.0;
    double ecal_detector_cone_energy_r03 = -1.0;
    double hcal_cone_energy_r03 = -1.0;
    int seed_valid = 0;
    double seed_energy = -999.0, seed_eta = -999.0, seed_phi = -999.0;
    std::vector<double> cone_energy_by_radius;
    std::vector<double> cone_fraction_by_radius;
    std::vector<int> cone_association_count_by_radius;
    std::vector<int> track_n_measurements;
    std::vector<int> track_ndf;
    std::vector<float> track_chi2;
    std::vector<int> pid_pdg;
    std::vector<int> pid_type;
    std::vector<float> pid_likelihood;
};

} // namespace

void eIDana(int Ee = 10, int Eh = 100, int beam_type = 1,
            int select_region = 0, int sr = 0, int file0 = -1,
            Long64_t max_events = -1, const char* local_input = "") {
    if (beam_type < 0 || beam_type > 5) {
        std::cerr << "Invalid beam_type " << beam_type << "; expected 0 through 5.\n";
        return;
    }

    AnaManager* ana_manager = new AnaManager("eIDstudy");
    ana_manager->SetBeamEnergy(Ee, Eh);
    ana_manager->Initialize(select_region, sr, file0, beam_type);
    std::vector<std::string> input_names;
    const std::string output_name = EIDStudyOutputName(Ee, Eh, beam_type,
                                                       select_region, sr, file0);
    if (local_input && local_input[0] != '\0') {
        input_names.emplace_back(local_input);
    } else {
        input_names = ana_manager->GetInputNames();
    }
    if (input_names.empty()) {
        std::cerr << "No valid input files resolved; eID study stopped.\n";
        delete ana_manager;
        return;
    }

    gSystem->mkdir("tmp", true);
    TFile output(output_name.c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "Could not create output file " << output_name << "\n";
        delete ana_manager;
        return;
    }

    TNamed schema_version("eIDStudySchemaVersion", std::to_string(kEIDStudySchemaVersion).c_str());
    schema_version.Write();
    TNamed beam_metadata("eIDStudyBeamEnergies", Form("%dx%d GeV", Ee, Eh));
    beam_metadata.Write();
    TNamed beam_type_metadata("eIDStudyBeamType", std::to_string(beam_type).c_str());
    beam_type_metadata.Write();
    TNamed campaign_metadata("eIDStudyCampaign", ana_manager->campaign.c_str());
    campaign_metadata.Write();
    TNamed candidate_definition("eIDStudyCandidateDefinition",
        "Every edm4eic::ReconstructedParticle in the ReconstructedParticles collection; no eID candidate cuts applied.");
    candidate_definition.Write();
    TNamed truth_definition("eIDStudyTruthElectronDefinition",
        "Final-state generator electrons with PDG 11 and status 1, selected by minimum-depth parent ancestry to a status-4 PDG-11 beam electron.");
    truth_definition.Write();
    std::ostringstream isolation_radii_text;
    for (size_t i = 0; i < kEIDStudyIsolationRadii.size(); ++i) {
        if (i) isolation_radii_text << ',';
        isolation_radii_text << kEIDStudyIsolationRadii[i];
    }
    TNamed isolation_radii("eIDStudyIsolationRadii", isolation_radii_text.str().c_str());
    isolation_radii.Write();
    TNamed isolation_definition("eIDStudyIsolationDefinition",
        "Leading associated cluster seeds each cone. DeltaR is measured in eta-phi. Denominator sums every cluster association on every ReconstructedParticle inside DeltaR<R, including self and duplicate associations. Per-radius sums and counts are stored for the listed radius grid.");
    isolation_definition.Write();
    TNamed weight_definition("eIDStudyWeightDefinition",
        "Raw edm4hep::EventHeader.weight and EventHeader.weights are copied when present; zero scalar weights are not interpreted as unit weights or as usable physics normalization. Use sample cross-section and generated-event normalization separately.");
    weight_definition.Write();

    TTree event_tree("Events", "One row per input event, including events with no reconstructed candidates");
    EIDStudyEventRow event_row;
    event_tree.Branch("source_file", &event_row.source_file);
    event_tree.Branch("source_entry", &event_row.source_entry);
    event_tree.Branch("source_file_entries", &event_row.source_file_entries);
    event_tree.Branch("event_key", &event_row.event_key);
    event_tree.Branch("source_q2_min", &event_row.source_q2_min);
    event_tree.Branch("n_candidates", &event_row.n_candidates);
    event_tree.Branch("n_truth_electrons", &event_row.n_truth_electrons);
    event_tree.Branch("event_header_valid", &event_row.event_header_valid);
    event_tree.Branch("event_weight", &event_row.event_weight);
    event_tree.Branch("event_weights", &event_row.event_weights);
    event_tree.Branch("truth_e_valid", &event_row.truth_e_valid);
    event_tree.Branch("truth_px", &event_row.truth_px);
    event_tree.Branch("truth_py", &event_row.truth_py);
    event_tree.Branch("truth_pz", &event_row.truth_pz);
    event_tree.Branch("truth_energy", &event_row.truth_energy);
    event_tree.Branch("truth_xB", &event_row.truth_xB);
    event_tree.Branch("truth_Q2", &event_row.truth_Q2);
    event_tree.Branch("truth_W2", &event_row.truth_W2);
    event_tree.Branch("truth_y", &event_row.truth_y);
    event_tree.Branch("truth_nu", &event_row.truth_nu);

    TTree candidate_tree("Candidates", "All ReconstructedParticles entries before eID selection cuts");
    EIDStudyCandidateRow candidate_row;
    candidate_tree.Branch("event_key", &candidate_row.event_key);
    candidate_tree.Branch("candidate_index", &candidate_row.candidate_index);
    candidate_tree.Branch("reco_object_index", &candidate_row.reco_object_index);
    candidate_tree.Branch("reco_pdg", &candidate_row.reco_pdg);
    candidate_tree.Branch("charge", &candidate_row.charge);
    candidate_tree.Branch("truth_match_valid", &candidate_row.truth_match_valid);
    candidate_tree.Branch("truth_pdg", &candidate_row.truth_pdg);
    candidate_tree.Branch("is_truth_scattered_electron", &candidate_row.is_truth_scattered_electron);
    candidate_tree.Branch("truth_px", &candidate_row.truth_px);
    candidate_tree.Branch("truth_py", &candidate_row.truth_py);
    candidate_tree.Branch("truth_pz", &candidate_row.truth_pz);
    candidate_tree.Branch("truth_energy", &candidate_row.truth_energy);
    candidate_tree.Branch("n_tracks", &candidate_row.n_tracks);
    candidate_tree.Branch("n_clusters", &candidate_row.n_clusters);
    candidate_tree.Branch("track_n_measurements", &candidate_row.track_n_measurements);
    candidate_tree.Branch("track_ndf", &candidate_row.track_ndf);
    candidate_tree.Branch("track_chi2", &candidate_row.track_chi2);
    candidate_tree.Branch("px", &candidate_row.px);
    candidate_tree.Branch("py", &candidate_row.py);
    candidate_tree.Branch("pz", &candidate_row.pz);
    candidate_tree.Branch("energy", &candidate_row.energy);
    candidate_tree.Branch("cluster_energy_sum", &candidate_row.cluster_energy_sum);
    candidate_tree.Branch("ecal_detector_cone_energy_r03", &candidate_row.ecal_detector_cone_energy_r03);
    candidate_tree.Branch("hcal_cone_energy_r03", &candidate_row.hcal_cone_energy_r03);
    candidate_tree.Branch("seed_valid", &candidate_row.seed_valid);
    candidate_tree.Branch("seed_energy", &candidate_row.seed_energy);
    candidate_tree.Branch("seed_eta", &candidate_row.seed_eta);
    candidate_tree.Branch("seed_phi", &candidate_row.seed_phi);
    candidate_tree.Branch("cone_energy_by_radius", &candidate_row.cone_energy_by_radius);
    candidate_tree.Branch("cone_fraction_by_radius", &candidate_row.cone_fraction_by_radius);
    candidate_tree.Branch("cone_association_count_by_radius", &candidate_row.cone_association_count_by_radius);
    candidate_tree.Branch("pid_pdg", &candidate_row.pid_pdg);
    candidate_tree.Branch("pid_type", &candidate_row.pid_type);
    candidate_tree.Branch("pid_likelihood", &candidate_row.pid_likelihood);

    Long64_t total_events = 0;
    Long64_t total_candidates = 0;
    Long64_t total_cluster_associations = 0;
    Long64_t compared_iso_candidates = 0;
    double largest_iso_abs_difference = 0.0;

    for (const std::string& input_name : input_names) {
        podio::ROOTReader reader;
        reader.openFiles({input_name});
        const size_t n_entries = reader.getEntries("events");
        std::cout << "Processing " << input_name << " (" << n_entries << " events)\n";

        for (size_t entry = 0; entry < n_entries; ++entry) {
            if (max_events >= 0 && total_events >= max_events) break;
            auto raw_event = reader.readEntry("events", entry);
            if (!raw_event) {
                std::cerr << "Failed reading entry " << entry << " from " << input_name << "\n";
                continue;
            }
            podio::Frame event(std::move(raw_event));
            const auto& reco_particles = static_cast<const edm4eic::ReconstructedParticleCollection&>(
                *(event.get("ReconstructedParticles")));
            const auto& mc_particles = static_cast<const edm4hep::MCParticleCollection&>(
                *(event.get("MCParticles")));
            const auto& reco_mc_associations = static_cast<const edm4eic::MCRecoParticleAssociationCollection&>(
                *(event.get("ReconstructedParticleAssociations")));
            struct DetectorCluster { double eta, phi, energy; };
            std::vector<DetectorCluster> detector_ecal_clusters, detector_hcal_clusters;
            const auto append_detector_clusters = [](const edm4eic::ClusterCollection& clusters,
                                                      std::vector<DetectorCluster>& destination) {
                for (const auto& cluster : clusters) {
                    const auto& position = cluster.getPosition();
                    destination.push_back({edm4hep::utils::eta(position),
                                           edm4hep::utils::angleAzimuthal(position),
                                           cluster.getEnergy()});
                }
            };
            append_detector_clusters(static_cast<const edm4eic::ClusterCollection&>(*(event.get("EcalBarrelScFiClusters"))), detector_ecal_clusters);
            append_detector_clusters(static_cast<const edm4eic::ClusterCollection&>(*(event.get("EcalEndcapNClusters"))), detector_ecal_clusters);
            append_detector_clusters(static_cast<const edm4eic::ClusterCollection&>(*(event.get("EcalEndcapPClusters"))), detector_ecal_clusters);
            append_detector_clusters(static_cast<const edm4eic::ClusterCollection&>(*(event.get("HcalBarrelClusters"))), detector_hcal_clusters);
            append_detector_clusters(static_cast<const edm4eic::ClusterCollection&>(*(event.get("HcalEndcapNClusters"))), detector_hcal_clusters);
            append_detector_clusters(static_cast<const edm4eic::ClusterCollection&>(*(event.get("LFHCALClusters"))), detector_hcal_clusters);
            const auto truth_electrons = EIDStudyTruthElectrons(mc_particles);

            event_row = EIDStudyEventRow{};
            event_row.source_file = input_name;
            event_row.source_entry = static_cast<Long64_t>(entry);
            event_row.source_file_entries = static_cast<Long64_t>(n_entries);
            event_row.event_key = EIDStudyEventKey(input_name, static_cast<Long64_t>(entry));
            event_row.source_q2_min = EIDStudyQ2Minimum(input_name);
            event_row.n_candidates = static_cast<int>(reco_particles.size());
            event_row.n_truth_electrons = static_cast<int>(truth_electrons.size());
            auto* event_header_base = event.get("EventHeader");
            if (event_header_base) {
                const auto& event_headers = static_cast<const edm4hep::EventHeaderCollection&>(*event_header_base);
                if (!event_headers.empty()) {
                    event_row.event_weight = event_headers[0].getWeight();
                    event_row.event_header_valid = 1;
                    for (double weight : event_headers[0].getWeights())
                        event_row.event_weights.push_back(weight);
                }
            }
            if (!truth_electrons.empty()) {
                const auto& truth_e = truth_electrons[0];
                const auto p = truth_e.getMomentum();
                event_row.truth_e_valid = 1;
                event_row.truth_px = p.x;
                event_row.truth_py = p.y;
                event_row.truth_pz = p.z;
                event_row.truth_energy = truth_e.getEnergy();
                EIDStudyKinematics(Ee, Eh, truth_e, event_row.truth_xB, event_row.truth_Q2,
                                   event_row.truth_W2, event_row.truth_y, event_row.truth_nu);
            }
            event_tree.Fill();

            // Record all associated clusters first so the study tree exactly
            // preserves the denominator loop's per-association counting.
            std::vector<double> all_cluster_eta;
            std::vector<double> all_cluster_phi;
            std::vector<double> all_cluster_energy;
            for (size_t candidate_index = 0; candidate_index < reco_particles.size(); ++candidate_index) {
                const auto& particle = reco_particles[candidate_index];
                for (const auto& cluster : particle.getClusters()) {
                    const auto& position = cluster.getPosition();
                    const double eta = edm4hep::utils::eta(position);
                    const double phi = edm4hep::utils::angleAzimuthal(position);
                    const double energy = cluster.getEnergy();
                    ++total_cluster_associations;
                    all_cluster_eta.push_back(eta);
                    all_cluster_phi.push_back(phi);
                    all_cluster_energy.push_back(energy);
                }
            }

            for (size_t candidate_index = 0; candidate_index < reco_particles.size(); ++candidate_index) {
                const auto& particle = reco_particles[candidate_index];
                const auto momentum = particle.getMomentum();
                candidate_row = EIDStudyCandidateRow{};
                candidate_row.event_key = event_row.event_key;
                candidate_row.candidate_index = static_cast<int>(candidate_index);
                candidate_row.reco_object_index = particle.getObjectID().index;
                candidate_row.reco_pdg = particle.getPDG();
                candidate_row.charge = particle.getCharge();
                candidate_row.n_tracks = static_cast<int>(particle.getTracks().size());
                candidate_row.n_clusters = static_cast<int>(particle.getClusters().size());
                candidate_row.px = momentum.x;
                candidate_row.py = momentum.y;
                candidate_row.pz = momentum.z;
                candidate_row.energy = particle.getEnergy();

                for (const auto& track : particle.getTracks()) {
                    candidate_row.track_n_measurements.push_back(track.measurements_size());
                    candidate_row.track_ndf.push_back(track.getNdf());
                    candidate_row.track_chi2.push_back(track.getChi2());
                }

                edm4hep::MCParticle truth_match;
                for (const auto& association : reco_mc_associations) {
                    if (association.getRec() == particle) {
                        truth_match = association.getSim();
                        if (truth_match.isAvailable()) {
                            candidate_row.truth_match_valid = 1;
                            candidate_row.truth_pdg = truth_match.getPDG();
                            candidate_row.is_truth_scattered_electron =
                                !truth_electrons.empty() &&
                                truth_match.getObjectID().index == truth_electrons[0].getObjectID().index;
                            const auto truth_p = truth_match.getMomentum();
                            candidate_row.truth_px = truth_p.x;
                            candidate_row.truth_py = truth_p.y;
                            candidate_row.truth_pz = truth_p.z;
                            candidate_row.truth_energy = truth_match.getEnergy();
                        }
                        break;
                    }
                }

                const edm4eic::Cluster* leading_cluster = nullptr;
                for (const auto& cluster : particle.getClusters()) {
                    candidate_row.cluster_energy_sum += cluster.getEnergy();
                    if (cluster.getEnergy() > (leading_cluster ? leading_cluster->getEnergy() : 0.0))
                        leading_cluster = &cluster;
                }
                if (leading_cluster) {
                    const auto& position = leading_cluster->getPosition();
                    candidate_row.seed_valid = 1;
                    candidate_row.seed_energy = leading_cluster->getEnergy();
                    candidate_row.seed_eta = edm4hep::utils::eta(position);
                    candidate_row.seed_phi = edm4hep::utils::angleAzimuthal(position);

                    // Legacy eID.C forms E/(E+H) using the candidate-associated
                    // ECal energy and all HCal clusters within DeltaR<0.3. It
                    // also has a gap-modified form using detector-wide ECal
                    // clusters in the same cone. Persist both missing sums.
                    candidate_row.ecal_detector_cone_energy_r03 = 0.0;
                    candidate_row.hcal_cone_energy_r03 = 0.0;
                    const auto sum_detector_cone = [&](const std::vector<DetectorCluster>& clusters) {
                        double sum = 0.0;
                        for (const auto& cluster : clusters) {
                            const double d_eta = cluster.eta - candidate_row.seed_eta;
                            const double d_phi = EIDStudyDeltaPhi(cluster.phi, candidate_row.seed_phi);
                            if (std::hypot(d_eta, d_phi) < 0.3) sum += cluster.energy;
                        }
                        return sum;
                    };
                    candidate_row.ecal_detector_cone_energy_r03 = sum_detector_cone(detector_ecal_clusters);
                    candidate_row.hcal_cone_energy_r03 = sum_detector_cone(detector_hcal_clusters);

                    // Compute each DeltaR once. upper_bound finds the first
                    // radius for which DeltaR < R; prefix sums then populate
                    // all larger cones while retaining repeated associations.
                    std::vector<double> energy_differences(kEIDStudyIsolationRadii.size() + 1, 0.0);
                    std::vector<int> count_differences(kEIDStudyIsolationRadii.size() + 1, 0);
                    for (size_t i = 0; i < all_cluster_energy.size(); ++i) {
                        const double d_eta = all_cluster_eta[i] - candidate_row.seed_eta;
                        const double d_phi = EIDStudyDeltaPhi(all_cluster_phi[i], candidate_row.seed_phi);
                        const double delta_r = std::hypot(d_eta, d_phi);
                        const auto first_radius = std::upper_bound(
                            kEIDStudyIsolationRadii.begin(), kEIDStudyIsolationRadii.end(), delta_r);
                        const size_t radius_index = static_cast<size_t>(
                            std::distance(kEIDStudyIsolationRadii.begin(), first_radius));
                        if (radius_index < kEIDStudyIsolationRadii.size()) {
                            energy_differences[radius_index] += all_cluster_energy[i];
                            ++count_differences[radius_index];
                        }
                    }
                    double cone_energy = 0.0;
                    int cone_associations = 0;
                    size_t r07_index = 0;
                    for (size_t radius_index = 0;
                         radius_index < kEIDStudyIsolationRadii.size(); ++radius_index) {
                        cone_energy += energy_differences[radius_index];
                        cone_associations += count_differences[radius_index];
                        candidate_row.cone_energy_by_radius.push_back(cone_energy);
                        candidate_row.cone_association_count_by_radius.push_back(cone_associations);
                        candidate_row.cone_fraction_by_radius.push_back(
                            cone_energy > 0.0 ? candidate_row.cluster_energy_sum / cone_energy : -1.0);
                        if (kEIDStudyIsolationRadii[radius_index] == 0.7)
                            r07_index = radius_index;
                    }

                    // Validate against a literal copy of the old nested loop on
                    // a bounded subset. Running this O(candidates*clusters)
                    // reference for every row is too expensive at sample scale.
                    if (candidate_row.n_tracks > 0 && candidate_row.n_clusters > 0 &&
                        compared_iso_candidates < kIsolationReferenceCheckLimit) {
                        double current_cone_energy = 0.0;
                        for (const auto& other_particle : reco_particles) {
                            for (const auto& other_cluster : other_particle.getClusters()) {
                                const auto& other_position = other_cluster.getPosition();
                                const double other_eta = edm4hep::utils::eta(other_position);
                                const double other_phi = edm4hep::utils::angleAzimuthal(other_position);
                                const double d_eta = other_eta - candidate_row.seed_eta;
                                const double d_phi = EIDStudyDeltaPhi(other_phi, candidate_row.seed_phi);
                                if (std::hypot(d_eta, d_phi) < 0.7)
                                    current_cone_energy += other_cluster.getEnergy();
                            }
                        }
                        if (current_cone_energy > 0.0) {
                            largest_iso_abs_difference = std::max(largest_iso_abs_difference,
                                std::abs(candidate_row.cone_fraction_by_radius[r07_index] -
                                         candidate_row.cluster_energy_sum / current_cone_energy));
                            ++compared_iso_candidates;
                        }
                    }
                } else {
                    candidate_row.cone_energy_by_radius.assign(kEIDStudyIsolationRadii.size(), 0.0);
                    candidate_row.cone_fraction_by_radius.assign(kEIDStudyIsolationRadii.size(), -1.0);
                    candidate_row.cone_association_count_by_radius.assign(kEIDStudyIsolationRadii.size(), 0);
                }

                for (const auto& pid : particle.getParticleIDs()) {
                    candidate_row.pid_pdg.push_back(pid.getPDG());
                    candidate_row.pid_type.push_back(pid.getType());
                    candidate_row.pid_likelihood.push_back(pid.getLikelihood());
                }
                candidate_tree.Fill();
                ++total_candidates;
            }
            ++total_events;
        }
        if (max_events >= 0 && total_events >= max_events) break;
    }

    output.cd();
    event_tree.Write();
    candidate_tree.Write();
    output.Close();

    std::cout << "Wrote " << output_name << "\n"
              << "Events: " << total_events << ", candidates: " << total_candidates
              << ", cluster associations used for isolation: " << total_cluster_associations << "\n"
              << "R=0.7 isolation benchmark: " << compared_iso_candidates
              << " candidates, maximum absolute difference " << largest_iso_abs_difference << "\n";
    delete ana_manager;
}

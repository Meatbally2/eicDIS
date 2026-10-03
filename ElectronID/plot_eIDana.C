// Readback validation and first-pass QA plots for the compact eID tree.
// Run from ElectronID/ with eic-shell:
//   root -b -q 'plot_eIDana.C("tmp/ep_10x100_eid_ana.root")'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

void plot_eIDana(const char* input_name = "tmp/ep_10x100_eid_ana.root",
                 const char* output_pdf = "") {
    TFile input(input_name, "READ");
    if (input.IsZombie()) {
        std::cerr << "Cannot open " << input_name << '\n';
        return;
    }

    auto* events = dynamic_cast<TTree*>(input.Get("Events"));
    auto* candidates = dynamic_cast<TTree*>(input.Get("Candidates"));
    auto* radii_metadata = dynamic_cast<TNamed*>(input.Get("eIDStudyIsolationRadii"));
    if (!events || !candidates || !radii_metadata) {
        std::cerr << "Missing expected Events/Candidates tree or isolation-radius metadata\n";
        return;
    }

    const char* event_branches[] = {"source_file", "source_entry", "source_file_entries", "source_q2_min",
                                    "event_key", "n_candidates", "event_header_valid", "event_weight",
                                    "event_weights", "truth_Q2", "truth_xB", "truth_y", "truth_W2"};
    const char* candidate_branches[] = {"event_key", "candidate_index", "reco_pdg", "charge",
                                        "truth_match_valid", "truth_pdg", "is_truth_scattered_electron",
                                        "n_tracks", "n_clusters", "px", "py", "pz", "energy",
                                        "cluster_energy_sum", "seed_valid", "seed_eta", "seed_phi",
                                        "cone_energy_by_radius", "cone_fraction_by_radius",
                                        "cone_association_count_by_radius", "pid_pdg", "pid_type",
                                        "pid_likelihood"};
    for (const char* branch : event_branches) {
        if (!events->GetBranch(branch)) {
            std::cerr << "Events is missing branch " << branch << '\n';
            return;
        }
    }
    for (const char* branch : candidate_branches) {
        if (!candidates->GetBranch(branch)) {
            std::cerr << "Candidates is missing branch " << branch << '\n';
            return;
        }
    }

    std::vector<double> radii;
    std::stringstream radius_stream(radii_metadata->GetTitle());
    std::string radius_token;
    while (std::getline(radius_stream, radius_token, ','))
        radii.push_back(std::stod(radius_token));
    size_t r07_index = radii.size();
    for (size_t i = 0; i < radii.size(); ++i)
        if (std::abs(radii[i] - 0.7) < 1e-9) r07_index = i;
    if (r07_index == radii.size()) {
        std::cerr << "The saved radius grid does not include R=0.7\n";
        return;
    }

    std::unordered_map<ULong64_t, Long64_t> expected_candidates;
    TTreeReader event_key_reader(events);
    TTreeReaderValue<ULong64_t> event_key(event_key_reader, "event_key");
    TTreeReaderValue<int> n_candidates_event(event_key_reader, "n_candidates");
    while (event_key_reader.Next())
        expected_candidates.emplace(*event_key, *n_candidates_event);

    Long64_t unmatched_candidates = 0;
    Long64_t wrong_radius_vectors = 0;
    Long64_t truth_electron_candidates = 0;
    std::unordered_map<ULong64_t, Long64_t> observed_candidates;
    TTreeReader candidate_reader(candidates);
    TTreeReaderValue<ULong64_t> candidate_event_key(candidate_reader, "event_key");
    TTreeReaderValue<int> n_tracks(candidate_reader, "n_tracks");
    TTreeReaderValue<int> n_clusters(candidate_reader, "n_clusters");
    TTreeReaderValue<int> truth_electron(candidate_reader, "is_truth_scattered_electron");
    TTreeReaderValue<int> seed_valid(candidate_reader, "seed_valid");
    TTreeReaderValue<double> px(candidate_reader, "px");
    TTreeReaderValue<double> py(candidate_reader, "py");
    TTreeReaderValue<double> pz(candidate_reader, "pz");
    TTreeReaderValue<double> cluster_energy(candidate_reader, "cluster_energy_sum");
    TTreeReaderValue<std::vector<double>> cone_fraction(candidate_reader, "cone_fraction_by_radius");
    TTreeReaderValue<std::vector<double>> cone_energy(candidate_reader, "cone_energy_by_radius");
    TTreeReaderValue<std::vector<int>> cone_count(candidate_reader, "cone_association_count_by_radius");
    TH1D h_candidate_pt("h_candidate_pt", ";candidate p_{T} (GeV);Candidates", 100, 0, 20);
    TH1D h_candidate_eta("h_candidate_eta", ";candidate #eta;Candidates", 100, -5, 5);
    TH1D h_n_tracks("h_n_tracks", ";associated tracks;Candidates", 8, -0.5, 7.5);
    TH1D h_n_clusters("h_n_clusters", ";associated clusters;Candidates", 8, -0.5, 7.5);
    TH1D h_cluster_energy("h_cluster_energy", ";candidate cluster energy sum (GeV);Candidates", 100, 0, 20);
    TH1D h_isolation("h_isolation", ";R=0.7 candidate/cone cluster energy;Candidates", 100, 0, 1.5);
    while (candidate_reader.Next()) {
        const auto event_it = expected_candidates.find(*candidate_event_key);
        if (event_it == expected_candidates.end()) {
            ++unmatched_candidates;
        } else {
            ++observed_candidates[*candidate_event_key];
        }
        if (cone_fraction->size() != radii.size() || cone_energy->size() != radii.size() ||
            cone_count->size() != radii.size()) {
            ++wrong_radius_vectors;
        } else if (*seed_valid && cone_fraction->at(r07_index) >= 0) {
            h_isolation.Fill(cone_fraction->at(r07_index));
        }
        const double pt = std::hypot(*px, *py);
        h_candidate_pt.Fill(pt);
        if (pt > 0) h_candidate_eta.Fill(std::asinh(*pz / pt));
        h_n_tracks.Fill(*n_tracks);
        h_n_clusters.Fill(*n_clusters);
        h_cluster_energy.Fill(*cluster_energy);
        truth_electron_candidates += (*truth_electron != 0);
    }
    Long64_t candidate_count_mismatches = 0;
    for (const auto& item : expected_candidates)
        if (observed_candidates[item.first] != item.second) ++candidate_count_mismatches;

    TTreeReader truth_reader(events);
    TTreeReaderValue<double> truth_q2(truth_reader, "truth_Q2");
    TTreeReaderValue<double> truth_y(truth_reader, "truth_y");
    TH1D h_truth_q2("h_truth_q2", ";truth Q^{2} (GeV^{2});Events", 100, 0, 1000);
    TH1D h_truth_y("h_truth_y", ";truth y;Events", 100, 0, 1);
    while (truth_reader.Next()) {
        if (*truth_q2 > 0) h_truth_q2.Fill(*truth_q2);
        if (*truth_y >= 0 && *truth_y <= 1) h_truth_y.Fill(*truth_y);
    }

    TH1D h_event_weight("h_event_weight", ";EventHeader scalar weight (raw);Events", 100, -10, 10);
    Long64_t event_headers = 0;
    Long64_t nonzero_scalar_weights = 0;
    Long64_t events_with_weight_vectors = 0;
    double sum_event_weight = 0.0;
    double sum_first_aux_weight = 0.0;
    TTreeReader weight_reader(events);
    TTreeReaderValue<int> event_header_valid(weight_reader, "event_header_valid");
    TTreeReaderValue<double> event_weight(weight_reader, "event_weight");
    TTreeReaderValue<std::vector<double>> event_weights(weight_reader, "event_weights");
    while (weight_reader.Next()) {
        if (*event_header_valid) {
            h_event_weight.Fill(*event_weight);
            sum_event_weight += *event_weight;
            ++event_headers;
            nonzero_scalar_weights += (*event_weight != 0.0);
        }
        if (!event_weights->empty()) {
            sum_first_aux_weight += event_weights->front();
            ++events_with_weight_vectors;
        }
    }

    TCanvas canvas("c_eIDana_qa", "eID inspection ntuple QA", 1500, 1200);
    canvas.Divide(3, 3);
    canvas.cd(1); h_truth_q2.Draw();
    canvas.cd(2); h_truth_y.Draw();
    canvas.cd(3); h_candidate_pt.Draw();
    canvas.cd(4); h_candidate_eta.Draw();
    canvas.cd(5); h_n_tracks.Draw();
    canvas.cd(6); h_n_clusters.Draw();
    canvas.cd(7); h_cluster_energy.Draw();
    canvas.cd(8); h_isolation.Draw();
    canvas.cd(9); h_event_weight.Draw();

    std::string pdf = output_pdf;
    if (pdf.empty()) pdf = std::string(input_name) + "_qa.pdf";
    canvas.SaveAs(pdf.c_str());

    std::cout << "Readback validation for " << input_name << '\n'
              << "  Events: " << events->GetEntries() << '\n'
              << "  Candidates: " << candidates->GetEntries() << " (truth scattered-electron matches: "
              << truth_electron_candidates << ")\n"
              << "  Radius grid: " << radii_metadata->GetTitle() << '\n'
              << "  Candidate keys absent from Events: " << unmatched_candidates << '\n'
              << "  Event candidate-count mismatches: " << candidate_count_mismatches << '\n'
              << "  Candidates with malformed radius vectors: " << wrong_radius_vectors << '\n'
              << "  EventHeader rows: " << event_headers << ", nonzero scalar weights: "
              << nonzero_scalar_weights << ", scalar-weight sum: " << sum_event_weight
              << ", events with auxiliary weights: " << events_with_weight_vectors
              << ", sum of first auxiliary weight: " << sum_first_aux_weight << '\n'
              << "  Plots: " << pdf << '\n';
}

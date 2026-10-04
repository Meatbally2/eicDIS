// First-pass weighted study plots for the compact eID inspection ntuple.
// Run from ElectronID/ with eic-shell:
//   root -b -q 'study_eIDana.C("tmp/ep_10x100_eid_ana.root")'
// Multiple inspection files may be passed as a semicolon-separated string.
// The normalization table is for 10x100 ep, matching eff.C:
//   weight = L_target [fb^-1] * sigma(minQ2 region) / N_generated(region).

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TImage.h>
#include <TPad.h>
#include <TLatex.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TSystem.h>
#include <TExec.h>
#include <TProfile.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <limits>
#include <vector>

#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/getBoost.h"

namespace {
constexpr int kRegions = 5; // minQ2=1, 10, 100, 1000, combined
constexpr int kInputRegions = 4;
constexpr int kOverall = 4;
constexpr int kClasses = 4; // scattered e, other e, pi-, other matched truth
constexpr int kMinTrackPoints = 4;
const char* kRegionNames[kRegions] = {"minQ2_1", "minQ2_10", "minQ2_100", "minQ2_1000", "overall"};
const char* kClassNames[kClasses] = {"scattered_e", "other_e", "pi_minus", "other"};
const int kRegionQ2Min[kInputRegions] = {1, 10, 100, 1000};
const double kSigmaFb[kInputRegions] = {5.56009e8, 3.99867e7, 1.34265e6, 6.82136e3};
const int kClassColors[kClasses] = {kRed, kViolet, kBlue, kGray + 2};
std::string gCampaign="unknown";
TImage* gEIDLogo = nullptr;

int RegionIndex(double q2min) {
    for (int i = 0; i < kInputRegions; ++i)
        if (std::abs(q2min - kRegionQ2Min[i]) < 0.5) return i;
    return -1;
}

int TruthClass(int valid, int pdg, int scattered) {
    if (!valid) return 3;
    if (scattered && std::abs(pdg) == 11) return 0;
    if (std::abs(pdg) == 11) return 1;
    if (pdg == -211) return 2;
    return 3;
}

struct EventInfo {
    int region = -1;
    double weight = 0;
};

struct Candidate {
    int index = -1, reco_pdg = 0, truth_pdg = 0, truth_valid = 0, scattered = 0;
    int n_tracks = 0, n_clusters = 0, first_points = 0, seed_valid = 0;
    double charge = 0, px = 0, py = 0, pz = 0, energy = 0;
    double cluster_energy = 0, seed_eta = 0, seed_phi = 0, isolation = -1;
    double ecal_cone = -1, hcal_cone = -1, truth_px = -999, truth_py = -999, truth_pz = -999, truth_energy = -999;
    double le = 0, lpi = 0, lk = 0, lp = 0;
    double P() const { return std::sqrt(px*px+py*py+pz*pz); }
    double Pt() const { return std::hypot(px,py); }
    double Theta() const { return std::atan2(Pt(),pz)*180.0/M_PI; }
    double ClusterTheta() const { return 2.0*std::atan(std::exp(-seed_eta))*180.0/M_PI; }
    double EoEH() const { const double d=cluster_energy+hcal_cone; return d>0 ? cluster_energy/d : 0; }
    double ModEoEH() const { const double d=ecal_cone+hcal_cone; return d>0 ? ecal_cone/d : 0; }
    double PIDe() const { const double d=le+lpi+lk+lp; return d>0 ? le/d : 0; }
    double PIDh() const { const double h=std::max({lpi,lk,lp}); return le+h>0 ? h/(le+h) : 0; }
    double TruthP() const { return std::sqrt(truth_px*truth_px+truth_py*truth_py+truth_pz*truth_pz); }
    double TruthEta() const { return std::asinh(truth_pz/std::hypot(truth_px,truth_py)); }
    bool Base() const { return n_tracks>0 && n_clusters>0 && charge<0 && first_points>=kMinTrackPoints && isolation>=0.9; }
    bool Tight() const { return Base() && EoEH()>0.85; }
    bool Gap() const { return Theta()>158 && Theta()<162 || seed_valid && ClusterTheta()>22 && ClusterTheta()<33; }
};

struct EventCandidates { std::vector<Candidate> rows; };

void FillRatio(TH1D& passed, TH1D& total, double x, bool success, double weight) {
    total.Fill(x, weight);
    if (success) passed.Fill(x, weight);
}

std::unique_ptr<TH1D> Ratio(const TH1D& passed, const TH1D& total, const std::string& name) {
    auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(passed.Clone(name.c_str())));
    result->Divide(&passed, &total, 1.0, 1.0);
    result->SetStats(false);
    result->SetMinimum(0);
    result->SetMaximum(1.5);
    return result;
}

std::unique_ptr<TH1D> MakeCentral68Width(const TH2D& response, const std::string& name) {
    const auto* xaxis = response.GetXaxis();
    auto width = std::make_unique<TH1D>(name.c_str(), ";truth variable;Half-width of central 68% of E/p", 
                                        xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax());
    width->SetDirectory(nullptr);
    width->SetStats(false);
    const double probabilities[] = {0.16, 0.84};
    for (int bin = 1; bin <= xaxis->GetNbins(); ++bin) {
        std::unique_ptr<TH1D> projection(response.ProjectionY(Form("%s_proj_%d", name.c_str(), bin), bin, bin));
        if (!projection || projection->Integral() <= 0) {
            width->SetBinContent(bin, -1.0);
            continue;
        }
        double quantiles[2] = {};
        projection->GetQuantiles(2, quantiles, probabilities);
        const double sigma68 = 0.5 * (quantiles[1] - quantiles[0]);
        width->SetBinContent(bin, sigma68);
        const double effective_entries = projection->GetEffectiveEntries();
        if (effective_entries > 1)
            width->SetBinError(bin, sigma68 / std::sqrt(2.0 * (effective_entries - 1.0)));
    }
    return width;
}

struct Hists {
    std::array<std::array<std::unique_ptr<TH1D>, kClasses>, kRegions> n_track_points;
    std::array<std::array<std::unique_ptr<TH1D>, kClasses>, kRegions> e_over_p;
    std::array<std::array<std::unique_ptr<TH1D>, kClasses>, kRegions> isolation;
    std::array<std::array<std::unique_ptr<TH1D>, kClasses>, kRegions> cluster_energy;
    std::array<std::array<std::unique_ptr<TH1D>, kClasses>, kRegions> e_over_eh, pid_e, pid_h;
    std::array<std::array<std::unique_ptr<TH1D>, 2>, kRegions> gap_eop_b, gap_eop_f, gap_eoeh_b, gap_eoeh_f;
    std::array<std::array<std::unique_ptr<TH1D>, 2>, kRegions> eminus_pz;
    std::array<std::array<std::unique_ptr<TH1D>, 3>, kRegions> selected_mult;
    std::array<std::unique_ptr<TH2D>, kRegions> n_clusters_tracks;
    // stage: all negative tracks (0), selected eID candidate (1), PID present (2).
    // category: electron veto, e, pi, K, p.
    std::array<std::array<std::array<std::array<std::unique_ptr<TH1D>, 2>, 5>, 3>, kRegions> purity_eta, purity_p;
    // categorical PID: general success, purity with PID, selected-candidate purity.
    std::array<std::array<std::array<std::unique_ptr<TH1D>, 2>, 3>, kRegions> pid_categorical;
    std::array<std::array<std::unique_ptr<TH1D>, 2>, kRegions> current_purity;
    std::array<std::array<std::unique_ptr<TH1D>, 2>, kRegions> current_purity_p;
    std::array<std::unique_ptr<TH1D>, kRegions> candidate_pt;
    std::array<std::unique_ptr<TH1D>, kRegions> candidate_eta;
    std::array<std::unique_ptr<TH1D>, kRegions> n_candidates;
    std::array<std::unique_ptr<TH1D>, kRegions> truth_q2;
    std::array<std::unique_ptr<TH2D>, kRegions> pt_theta;
    // Each response map is truth-variable versus reconstructed E/p. Stage 0
    // contains matched electrons with a track and cluster; stage 1 contains
    // the final event-level candidate selected by the baseline eID replay.
    std::array<std::array<std::array<std::unique_ptr<TH2D>, kRegions>, 2>, 4> eop_response;
    std::array<std::array<std::array<std::unique_ptr<TProfile>, kRegions>, 2>, 4> eop_mean;
    std::array<std::array<std::array<std::unique_ptr<TH1D>, kRegions>, 2>, 4> eop_sigma68;
    std::array<std::array<std::unique_ptr<TH2D>, kRegions>, 4> eop_kinematics_electron;
    std::array<std::array<std::unique_ptr<TH2D>, kRegions>, 4> eop_kinematics_negative_hadron;
    std::array<std::array<std::unique_ptr<TH2D>, kRegions>, 4> eoeh_kinematics_electron;
    std::array<std::array<std::unique_ptr<TH2D>, kRegions>, 4> eoeh_kinematics_negative_hadron;

    Hists() {
        for (int r = 0; r < kRegions; ++r) {
            const std::string tag = kRegionNames[r];
            candidate_pt[r] = std::make_unique<TH1D>(("h_pt_" + tag).c_str(), ";candidate p_{T} (GeV);Weighted candidates", 100, 0, 20);
            candidate_eta[r] = std::make_unique<TH1D>(("h_eta_" + tag).c_str(), ";candidate #eta;Weighted candidates", 100, -5, 5);
            n_candidates[r] = std::make_unique<TH1D>(("h_mult_" + tag).c_str(), ";candidates per event;Weighted events", 20, -0.5, 19.5);
            truth_q2[r] = std::make_unique<TH1D>(("h_truth_q2_" + tag).c_str(), ";truth Q^{2} (GeV^{2});Weighted events", 100, 0, 1000);
            pt_theta[r] = std::make_unique<TH2D>(("h_pt_theta_" + tag).c_str(), ";#theta (degrees);p_{T} (GeV)", 180, 0, 180, 50, 0, 50);
            const char* response_names[] = {"theta", "eta", "pt", "energy"};
            const char* response_titles[] = {"truth #theta (degrees)", "truth #eta", "truth p_{T} (GeV)", "truth E_{e} (GeV)"};
            const int response_bins[] = {18, 20, 20, 20};
            const double response_min[] = {0, -5, 0, 0};
            const double response_max[] = {180, 5, 10, 10};
            for (int v = 0; v < 4; ++v) for (int stage = 0; stage < 2; ++stage) {
                const char* stage_name = stage == 0 ? "matched" : "selected";
                eop_response[v][stage][r] = std::make_unique<TH2D>(
                    ("h_eop_vs_truth_" + std::string(response_names[v]) + "_" + stage_name + "_" + tag).c_str(),
                    Form(";%s;E/p", response_titles[v]), response_bins[v], response_min[v], response_max[v], 400, 0, 10);
                eop_response[v][stage][r]->Sumw2();
            }
            for (int v = 0; v < 4; ++v) {
                const char* response_names[] = {"theta", "eta", "pt", "energy"};
                const char* particle_names[] = {"electron", "negative_hadron"};
                auto make_eop_map = [&](const char* particle) {
                    const char* x_title = v == 3
                        ? (std::string(particle) == "electron" ? "truth electron energy (GeV)" : "truth hadron energy (GeV)")
                        : response_titles[v];
                    auto map = std::make_unique<TH2D>(
                        ("h_eop_map_" + std::string(response_names[v]) + "_" + particle + "_" + tag).c_str(),
                        Form(";%s;E/p;Fraction", x_title), response_bins[v], response_min[v], response_max[v], 400, 0, 10);
                    map->Sumw2();
                    return map;
                };
                eop_kinematics_electron[v][r] = make_eop_map(particle_names[0]);
                eop_kinematics_negative_hadron[v][r] = make_eop_map(particle_names[1]);

                auto make_eoeh_map = [&](const char* particle) {
                    const char* x_title = v == 3
                        ? (std::string(particle) == "electron" ? "truth electron energy (GeV)" : "truth hadron energy (GeV)")
                        : response_titles[v];
                    auto map = std::make_unique<TH2D>(
                        ("h_eoeh_map_" + std::string(response_names[v]) + "_" + particle + "_" + tag).c_str(),
                        Form(";%s;E/(E+H);Fraction per truth bin", x_title),
                        response_bins[v], response_min[v], response_max[v], 100, 0, 1);
                    map->Sumw2();
                    return map;
                };
                eoeh_kinematics_electron[v][r] = make_eoeh_map(particle_names[0]);
                eoeh_kinematics_negative_hadron[v][r] = make_eoeh_map(particle_names[1]);
            }
            n_clusters_tracks[r] = std::make_unique<TH2D>(("h_n_clusters_n_tracks_"+tag).c_str(), ";N_{tracks};N_{clusters}",5,-0.5,4.5,5,-0.5,4.5);
            const char* gap_names[] = {"gap_eop_b","gap_eop_f","gap_eoeh_b","gap_eoeh_f"};
            std::array<std::unique_ptr<TH1D>,2>* gaps[] = {&gap_eop_b[r],&gap_eop_f[r],&gap_eoeh_b[r],&gap_eoeh_f[r]};
            for(int g=0;g<4;++g) for(int m=0;m<2;++m)
                (*gaps[g])[m] = std::make_unique<TH1D>((std::string("h_")+gap_names[g]+(m?"_mod_":"_base_")+tag).c_str(),
                    g<2?";E/p;Weighted candidates":";E/(E+H);Weighted candidates",g<2?100:110,0,g<2?2:1.1);
            for(int m=0;m<2;++m) {
                eminus_pz[r][m]=std::make_unique<TH1D>((std::string("h_")+(m?"Cal":"Track")+"EminusPz_"+tag).c_str(),";#Sigma(E-P_{z}) (GeV);Weighted events",200,0,50);
                current_purity[r][m]=std::make_unique<TH1D>((std::string("h_current_purity_")+(m?"total_":"passed_")+tag).c_str(),m?";#eta;Weighted candidates":";#eta;Weighted candidates",20,-5,5);
                current_purity_p[r][m]=std::make_unique<TH1D>((std::string("h_current_purity_p_")+(m?"total_":"passed_")+tag).c_str(),";p (GeV);Purity",150,0,150);
            }
            for(int m=0;m<3;++m) selected_mult[r][m]=std::make_unique<TH1D>(("h_selected_mult_"+std::to_string(m)+"_"+tag).c_str(),";Selected electron candidates;Weighted events",10,-0.5,9.5);
            for(int stage=0;stage<3;++stage) {
                for(int category=0;category<5;++category) for(int pass=0;pass<2;++pass) {
                    const std::string prefix="h_purity_"+std::to_string(stage)+"_"+std::to_string(category)+"_"+std::to_string(pass)+"_"+tag;
                    purity_eta[r][stage][category][pass]=std::make_unique<TH1D>((prefix+"_eta").c_str(),";#eta;Purity",20,-5,5);
                    purity_p[r][stage][category][pass]=std::make_unique<TH1D>((prefix+"_p").c_str(),";p (GeV);Purity",150,0,150);
                }
                for(int pass=0;pass<2;++pass) pid_categorical[r][stage][pass]=std::make_unique<TH1D>(("h_pid_cat_"+std::to_string(stage)+"_"+std::to_string(pass)+"_"+tag).c_str(),";PID category;Purity",5,-0.5,4.5);
            }
            for (int c = 0; c < kClasses; ++c) {
                const std::string ct = tag + "_" + kClassNames[c];
                n_track_points[r][c] = std::make_unique<TH1D>(("h_nTPts_" + ct).c_str(), ";first-track measurements;Weighted candidates", 14, -0.5, 13.5);
                e_over_p[r][c] = std::make_unique<TH1D>(("h_EoP_" + ct).c_str(), ";cluster energy sum / |p|;Weighted candidates", 100, 0, 2);
                isolation[r][c] = std::make_unique<TH1D>(("h_isoE_" + ct).c_str(), ";R=0.7 candidate/cone cluster energy;Weighted candidates", 110, 0, 1.1);
                cluster_energy[r][c] = std::make_unique<TH1D>(("h_clusterE_" + ct).c_str(), ";candidate cluster energy sum (GeV);Weighted candidates", 100, 0, 20);
                e_over_eh[r][c] = std::make_unique<TH1D>(("h_EoEH_"+ct).c_str(),";E/(E+H);Weighted candidates",110,0,1.1);
                pid_e[r][c] = std::make_unique<TH1D>(("h_PIDe_"+ct).c_str(),";L_{e}/#Sigma L;Weighted candidates",100,0,1);
                pid_h[r][c] = std::make_unique<TH1D>(("h_PIDh_"+ct).c_str(),";L_{h}/(L_{h}+L_{e});Weighted candidates",100,0,1);
            }
        }
    }
};

constexpr int kCutFlowBins = 8;
constexpr int kThresholdBins = 11;
constexpr double kIsolationThresholds[kThresholdBins] = {0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
constexpr double kEoEHThresholds[kThresholdBins] = {0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
constexpr int kTrackPointThresholds[kThresholdBins] = {0,1,2,3,4,5,6,7,8,9,10};

struct CutStudyHists {
    std::array<std::unique_ptr<TH1D>, kRegions> signal_event_flow;
    std::array<std::unique_ptr<TH1D>, kRegions> background_candidate_flow;
    // scan index: track measurements, isolation, E/(E+H).
    std::array<std::array<std::unique_ptr<TH1D>, kRegions>, 3> signal_pass, signal_total;
    std::array<std::array<std::unique_ptr<TH1D>, kRegions>, 3> background_pass, background_total;

    CutStudyHists() {
        const char* scan_names[] = {"track_points", "isolation", "eoeh"};
        const char* scan_titles[] = {"N_{track points} threshold", "Isolation fraction threshold", "E/(E+H) threshold"};
        for (int r = 0; r < kRegions; ++r) {
            const std::string tag = kRegionNames[r];
            signal_event_flow[r] = std::make_unique<TH1D>(("h_cutflow_signal_events_"+tag).c_str(), ";Selection stage;Weighted truth-electron events", kCutFlowBins, 0.5, kCutFlowBins+0.5);
            background_candidate_flow[r] = std::make_unique<TH1D>(("h_cutflow_background_candidates_"+tag).c_str(), ";Selection stage;Weighted background candidates", kCutFlowBins, 0.5, kCutFlowBins+0.5);
            signal_event_flow[r]->Sumw2();
            background_candidate_flow[r]->Sumw2();
            for (int s = 0; s < 3; ++s) {
                const int bins = s == 0 ? kThresholdBins : kThresholdBins;
                const double xmin = s == 0 ? -0.5 : (s == 1 ? 0.475 : 0.475);
                const double xmax = s == 0 ? 10.5 : 1.025;
                const std::string base = std::string("h_")+scan_names[s]+"_"+tag;
                signal_pass[s][r] = std::make_unique<TH1D>((base+"_signal_pass").c_str(), Form(";%s;Weighted events",scan_titles[s]), bins, xmin, xmax);
                signal_total[s][r] = std::make_unique<TH1D>((base+"_signal_total").c_str(), Form(";%s;Weighted events",scan_titles[s]), bins, xmin, xmax);
                background_pass[s][r] = std::make_unique<TH1D>((base+"_background_pass").c_str(), Form(";%s;Weighted candidates",scan_titles[s]), bins, xmin, xmax);
                background_total[s][r] = std::make_unique<TH1D>((base+"_background_total").c_str(), Form(";%s;Weighted candidates",scan_titles[s]), bins, xmin, xmax);
                signal_pass[s][r]->Sumw2(); signal_total[s][r]->Sumw2();
                background_pass[s][r]->Sumw2(); background_total[s][r]->Sumw2();
            }
        }
    }
};

std::unique_ptr<TH1D> MakeSurvivalRatio(const TH1D& passed, int denominator_bin, const std::string& name) {
    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D*>(passed.Clone(name.c_str())));
    const double denominator = passed.GetBinContent(denominator_bin);
    const double denominator_var = std::pow(passed.GetBinError(denominator_bin), 2);
    for (int bin = 1; bin <= passed.GetNbinsX(); ++bin) {
        const double numerator = passed.GetBinContent(bin);
        const double numerator_var = std::pow(passed.GetBinError(bin), 2);
        const double value = denominator > 0 ? numerator / denominator : 0;
        // Weighted subset ratio uncertainty, treating the passed and failed
        // contributions as independent sums. Candidate-level correlations
        // within an event are not included.
        const double fail_var = std::max(0.0, denominator_var - numerator_var);
        const double variance = denominator > 0
            ? ((1.0-value)*(1.0-value)*numerator_var + value*value*fail_var) / (denominator*denominator)
            : 0;
        ratio->SetBinContent(bin, value);
        ratio->SetBinError(bin, std::sqrt(std::max(0.0, variance)));
    }
    ratio->SetStats(false);
    ratio->SetMinimum(0);
    ratio->SetMaximum(1.05);
    return ratio;
}

std::unique_ptr<TH1D> MakeEfficiencyRatio(const TH1D& passed, const TH1D& total, const std::string& name) {
    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D*>(passed.Clone(name.c_str())));
    for (int bin = 1; bin <= passed.GetNbinsX(); ++bin) {
        const double denominator = total.GetBinContent(bin);
        const double denominator_var = std::pow(total.GetBinError(bin), 2);
        const double numerator = passed.GetBinContent(bin);
        const double numerator_var = std::pow(passed.GetBinError(bin), 2);
        const double value = denominator > 0 ? numerator / denominator : 0;
        const double fail_var = std::max(0.0, denominator_var - numerator_var);
        const double variance = denominator > 0
            ? ((1.0-value)*(1.0-value)*numerator_var + value*value*fail_var) / (denominator*denominator)
            : 0;
        ratio->SetBinContent(bin, value);
        ratio->SetBinError(bin, std::sqrt(std::max(0.0, variance)));
    }
    ratio->SetStats(false);
    ratio->SetMinimum(0);
    ratio->SetMaximum(1.05);
    return ratio;
}

void AddEIDPlotLabels(TCanvas& canvas, int region, double lumi_fb) {
    DrawManager labels("ep", "10x100 GeV", gCampaign);
    labels.SetEPIC("Work in Progress");
    labels.SetLumi(lumi_fb);
    TCanvas* canvas_ptr = &canvas;
    labels.LableAndCollect(canvas_ptr);

    (void)region;
    canvas.Modified();
    canvas.Update();
}

void SaveStudyCanvas(TCanvas& canvas, const std::string& png_dir,
                     int& page, int region, const std::string& plot_name = "plot") {
    const std::string png_name = Form("%s/%03d_%s_%s.png", png_dir.c_str(), page++,
                                      plot_name.c_str(), kRegionNames[region]);
    canvas.Print(png_name.c_str());
}

void AddGapCanvasLogo(TCanvas& canvas, int region) {
    if (!gEIDLogo) gEIDLogo = TImage::Open("../GlobalUtil/EPIC-logo_black_transparent.png");
    if (!gEIDLogo) return;
    const double width = canvas.GetWw();
    const double height = canvas.GetWh();
    double scale = std::min(width / 1398.0, height / 575.0);
    if (scale > 1.0) scale /= 1.9;
    if (scale < 1.0) scale = std::sqrt(scale);
    const double logo_height = 0.15 * scale;
    const double logo_width = logo_height * double(gEIDLogo->GetWidth()) / gEIDLogo->GetHeight() / (width / height);
    canvas.cd();
    auto* logo_pad = new TPad(Form("gap_logo_pad_%d", region), "",
                              0.19, 0.93 - logo_height, 0.19 + logo_width, 0.93);
    logo_pad->SetFillStyle(0); logo_pad->SetFillColor(0); logo_pad->SetFrameFillStyle(0);
    logo_pad->SetBorderMode(0); logo_pad->SetBorderSize(0); logo_pad->SetFrameBorderMode(0);
    logo_pad->SetLineWidth(0); logo_pad->SetLineColor(0); logo_pad->SetFrameLineWidth(0);
    logo_pad->SetFrameLineColor(0); logo_pad->SetLeftMargin(0); logo_pad->SetBottomMargin(0);
    logo_pad->SetRightMargin(0); logo_pad->SetTopMargin(0);
    logo_pad->Draw(); logo_pad->cd();
    gEIDLogo->SetConstRatio(kTRUE); gEIDLogo->Draw();
    canvas.cd();
}
}

void study_eIDana(const char* input_name = "tmp/ep_10x100_eid_ana.root",
                  double target_lumi_fb = 1.0,
                  const char* output_root = "",
                  const char* output_png_dir = "") {
    // Accept compatible inspection files separated by semicolons, e.g. one
    // file for each generated minQ2 region. The default remains one input file.
    std::vector<std::string> input_files;
    std::stringstream input_stream(input_name);
    std::string input_token;
    while (std::getline(input_stream, input_token, ';'))
        if (!input_token.empty()) input_files.push_back(input_token);
    if (input_files.empty()) { std::cerr << "No input files supplied\n"; return; }

    TChain events("Events"), candidates("Candidates");
    std::string isolation_metadata;
    for (const auto& path : input_files) {
        if (!events.Add(path.c_str()) || !candidates.Add(path.c_str())) {
            std::cerr << "Cannot add inspection file " << path << '\n';
            return;
        }
        TFile check_file(path.c_str(), "READ");
        auto* file_radii = dynamic_cast<TNamed*>(check_file.Get("eIDStudyIsolationRadii"));
        auto* file_schema = dynamic_cast<TNamed*>(check_file.Get("eIDStudySchemaVersion"));
        auto* file_campaign = dynamic_cast<TNamed*>(check_file.Get("eIDStudyCampaign"));
        if (check_file.IsZombie() || !file_radii) {
            std::cerr << "Missing isolation-radius metadata in " << path << '\n';
            return;
        }
        if(!file_schema||std::atoi(file_schema->GetTitle())<7) {
            std::cerr<<"study_eIDana requires inspection schema 7 or newer: "<<path<<'\n';return;
        }
        if(file_campaign) {
            if(gCampaign=="unknown") gCampaign=file_campaign->GetTitle();
            else if(gCampaign!=file_campaign->GetTitle()) gCampaign="mixed";
        }
        if (isolation_metadata.empty()) isolation_metadata = file_radii->GetTitle();
        else if (isolation_metadata != file_radii->GetTitle()) {
            std::cerr << "Incompatible isolation-radius grids in input files\n";
            return;
        }
    }

    std::vector<double> radii;
    std::stringstream radius_stream(isolation_metadata);
    std::string token;
    while (std::getline(radius_stream, token, ',')) radii.push_back(std::stod(token));
    size_t r07 = radii.size();
    for (size_t i = 0; i < radii.size(); ++i)
        if (std::abs(radii[i] - 0.7) < 1e-9) r07 = i;
    if (r07 == radii.size()) {
        std::cerr << "Isolation grid does not contain R=0.7\n";
        return;
    }

    // Count each distinct source file once. These are complete file entry counts,
    // while the Events tree may contain only a processed prefix of each file.
    std::map<std::pair<int, std::string>, Long64_t> generated_by_file;
    std::unordered_map<ULong64_t, EventInfo> event_info;
    TTreeReader er(&events);
    TTreeReaderValue<std::string> source_file(er, "source_file");
    TTreeReaderValue<Long64_t> source_file_entries(er, "source_file_entries");
    TTreeReaderValue<double> source_q2_min(er, "source_q2_min");
    TTreeReaderValue<ULong64_t> event_key(er, "event_key");
    TTreeReaderValue<int> n_candidates(er, "n_candidates");
    TTreeReaderValue<double> truth_q2(er, "truth_Q2");
    TTreeReaderValue<int> truth_valid(er,"truth_e_valid");
    TTreeReaderValue<double> truth_px(er,"truth_px"), truth_py(er,"truth_py"), truth_pz(er,"truth_pz");
    struct EventRow { ULong64_t key; int region; int n; double q2; int truth_valid; double px,py,pz; };
    std::vector<EventRow> event_rows;
    Long64_t skipped_regions = 0;
    while (er.Next()) {
        int region = RegionIndex(*source_q2_min);
        if (region < 0) { ++skipped_regions; continue; }
        auto file_key = std::make_pair(region, *source_file);
        auto inserted = generated_by_file.emplace(file_key, *source_file_entries);
        if (!inserted.second && inserted.first->second != *source_file_entries) {
            std::cerr << "Inconsistent source_file_entries for " << *source_file << '\n';
            return;
        }
        event_rows.push_back({*event_key, region, *n_candidates, *truth_q2,*truth_valid,*truth_px,*truth_py,*truth_pz});
        event_info.emplace(*event_key, EventInfo{region, 0.0});
    }
    std::array<Long64_t, kInputRegions> n_generated{};
    for (const auto& item : generated_by_file) n_generated[item.first.first] += item.second;
    std::array<double, kInputRegions> event_weights{};
    for (int r = 0; r < kInputRegions; ++r)
        if (n_generated[r] > 0)
            event_weights[r] = target_lumi_fb * kSigmaFb[r] / n_generated[r];
    for (auto& row : event_rows) {
        double w = event_weights[row.region];
        event_info[row.key].weight = w;
    }

    // Match eff.C: install ePIC defaults before booking histograms so their
    // axes inherit the user's font, text size, and margins.
    set_ePIC_style();
    Hists h;
    for (const auto& row : event_rows) {
        const double w = event_info[row.key].weight;
        for (int r : {row.region, kOverall}) {
            h.n_candidates[r]->Fill(row.n, w);
            if (row.q2 > 0) h.truth_q2[r]->Fill(row.q2, w);
            if(row.truth_valid) h.pt_theta[r]->Fill(std::atan2(std::hypot(row.px,row.py),row.pz)*180.0/M_PI,std::hypot(row.px,row.py),w);
        }
    }

    TTreeReader cr(&candidates);
    TTreeReaderValue<ULong64_t> candidate_event_key(cr, "event_key");
    TTreeReaderValue<int> truth_match_valid(cr, "truth_match_valid");
    TTreeReaderValue<int> truth_pdg(cr, "truth_pdg");
    TTreeReaderValue<int> truth_scattered(cr, "is_truth_scattered_electron");
    TTreeReaderValue<int> candidate_index(cr,"candidate_index"), reco_pdg(cr,"reco_pdg"), n_clusters(cr,"n_clusters");
    TTreeReaderValue<int> n_tracks(cr, "n_tracks");
    TTreeReaderValue<double> charge(cr, "charge");
    TTreeReaderValue<double> px(cr, "px");
    TTreeReaderValue<double> py(cr, "py");
    TTreeReaderValue<double> pz(cr, "pz");
    TTreeReaderValue<double> energy(cr,"energy"), seed_eta(cr,"seed_eta"), seed_phi(cr,"seed_phi");
    TTreeReaderValue<double> ecal_cone(cr,"ecal_detector_cone_energy_r03"), hcal_cone(cr,"hcal_cone_energy_r03");
    TTreeReaderValue<double> candidate_truth_px(cr,"truth_px"),candidate_truth_py(cr,"truth_py"),candidate_truth_pz(cr,"truth_pz");
    TTreeReaderValue<double> candidate_truth_energy(cr,"truth_energy");
    TTreeReaderValue<std::vector<int>> pid_pdg(cr,"pid_pdg");
    TTreeReaderValue<std::vector<float>> pid_likelihood(cr,"pid_likelihood");
    TTreeReaderValue<double> cluster_energy(cr, "cluster_energy_sum");
    TTreeReaderValue<int> seed_valid(cr, "seed_valid");
    TTreeReaderValue<std::vector<int>> track_measurements(cr, "track_n_measurements");
    TTreeReaderValue<std::vector<double>> cone_fraction(cr, "cone_fraction_by_radius");
    Long64_t missing_keys = 0;
    std::unordered_map<ULong64_t, EventCandidates> candidate_events;
    while (cr.Next()) {
        const auto evt = event_info.find(*candidate_event_key);
        if (evt == event_info.end()) { ++missing_keys; continue; }
        Candidate c;
        c.index=*candidate_index; c.reco_pdg=*reco_pdg; c.truth_pdg=*truth_pdg;
        c.truth_valid=*truth_match_valid; c.scattered=*truth_scattered;
        c.n_tracks=*n_tracks; c.n_clusters=*n_clusters;
        c.first_points=track_measurements->empty()?0:track_measurements->front();
        c.charge=*charge; c.px=*px; c.py=*py; c.pz=*pz; c.energy=*energy;
        c.cluster_energy=*cluster_energy; c.seed_valid=*seed_valid;
        c.seed_eta=*seed_eta; c.seed_phi=*seed_phi;
        c.ecal_cone=*ecal_cone; c.hcal_cone=*hcal_cone;
        c.truth_px=*candidate_truth_px; c.truth_py=*candidate_truth_py; c.truth_pz=*candidate_truth_pz;
        c.truth_energy=*candidate_truth_energy;
        c.isolation=cone_fraction->size()==radii.size()?cone_fraction->at(r07):-1;
        for(size_t i=0;i<pid_pdg->size()&&i<pid_likelihood->size();++i) {
            const double l=pid_likelihood->at(i);
            switch(std::abs(pid_pdg->at(i))) {
                case 11:c.le=std::max(c.le,l);break;
                case 211:c.lpi=std::max(c.lpi,l);break;
                case 321:c.lk=std::max(c.lk,l);break;
                case 2212:c.lp=std::max(c.lp,l);break;
            }
        }
        candidate_events[*candidate_event_key].rows.push_back(c);
        const int region=evt->second.region;
        const double w=evt->second.weight;
        const double eta=c.Pt()>0?std::asinh(c.pz/c.Pt()):-999;
        if (c.truth_valid && c.charge < 0 && c.n_tracks > 0 && c.n_clusters > 0 && c.P() > 0 &&
            c.TruthP() > 0 && c.truth_energy > 0 && c.cluster_energy >= 0) {
            const double truth_pt = std::hypot(c.truth_px, c.truth_py);
            const double truth_theta = std::atan2(truth_pt, c.truth_pz) * 180.0 / M_PI;
            const double truth_eta = c.TruthEta();
            const double eop = c.cluster_energy / c.P();
            const double values[] = {truth_theta, truth_eta, truth_pt, c.truth_energy};
            const int eop_class = c.scattered ? 0 :
                (c.truth_pdg == -211 || c.truth_pdg == -321 || c.truth_pdg == -2212) ? 1 : -1;
            for (int r : {evt->second.region, kOverall}) {
                for (int v = 0; v < 4; ++v) {
                    if (c.scattered)
                        h.eop_response[v][0][r]->Fill(values[v], eop, w);
                    if (eop_class == 0) h.eop_kinematics_electron[v][r]->Fill(values[v], eop, w);
                    else if (eop_class == 1) h.eop_kinematics_negative_hadron[v][r]->Fill(values[v], eop, w);
                }
            }
        }
        if (c.truth_valid && c.charge < 0 && c.n_tracks > 0 && c.n_clusters > 0 && c.P() > 0 &&
            c.TruthP() > 0 && c.truth_energy > 0 && c.cluster_energy >= 0 && c.hcal_cone >= 0 &&
            c.cluster_energy + c.hcal_cone > 0) {
            const double truth_pt = std::hypot(c.truth_px, c.truth_py);
            const double truth_theta = std::atan2(truth_pt, c.truth_pz) * 180.0 / M_PI;
            const double truth_eta = c.TruthEta();
            const double eoeh = c.EoEH();
            const double values[] = {truth_theta, truth_eta, truth_pt, c.truth_energy};
            const int eoeh_class = c.scattered ? 0 :
                (c.truth_pdg == -211 || c.truth_pdg == -321 || c.truth_pdg == -2212) ? 1 : -1;
            for (int r : {evt->second.region, kOverall}) {
                for (int v = 0; v < 4; ++v) {
                    if (eoeh_class == 0) h.eoeh_kinematics_electron[v][r]->Fill(values[v], eoeh, w);
                    else if (eoeh_class == 1) h.eoeh_kinematics_negative_hadron[v][r]->Fill(values[v], eoeh, w);
                }
            }
        }
        for (int r : {region, kOverall}) {
            h.candidate_pt[r]->Fill(c.Pt(), w);
            if (eta > -900) h.candidate_eta[r]->Fill(eta, w);
            // The legacy eID.C feature-comparison loop requires a negative track.
            if(c.n_tracks<=0||c.charge>=0) continue;
            const int cls=TruthClass(c.truth_valid,c.truth_pdg,c.scattered);
            if(cls<0) continue;
            h.n_track_points[r][cls]->Fill(c.first_points,w);
            h.cluster_energy[r][cls]->Fill(c.cluster_energy,w);
            h.e_over_p[r][cls]->Fill(c.n_clusters>0&&c.P()>0?c.cluster_energy/c.P():-1,w);
            h.isolation[r][cls]->Fill(c.n_clusters>0?c.isolation:-1,w);
            h.e_over_eh[r][cls]->Fill(c.n_clusters>0?c.EoEH():-1,w);
            h.pid_e[r][cls]->Fill(c.PIDe(),w);
            h.pid_h[r][cls]->Fill(c.PIDh(),w);
        }
    }

    CutStudyHists cut_study;
    const char* cutflow_labels[kCutFlowBins] = {
        "Truth e", "Matched", "Trk+clu", "q < 0",
        "NTP #geq 4", "Iso #geq 0.9", "E/(E+H)", "Selected e+gap"
    };
    for (int r = 0; r < kRegions; ++r) {
        for (int b = 1; b <= kCutFlowBins; ++b) {
            cut_study.signal_event_flow[r]->GetXaxis()->SetBinLabel(b, cutflow_labels[b-1]);
            cut_study.background_candidate_flow[r]->GetXaxis()->SetBinLabel(b, cutflow_labels[b-1]);
        }
    }
    for (const auto& row : event_rows) {
        const auto event_it = candidate_events.find(row.key);
        const std::vector<Candidate>* rows = event_it == candidate_events.end() ? nullptr : &event_it->second.rows;
        const double w = event_info[row.key].weight;
        bool has_match=false, has_track_cluster=false, has_negative=false, has_points=false, has_isolated=false, has_tight=false;
        bool signal_scan_passed[3][kThresholdBins] = {};
        const Candidate* best_tight=nullptr;
        const Candidate* best_gap=nullptr;
        if (rows) for (const auto& c : *rows) {
            const bool track_cluster = c.n_tracks>0 && c.n_clusters>0;
            const bool negative = c.charge<0;
            const bool points = c.first_points>=kMinTrackPoints;
            const bool isolated = c.isolation>=0.9;
            const bool eoeh = c.EoEH()>0.85;
            if (c.scattered) {
                has_match = true;
                has_track_cluster = has_track_cluster || track_cluster;
                has_negative = has_negative || (track_cluster && negative);
                has_points = has_points || (track_cluster && negative && points);
                has_isolated = has_isolated || (track_cluster && negative && points && isolated);
                has_tight = has_tight || (track_cluster && negative && points && isolated && eoeh);
            }

            if (c.Tight()) {
                if (!best_tight || c.Pt()>best_tight->Pt()) best_tight=&c;
            } else if (c.Base() && c.Gap()) {
                if (!best_gap || c.Pt()>best_gap->Pt()) best_gap=&c;
            }

            const bool background = !c.scattered;
            if (background && track_cluster && negative) {
                for (int r : {row.region,kOverall}) {
                    cut_study.background_candidate_flow[r]->Fill(4,w);
                    if (points) cut_study.background_candidate_flow[r]->Fill(5,w);
                    if (points && isolated) cut_study.background_candidate_flow[r]->Fill(6,w);
                    if (points && isolated && eoeh) cut_study.background_candidate_flow[r]->Fill(7,w);
                }
            }

            if (c.scattered) for (int scan=0; scan<3; ++scan) for (int i=0; i<kThresholdBins; ++i) {
                const double threshold = scan==0 ? kTrackPointThresholds[i] :
                                         scan==1 ? kIsolationThresholds[i] : kEoEHThresholds[i];
                const bool fixed = scan==0 ? (track_cluster && negative && isolated && eoeh) :
                                   scan==1 ? (track_cluster && negative && points && eoeh) :
                                             (track_cluster && negative && points && isolated);
                if (!fixed) continue;
                const bool passes = scan==0 ? c.first_points>=threshold :
                                    scan==1 ? c.isolation>=threshold : c.EoEH()>threshold;
                signal_scan_passed[scan][i] = signal_scan_passed[scan][i] || passes;
            }

            for (int r : {row.region,kOverall}) {
                for (int scan=0; scan<3; ++scan) for (int i=0; i<kThresholdBins; ++i) {
                    const double threshold = scan==0 ? kTrackPointThresholds[i] :
                                             scan==1 ? kIsolationThresholds[i] : kEoEHThresholds[i];
                    const bool fixed = scan==0 ? (track_cluster && negative && isolated && eoeh) :
                                       scan==1 ? (track_cluster && negative && points && eoeh) :
                                                 (track_cluster && negative && points && isolated);
                    if (!fixed) continue;
                    const bool passes = scan==0 ? c.first_points>=threshold :
                                        scan==1 ? c.isolation>=threshold : c.EoEH()>threshold;
                    const double x = threshold;
                    if (background) {
                        cut_study.background_total[scan][r]->Fill(x,w);
                        if (passes) cut_study.background_pass[scan][r]->Fill(x,w);
                    }
                }
            }
        }

        const Candidate* best = best_tight ? best_tight : best_gap;
        if (best && best->scattered && best->n_tracks > 0 && best->n_clusters > 0 &&
            best->P() > 0 && best->TruthP() > 0 && best->truth_energy > 0 && best->cluster_energy >= 0) {
            const double truth_pt = std::hypot(best->truth_px, best->truth_py);
            const double truth_theta = std::atan2(truth_pt, best->truth_pz) * 180.0 / M_PI;
            const double truth_eta = best->TruthEta();
            const double eop = best->cluster_energy / best->P();
            const double values[] = {truth_theta, truth_eta, truth_pt, best->truth_energy};
            for (int r : {row.region, kOverall})
                for (int v = 0; v < 4; ++v)
                    h.eop_response[v][1][r]->Fill(values[v], eop, w);
        }
        for (int r : {row.region,kOverall}) {
            if (row.truth_valid) {
                cut_study.signal_event_flow[r]->Fill(1,w);
                if (has_match) cut_study.signal_event_flow[r]->Fill(2,w);
                if (has_track_cluster) cut_study.signal_event_flow[r]->Fill(3,w);
                if (has_negative) cut_study.signal_event_flow[r]->Fill(4,w);
                if (has_points) cut_study.signal_event_flow[r]->Fill(5,w);
                if (has_isolated) cut_study.signal_event_flow[r]->Fill(6,w);
                if (has_tight) cut_study.signal_event_flow[r]->Fill(7,w);
                if (best && best->scattered) cut_study.signal_event_flow[r]->Fill(8,w);
            }

            for (int scan=0; scan<3; ++scan) for (int i=0; i<kThresholdBins; ++i) {
                const double threshold = scan==0 ? kTrackPointThresholds[i] :
                                         scan==1 ? kIsolationThresholds[i] : kEoEHThresholds[i];
                const double x = threshold;
                if (row.truth_valid) cut_study.signal_total[scan][r]->Fill(x,w);
                if (row.truth_valid && signal_scan_passed[scan][i]) cut_study.signal_pass[scan][r]->Fill(x,w);
            }
        }
    }

    const auto beam_boost=getBoost(10,100,MASS_ELECTRON,MASS_PROTON);
    for(const auto& row:event_rows) {
        const auto it=candidate_events.find(row.key);
        if(it==candidate_events.end()) continue;
        const auto& rows=it->second.rows;
        const double w=event_info[row.key].weight;
        std::vector<const Candidate*> tight,relaxed;
        const Candidate* first_truth=nullptr;
        double track_epz=0,cal_epz=0;
        for(const auto& c:rows) {
            if(c.scattered&&c.n_tracks>0&&c.n_clusters>0) for(int r:{row.region,kOverall}) {
                if(c.ClusterTheta()>22&&c.ClusterTheta()<33) {
                    h.gap_eop_b[r][1]->Fill(c.P()>0?c.ecal_cone/c.P():-1,w);
                    h.gap_eoeh_b[r][1]->Fill(c.ModEoEH(),w);
                }
                if(c.Theta()>158&&c.Theta()<162) {
                    h.gap_eop_f[r][1]->Fill(c.P()>0?c.ecal_cone/c.P():-1,w);
                    h.gap_eoeh_f[r][1]->Fill(c.ModEoEH(),w);
                }
            }
            if(c.scattered&&!first_truth) first_truth=&c;
            if(c.Tight()) tight.push_back(&c);
            else if(c.Base()&&c.Gap()) relaxed.push_back(&c);
            const bool track=c.n_tracks>0,cluster=c.n_clusters>0&&c.seed_valid;
            double t_epz=0,c_epz=0;
            if(track&&c.first_points<4) continue;
            if(track) {
                const auto boosted=beam_boost(PxPyPzEVector(c.px,c.py,c.pz,c.energy));
                t_epz=boosted.E()-boosted.Pz();
            }
            if(cluster) {
                const double pt=c.cluster_energy/std::cosh(c.seed_eta);
                const auto boosted=beam_boost(PxPyPzEVector(pt*std::cos(c.seed_phi),pt*std::sin(c.seed_phi),pt*std::sinh(c.seed_eta),c.cluster_energy));
                c_epz=boosted.E()-boosted.Pz();
            }
            // Match the branches of ElectronID::GetEminusPzSum, including
            // its second addition for track-only particles.
            if(track) {
                if(c.first_points>=4) track_epz+=t_epz;
            } else if(cluster) track_epz+=c_epz;
            if(cluster) cal_epz+=c_epz;
            else if(track) track_epz+=t_epz;
        }
        const auto& selected=tight.empty()?relaxed:tight;
        const Candidate* best=nullptr;
        for(const auto* c:selected) if(!best||c->Pt()>best->Pt()) best=c;
        for(int r:{row.region,kOverall}) {
            if(first_truth) h.n_clusters_tracks[r]->Fill(first_truth->n_tracks,first_truth->n_clusters,w);
            for(const auto& c:rows) {
                if(c.n_tracks<=0||c.charge>=0) continue;
                const bool is_e=c.scattered||c.truth_valid&&std::abs(c.truth_pdg)==11;
                const bool correct[]={!is_e,is_e,c.truth_valid&&std::abs(c.truth_pdg)==211,c.truth_valid&&std::abs(c.truth_pdg)==321,c.truth_valid&&std::abs(c.truth_pdg)==2212};
                const int reco=std::abs(c.reco_pdg);
                const int cat=reco==11?1:reco==211?2:reco==321?3:reco==2212?4:0;
                // Legacy h_pID_suc is filled for every reconstructed particle,
                // independently of track charge or whether a PID exists.
                const int true_cat=is_e?1:0;
                FillRatio(*h.pid_categorical[r][0][0],*h.pid_categorical[r][0][1],true_cat,is_e?reco==11:reco!=11,w);
                if(correct[2]) FillRatio(*h.pid_categorical[r][0][0],*h.pid_categorical[r][0][1],2,reco==211,w);
                if(correct[3]) FillRatio(*h.pid_categorical[r][0][0],*h.pid_categorical[r][0][1],3,reco==321,w);
                if(correct[4]) FillRatio(*h.pid_categorical[r][0][0],*h.pid_categorical[r][0][1],4,reco==2212,w);
                double p=c.P(),eta=c.Pt()>0?std::asinh(c.pz/c.Pt()):-999;
                if(c.truth_valid&&c.n_clusters>0&&c.TruthP()>0) {p=c.TruthP();eta=c.TruthEta();}
                const auto fill_stage=[&](int stage) {
                    const auto fill_cat=[&](int category){
                        if(eta>-900&&std::isfinite(eta)) FillRatio(*h.purity_eta[r][stage][category][0],*h.purity_eta[r][stage][category][1],eta,correct[category],w);
                        FillRatio(*h.purity_p[r][stage][category][0],*h.purity_p[r][stage][category][1],p,correct[category],w);
                    };
                    if(reco!=11) fill_cat(0);
                    if(cat>0) fill_cat(cat);
                };
                fill_stage(0);
                if(reco!=0) fill_stage(2);
                if(reco!=0) {
                    if(reco!=11) FillRatio(*h.pid_categorical[r][2][0],*h.pid_categorical[r][2][1],0,!is_e,w);
                    if(cat>0) FillRatio(*h.pid_categorical[r][2][0],*h.pid_categorical[r][2][1],cat,correct[cat],w);
                }
                if(c.scattered&&c.n_clusters>0&&row.truth_valid) {
                    const double truth_theta=std::atan2(std::hypot(row.px,row.py),row.pz)*180.0/M_PI;
                    if(truth_theta>22&&truth_theta<33) {
                        h.gap_eop_b[r][0]->Fill(c.P()>0?c.cluster_energy/c.P():-1,w);
                        h.gap_eoeh_b[r][0]->Fill(c.EoEH(),w);
                    }
                    if(truth_theta>158&&truth_theta<162) {
                        h.gap_eop_f[r][0]->Fill(c.P()>0?c.cluster_energy/c.P():-1,w);
                        h.gap_eoeh_f[r][0]->Fill(c.EoEH(),w);
                    }
                }
            }
            if(best) {
                h.eminus_pz[r][0]->Fill(track_epz,w);
                h.eminus_pz[r][1]->Fill(cal_epz,w);
                const int n=selected.size();
                h.selected_mult[r][0]->Fill(n,w);
                h.selected_mult[r][best->scattered?1:2]->Fill(n,w);
                const double eta=best->Pt()>0?std::asinh(best->pz/best->Pt()):-999;
                if(eta>-900) FillRatio(*h.current_purity[r][0],*h.current_purity[r][1],eta,best->scattered,w);
                FillRatio(*h.current_purity_p[r][0],*h.current_purity_p[r][1],best->P(),best->scattered,w);
                if(best->reco_pdg!=0) {
                    const int reco=std::abs(best->reco_pdg);
                    const bool is_e=best->scattered||best->truth_valid&&std::abs(best->truth_pdg)==11;
                    const bool correct[]={!is_e,is_e,best->truth_valid&&std::abs(best->truth_pdg)==211,best->truth_valid&&std::abs(best->truth_pdg)==321,best->truth_valid&&std::abs(best->truth_pdg)==2212};
                    const int cat=reco==11?1:reco==211?2:reco==321?3:reco==2212?4:0;
                    double p=best->P(),truth_eta=eta;
                    if(best->truth_valid&&best->TruthP()>0){p=best->TruthP();truth_eta=best->TruthEta();}
                    const auto fill_selected=[&](int category) {
                        if(truth_eta>-900&&std::isfinite(truth_eta)) FillRatio(*h.purity_eta[r][1][category][0],*h.purity_eta[r][1][category][1],truth_eta,correct[category],w);
                        FillRatio(*h.purity_p[r][1][category][0],*h.purity_p[r][1][category][1],p,correct[category],w);
                        FillRatio(*h.pid_categorical[r][1][0],*h.pid_categorical[r][1][1],category,correct[category],w);
                    };
                    if(reco!=11) {
                        fill_selected(0);
                        FillRatio(*h.pid_categorical[r][1][0],*h.pid_categorical[r][1][1],0,!is_e,w);
                    }
                    if(cat>0) fill_selected(cat);
                    if(cat>0) FillRatio(*h.pid_categorical[r][1][0],*h.pid_categorical[r][1][1],cat,correct[cat],w);
                }
            }
        }
    }

    for (int r = 0; r < kRegions; ++r) for (int v = 0; v < 4; ++v) for (int stage = 0; stage < 2; ++stage) {
        const std::string suffix = std::string("_") + kRegionNames[r] + (stage == 0 ? "_matched" : "_selected");
        h.eop_mean[v][stage][r].reset(h.eop_response[v][stage][r]->ProfileX(("p_eop_mean_truth_" + std::to_string(v) + suffix).c_str()));
        h.eop_mean[v][stage][r]->SetDirectory(nullptr);
        h.eop_mean[v][stage][r]->SetStats(false);
        h.eop_sigma68[v][stage][r] = MakeCentral68Width(
            *h.eop_response[v][stage][r], "h_eop_sigma68_truth_" + std::to_string(v) + suffix);
    }

    std::string root_name = output_root;
    std::string png_dir = output_png_dir;
    if (root_name.empty() || png_dir.empty()) {
        std::string stem = input_files.front();
        const auto ext = stem.rfind(".root");
        if (ext != std::string::npos && ext + 5 == stem.size()) stem.erase(ext);
        if (root_name.empty()) root_name = stem + "_study.root";
        if (png_dir.empty()) png_dir = stem + "_study_png";
    }
    gSystem->mkdir(png_dir.c_str(), true);
    TFile output(root_name.c_str(), "RECREATE");
    if (output.IsZombie()) { std::cerr << "Cannot create " << root_name << '\n'; return; }
    for (int r = 0; r < kRegions; ++r) {
        h.candidate_pt[r]->Write(); h.candidate_eta[r]->Write(); h.n_candidates[r]->Write();
        h.truth_q2[r]->Write(); h.pt_theta[r]->Write(); h.n_clusters_tracks[r]->Write();
        for (int v = 0; v < 4; ++v) for (int stage = 0; stage < 2; ++stage) {
            h.eop_response[v][stage][r]->Write();
            h.eop_mean[v][stage][r]->Write();
            h.eop_sigma68[v][stage][r]->Write();
        }
        for (int v = 0; v < 4; ++v) {
            h.eop_kinematics_electron[v][r]->Write();
            h.eop_kinematics_negative_hadron[v][r]->Write();
            h.eoeh_kinematics_electron[v][r]->Write();
            h.eoeh_kinematics_negative_hadron[v][r]->Write();
        }
        for(int m=0;m<2;++m) {
            h.gap_eop_b[r][m]->Write();h.gap_eop_f[r][m]->Write();
            h.gap_eoeh_b[r][m]->Write();h.gap_eoeh_f[r][m]->Write();
            h.eminus_pz[r][m]->Write();h.current_purity[r][m]->Write();h.current_purity_p[r][m]->Write();
        }
        for(int m=0;m<3;++m) {
            h.selected_mult[r][m]->Write();
            for(int pass=0;pass<2;++pass) h.pid_categorical[r][m][pass]->Write();
            for(int cat=0;cat<5;++cat) for(int pass=0;pass<2;++pass) {
                h.purity_eta[r][m][cat][pass]->Write();h.purity_p[r][m][cat][pass]->Write();
            }
        }
        for (int c = 0; c < kClasses; ++c) {
            h.n_track_points[r][c]->Write(); h.e_over_p[r][c]->Write();
            h.isolation[r][c]->Write(); h.cluster_energy[r][c]->Write();
            h.e_over_eh[r][c]->Write(); h.pid_e[r][c]->Write();h.pid_h[r][c]->Write();
        }
        cut_study.signal_event_flow[r]->Write();
        cut_study.background_candidate_flow[r]->Write();
        for (int scan=0; scan<3; ++scan) {
            cut_study.signal_pass[scan][r]->Write();
            cut_study.signal_total[scan][r]->Write();
            cut_study.background_pass[scan][r]->Write();
            cut_study.background_total[scan][r]->Write();
        }
    }
    std::ostringstream normalization;
    normalization << "10x100 ep; target luminosity=" << target_lumi_fb << " fb^-1; ";
    for (int r = 0; r < kInputRegions; ++r)
        normalization << "minQ2=" << kRegionQ2Min[r] << ": sigma=" << kSigmaFb[r]
                      << " fb, N_generated=" << n_generated[r] << ", processed="
                      << std::count_if(event_rows.begin(), event_rows.end(), [r](const EventRow& row){ return row.region == r; })
                      << ", per-event weight=" << event_weights[r] << "; ";
    TNamed norm_note("normalization", normalization.str().c_str());
    norm_note.Write();
    TNamed selection_note("legacy_selection",
        "ElectronID::FindScatteredElectron selection study: negative charge, track and cluster, first-track measurements >=4, R=0.7 isolation >=0.9, E/(E+H)>0.85; if no tight candidates, use candidates in track-theta 158-162 or cluster-theta 22-33 gaps; choose highest reconstructed pT for event-level purity.");
    selection_note.Write();
    TNamed cutflow_definition("cutflow_definition",
        "Signal cutflow counts weighted truth-electron events with at least one matching candidate passing each cumulative requirement; its denominator is Events rows with truth_e_valid. Background cutflow counts weighted non-scattered candidate rows and is normalized to negative-charge candidates with both track and cluster. The final signal stage requires the highest-pT candidate selected by the replayed eID logic to match the truth scattered electron.");
    cutflow_definition.Write();
    TNamed threshold_scan_definition("threshold_scan_definition",
        "Track-measurement, isolation, and E/(E+H) scans vary one active eID threshold while holding the other tight-selection requirements at their current values, including the four-measurement track requirement where it is not being scanned. Signal efficiency is event-level relative to all truth-electron Events rows; background acceptance is candidate-level relative to the sample passing the other fixed requirements. It is not the current event-level fake rate. E/p and PID veto are not active requirements in ElectronID::FindScatteredElectron and are therefore not treated as current cuts. Weighted ratio errors use independent passed/failed sumw2 and do not include correlations among candidates in one event.");
    threshold_scan_definition.Write();
    TNamed eop_resolution_definition("eop_resolution_definition",
        "E/p response maps use associated ECal cluster-energy sum divided by reconstructed track momentum for truth-matched scattered electrons with both track and cluster. Truth theta, eta, pT, and truth scattered-electron energy define the bins. The matched stage is before baseline eID cuts; the selected stage is the final highest-pT candidate from the event-level baseline replay, including its gap fallback, only when truth matched. The response is the weighted mean E/p; resolution is half the weighted central 68 percent interval of E/p per truth bin. Region weights are applied before computing both summaries. Resolution errors are the approximate Gaussian width uncertainty sigma68/sqrt(2*(N_eff-1)).");
    eop_resolution_definition.Write();
    TNamed eop_map_definition("eop_map_definition",
        "The 2D E/p maps show truth-matched scattered electrons and negative reconstructed tracks matched to pi-minus, K-minus, or antiproton, before the baseline eID cuts and requiring both track and cluster. E/p is cluster-associated ECal energy divided by reconstructed track momentum. Each displayed map is normalized to unit weighted integral; the ROOT map histograms retain weighted candidate yields. The hadron energy axis uses matched truth-particle energy.");
    eop_map_definition.Write();
    TNamed comparison_note("legacy_plot_comparison",
        "The pion overlay uses truth PDG -211 and unmatched candidates appear in Others. In eID.C the comparison class is abs(mc_pdg), making its -211 pion branch unreachable. The first truth-matched reconstructed particle for N_tracks vs N_clusters is taken in Candidates tree order, which may differ from association order.");
    comparison_note.Write();
    output.Close();

    TCanvas canvas("c_eID_study", "eID weighted study", 1000, 600);
    int page = 1;
    // Render the weighted combination only. Per-region histograms remain in
    // the ROOT output for targeted diagnostics; the combined 2D E/p maps show
    // the kinematic dependence without multiplying the plot count by Q2 region.
    for (int r = kOverall; r < kRegions; ++r) {
        if (r < kInputRegions && n_generated[r] == 0) {
            canvas.Clear();
            TLatex missing;
            missing.SetNDC();
            missing.SetTextAlign(22);
            missing.SetTextFont(42);
            missing.SetTextSize(0.045);
            missing.DrawLatex(0.5, 0.55, Form("%s: no inspection input supplied", kRegionNames[r]));
            missing.DrawLatex(0.5, 0.45, "This Q^{2} region is not included in the combined plots");
            AddEIDPlotLabels(canvas, r, target_lumi_fb);
            SaveStudyCanvas(canvas, png_dir, page, r, "missing_region_input");
            continue;
        }
        for (int feature = 0; feature < 7; ++feature) {
            canvas.Clear();
            canvas.SetLogy(feature>0&&feature<6);
            TH1D* hs[kClasses];
            for (int c = 0; c < kClasses; ++c) {
                hs[c] = feature == 0 ? h.n_track_points[r][c].get() :
                        feature == 1 ? h.e_over_p[r][c].get() :
                        feature == 2 ? h.isolation[r][c].get() :
                        feature == 3 ? h.e_over_eh[r][c].get() :
                        feature == 4 ? h.pid_e[r][c].get() :
                        feature == 5 ? h.pid_h[r][c].get() : h.cluster_energy[r][c].get();
                hs[c]->SetLineColor(kClassColors[c]); hs[c]->SetLineWidth(2);
                hs[c]->SetStats(false);
            }
            hs[0]->SetFillColor(kRed);
            hs[0]->SetFillStyle(3003);
            const char* xlabel[]={"N_{track points}","E/p","Isolation fraction (R=0.7)","E/(E+H)","L_{e}/#Sigma L","L_{h}/(L_{h}+L_{e})","E_{cluster} (GeV)"};
            const char* feature_names[]={"track_points_by_class","eop_by_class","isolation_by_class","e_over_eh_by_class","pide_by_class","pidh_by_class","cluster_energy_by_class"};
            hs[3]->SetTitle(Form(";%s;Expected candidates",xlabel[feature]));
            double ymax=0;
            for(auto* hist:hs) ymax=std::max(ymax,hist->GetMaximum());
            hs[3]->SetMaximum(ymax>0?1.2*ymax:1.0);
            hs[3]->SetMinimum(canvas.GetLogy()?0.1:0.0);
            hs[3]->Draw("hist");
            hs[2]->Draw("hist same");
            hs[1]->Draw("hist same");
            hs[0]->Draw("hist same");
            const double legend_x = feature == 1 ? 0.70 : 0.42;
            TLegend legend(legend_x, 0.60, legend_x + 0.25, 0.88);
            legend.SetBorderSize(0);
            legend.SetFillStyle(0);
            legend.SetTextFont(42);
            legend.AddEntry(hs[0], "Electrons", "l");
            legend.AddEntry(hs[1], "Other e's", "l");
            legend.AddEntry(hs[2], "Pions", "l");
            legend.AddEntry(hs[3], "Others", "l");
            legend.Draw();
            const std::vector<double> cuts=feature==0?std::vector<double>{3.5}:
                feature==1?std::vector<double>{0.8,1.2}:
                feature==2?std::vector<double>{0.9}:
                feature==3?std::vector<double>{0.85}:
                feature==5?std::vector<double>{0.62}:std::vector<double>{};
            for(double cut:cuts) {
                TLine line(cut,canvas.GetLogy()?0.1:0,cut,hs[3]->GetMaximum());
                line.SetLineStyle(7);line.DrawClone("same");
            }
            AddEIDPlotLabels(canvas, r, target_lumi_fb);
            SaveStudyCanvas(canvas, png_dir, page, r, feature_names[feature]);
        }
        {
            const char* variable_titles[2][4] = {
                {"truth #theta (degrees)", "truth #eta", "truth p_{T} (GeV)", "truth electron energy (GeV)"},
                {"truth #theta (degrees)", "truth #eta", "truth p_{T} (GeV)", "truth hadron energy (GeV)"}
            };
            const char* sample_labels[] = {"#bf{e^{-}}", "#bf{h^{-}}"};
            const char* variable_names[] = {"theta", "eta", "pt", "energy"};
            for (int v = 0; v < 4; ++v) for (int sample = 0; sample < 2; ++sample) {
                canvas.Clear();
                canvas.SetCanvasSize(1400, 600);
                canvas.SetLogz(true);
                canvas.SetLogx(false); canvas.SetLogy(false);
                canvas.SetLeftMargin(0.13); canvas.SetBottomMargin(0.13); canvas.SetRightMargin(0.16);
                auto* source = sample == 0 ? h.eop_kinematics_electron[v][r].get()
                                           : h.eop_kinematics_negative_hadron[v][r].get();
                std::unique_ptr<TH2D> display(static_cast<TH2D*>(source->Clone(
                    Form("display_eop_map_%d_%d_%d", r, v, sample))));
                display->SetDirectory(nullptr);
                display->SetStats(false);
                const double integral = display->Integral(0, display->GetNbinsX() + 1,
                                                         0, display->GetNbinsY() + 1);
                if (integral > 0) display->Scale(1.0 / integral);
                display->SetTitle(Form("%s;%s;E/p;Fraction",
                                       sample_labels[sample], variable_titles[sample][v]));
                display->GetYaxis()->SetRangeUser(0, 2);
                display->GetXaxis()->CenterTitle(); display->GetYaxis()->CenterTitle(); display->GetZaxis()->CenterTitle();
                display->GetXaxis()->SetTitleOffset(1.2); display->GetYaxis()->SetTitleOffset(1.2);
                display->GetZaxis()->SetTitleOffset(1.1);
                display->SetMinimum(1e-5);
                display->GetListOfFunctions()->Add(new TExec(
                    Form("eop_palette_%d_%d_%d", r, v, sample), "gStyle->SetPalette(kLightTemperature);"));
                display->Draw("colz");
                // Match eff.C: DrawManager adds its ePIC block to the same
                // canvas as the 2D histogram, after drawing the map.
                AddEIDPlotLabels(canvas, r, target_lumi_fb);
                TLatex sample_label;
                sample_label.SetNDC(); sample_label.SetTextFont(42); sample_label.SetTextSize(0.055); sample_label.SetTextAlign(12);
                sample_label.DrawLatex(0.72, 0.855, sample_labels[sample]);
                const std::string variable_name = v == 3
                    ? (sample == 0 ? "electron_energy" : "hadron_energy")
                    : variable_names[v];
                SaveStudyCanvas(canvas, png_dir, page, r,
                    std::string("eop_") + (sample == 0 ? "electron_vs_" : "negative_hadrons_vs_") + variable_name);
                canvas.SetCanvasSize(1000, 600);
                canvas.SetLogz(false);
            }
        }
        {
            const char* variable_titles[2][4] = {
                {"truth #theta (degrees)", "truth #eta", "truth p_{T} (GeV)", "truth electron energy (GeV)"},
                {"truth #theta (degrees)", "truth #eta", "truth p_{T} (GeV)", "truth hadron energy (GeV)"}
            };
            const char* sample_labels[] = {"#bf{e^{-}}", "#bf{h^{-}}"};
            const char* variable_names[] = {"theta", "eta", "pt", "energy"};
            for (int v = 0; v < 4; ++v) for (int sample = 0; sample < 2; ++sample) {
                canvas.Clear();
                canvas.SetCanvasSize(1400, 600);
                canvas.SetLogz(true);
                canvas.SetLogx(false); canvas.SetLogy(false);
                canvas.SetLeftMargin(0.13); canvas.SetBottomMargin(0.13); canvas.SetRightMargin(0.16);
                canvas.SetTopMargin(0.37); // reserve the ePIC header above the E/(E+H) peak near one
                auto* source = sample == 0 ? h.eoeh_kinematics_electron[v][r].get()
                                           : h.eoeh_kinematics_negative_hadron[v][r].get();
                std::unique_ptr<TH2D> display(static_cast<TH2D*>(source->Clone(
                    Form("display_eoeh_map_%d_%d_%d", r, v, sample))));
                display->SetDirectory(nullptr);
                display->SetStats(false);
                // Normalize each truth-kinematic bin independently. Each vertical
                // slice then shows the E/(E+H) distribution conditional on x.
                for (int xbin = 1; xbin <= display->GetNbinsX(); ++xbin) {
                    const double integral = display->Integral(xbin, xbin, 1, display->GetNbinsY());
                    if (integral <= 0) continue;
                    for (int ybin = 1; ybin <= display->GetNbinsY(); ++ybin) {
                        display->SetBinContent(xbin, ybin, display->GetBinContent(xbin, ybin) / integral);
                        display->SetBinError(xbin, ybin, display->GetBinError(xbin, ybin) / integral);
                    }
                }
                display->SetTitle(Form("%s;%s;E/(E+H);Fraction per truth bin",
                                       sample_labels[sample], variable_titles[sample][v]));
                display->GetYaxis()->SetRangeUser(0, 1);
                display->GetXaxis()->CenterTitle(); display->GetYaxis()->CenterTitle(); display->GetZaxis()->CenterTitle();
                display->GetXaxis()->SetTitleOffset(1.2); display->GetYaxis()->SetTitleOffset(1.2);
                display->GetZaxis()->SetTitleOffset(1.1);
                display->SetMinimum(1e-4);
                display->GetListOfFunctions()->Add(new TExec(
                    Form("eoeh_palette_%d_%d_%d", r, v, sample), "gStyle->SetPalette(kLightTemperature);"));
                display->Draw("colz");
                TLine cut_line(display->GetXaxis()->GetXmin(), 0.85,
                               display->GetXaxis()->GetXmax(), 0.85);
                cut_line.SetLineColor(kRed + 1); cut_line.SetLineStyle(7); cut_line.SetLineWidth(2);
                cut_line.DrawClone("same");
                AddEIDPlotLabels(canvas, r, target_lumi_fb);
                TLatex sample_label;
                sample_label.SetNDC(); sample_label.SetTextFont(42); sample_label.SetTextSize(0.055); sample_label.SetTextAlign(12);
                sample_label.DrawLatex(0.72, 0.855, sample_labels[sample]);
                const std::string variable_name = v == 3
                    ? (sample == 0 ? "electron_energy" : "hadron_energy")
                    : variable_names[v];
                SaveStudyCanvas(canvas, png_dir, page, r,
                    std::string("eoeh_") + (sample == 0 ? "electron_vs_" : "negative_hadrons_vs_") + variable_name);
                canvas.SetCanvasSize(1000, 600);
                canvas.SetLogz(false);
            }
            canvas.SetTopMargin(0.1);
        }
        canvas.SetLogy(false);
        const char* gap_titles[]={"Backward gap E/p","Forward gap E/p","Backward gap E/(E+H)","Forward gap E/(E+H)"};
        std::array<std::unique_ptr<TH1D>,2>* gaps[]={&h.gap_eop_b[r],&h.gap_eop_f[r],&h.gap_eoeh_b[r],&h.gap_eoeh_f[r]};
        canvas.Clear();canvas.SetCanvasSize(1400,800);canvas.cd();
        AddEIDPlotLabels(canvas, r, target_lumi_fb);
        const double pad_x1[]={0.08,0.55,0.08,0.55};
        const double pad_x2[]={0.48,0.95,0.48,0.95};
        const double pad_y1[]={0.33,0.33,0.02,0.02};
        const double pad_y2[]={0.62,0.62,0.31,0.31};
        for(int g=0;g<4;++g) {
            auto* pad=new TPad(Form("gap_pad_%d_%d",r,g),"",pad_x1[g],pad_y1[g],pad_x2[g],pad_y2[g]);
            pad->SetLeftMargin(0.16);pad->SetRightMargin(0.05);pad->SetBottomMargin(0.16);pad->SetTopMargin(0.08);
            pad->Draw();pad->cd();pad->SetLogy();
            auto* base=(*gaps[g])[0].get();auto* modified=(*gaps[g])[1].get();
            base->SetTitle(gap_titles[g]);base->SetLineColor(kGray+2);modified->SetLineColor(kRed);
            base->SetMaximum(std::max(base->GetMaximum(),modified->GetMaximum())*1.4+1);
            base->SetMinimum(0.1);base->Draw("hist");modified->Draw("hist same");
            canvas.cd();
        }
        canvas.cd();
        TLegend gap_legend(0.69,0.66,0.93,0.76);gap_legend.SetBorderSize(0);gap_legend.SetTextSize(0.025);
        gap_legend.AddEntry(h.gap_eop_b[r][0].get(),"Baseline E","l");
        gap_legend.AddEntry(h.gap_eop_b[r][1].get(),"Detector cone E","l");gap_legend.Draw();
        AddGapCanvasLogo(canvas, r);
        SaveStudyCanvas(canvas, png_dir, page, r, "gap_response_four_panels");canvas.Clear();canvas.SetCanvasSize(1000,600);canvas.SetLogy(false);
        {
            auto* track=h.eminus_pz[r][0].get();auto* cal=h.eminus_pz[r][1].get();
            track->SetLineColor(kBlue);cal->SetLineColor(kGray+2);track->SetFillColor(kBlue);track->SetFillStyle(3003);
            cal->SetMaximum(1.4*std::max(track->GetMaximum(),cal->GetMaximum())+1);cal->Draw("hist");track->Draw("hist same");
            TLegend legend(0.65,0.42,0.88,0.59);legend.SetBorderSize(0);legend.AddEntry(track,"Using E_{Track}","l");legend.AddEntry(cal,"Using E_{Cluster}","l");legend.Draw();
            TLine cut(20,0,20,cal->GetMaximum());cut.SetLineStyle(7);cut.DrawClone("same");
            AddEIDPlotLabels(canvas,r,target_lumi_fb);SaveStudyCanvas(canvas, png_dir, page, r, "eminus_pz_track_vs_cluster");
        }
        canvas.Clear();canvas.SetLogy(true);
        {
            const int colors[]={kGray+2,kBlue,kOrange+7};
            TLegend legend(0.60,0.39,0.88,0.60);legend.SetBorderSize(0);
            const char* labels[]={"All candidates","Scat. e has highest p_{T}","Others have highest p_{T}"};
            for(int m=0;m<3;++m) {
                auto* hist=h.selected_mult[r][m].get();hist->SetLineColor(colors[m]);hist->SetLineWidth(2);
                if(m==0){hist->SetMinimum(0.1);hist->SetMaximum(1.4*hist->GetMaximum()+1);hist->Draw("hist");}
                else hist->Draw("hist same");
                legend.AddEntry(hist,labels[m],"l");
            }
            legend.Draw();AddEIDPlotLabels(canvas,r,target_lumi_fb);SaveStudyCanvas(canvas, png_dir, page, r, "selected_candidate_multiplicity");
        }
        canvas.Clear();canvas.SetLogy(false);
        {
            auto hist=std::unique_ptr<TH2D>(static_cast<TH2D*>(h.n_clusters_tracks[r]->Clone(Form("display_clusters_%d",r))));
            if(hist->Integral()>0) hist->Scale(1.0/hist->Integral());
            hist->SetStats(false);hist->Draw("colz text");AddEIDPlotLabels(canvas,r,target_lumi_fb);SaveStudyCanvas(canvas, png_dir, page, r, "cluster_vs_track_multiplicity");
        }
        canvas.Clear();
        {
            const char* names[]={"Not e","e","#pi","K","p"};
            const int stage_order[]={1,2,0};
            const int colors[]={kP8Green,kP8Blue,kP8Red};
            const int markers[]={29,21,20};
            const char* stages[]={"Purity for e candidates","Purity if pID exists","Purity in general"};
            TLegend legend(0.60,0.39,0.88,0.60);legend.SetBorderSize(0);
            std::vector<std::unique_ptr<TH1D>> ratios;
            for(int series=0;series<3;++series) {
                const int stage=stage_order[series];
                ratios.push_back(Ratio(*h.pid_categorical[r][stage][0],*h.pid_categorical[r][stage][1],Form("pid_ratio_%d_%d",r,stage)));
                auto* hist=ratios.back().get();hist->SetLineColor(colors[series]);hist->SetMarkerColor(colors[series]);hist->SetMarkerStyle(markers[series]);
                for(int i=0;i<5;++i) hist->GetXaxis()->SetBinLabel(i+1,names[i]);
                if(series==0) {hist->SetFillColor(kP8Green);hist->SetFillStyle(3003);hist->SetMarkerSize(0.0);hist->Draw("hist");}
                else hist->Draw("p e1 same");
                legend.AddEntry(hist,stages[series],series==0?"l":"lp");
            }
            legend.Draw();AddEIDPlotLabels(canvas,r,target_lumi_fb);SaveStudyCanvas(canvas, png_dir, page, r, "pid_categorical_purity");
        }
        const char* stage_names[]={"Negative-track PID purity","Selected-candidate PID purity","All-particle PID purity"};
        const char* labels[]={"electron veto","electron","pion","kaon","proton"};
        const int colors[]={kP10Red,kP10Blue,kP10Brown,kP10Green,kP10Ash};
        const int markers[]={29,20,21,22,23};
        for(int stage=0;stage<3;++stage) for(int axis=0;axis<2;++axis) {
            canvas.Clear();canvas.SetLogy(false);
            TLegend legend(0.62,0.38,0.88,0.63);legend.SetBorderSize(0);
            std::vector<std::unique_ptr<TH1D>> ratios;
            auto& current_pair=axis?h.current_purity_p[r]:h.current_purity[r];
            ratios.push_back(Ratio(*current_pair[0],*current_pair[1],Form("current_ratio_%d_%d_%d",r,stage,axis)));
            auto* baseline=ratios.back().get();baseline->SetLineColor(kP10Yellow);baseline->SetFillColor(kP10Yellow);
            baseline->SetLineStyle(1);baseline->SetLineWidth(2);baseline->SetFillStyle(3003);baseline->SetMarkerSize(0.0);
            baseline->SetTitle(Form(";%s;Purity",axis?"p [GeV/c]":"#eta"));baseline->Draw("hist");
            legend.AddEntry(baseline,"baseline","l");
            for(int cat=0;cat<5;++cat) {
                auto& pair=axis?h.purity_p[r][stage][cat]:h.purity_eta[r][stage][cat];
                ratios.push_back(Ratio(*pair[0],*pair[1],Form("purity_ratio_%d_%d_%d_%d",r,stage,axis,cat)));
                auto* hist=ratios.back().get();hist->SetLineColor(colors[cat]);hist->SetMarkerColor(colors[cat]);hist->SetMarkerStyle(markers[cat]);hist->SetMarkerSize(1.0);hist->SetTitle(Form(";%s;Purity",axis?"p [GeV/c]":"#eta"));
                hist->Draw("p e1 same");legend.AddEntry(hist,labels[cat],"lp");
            }
            legend.Draw();AddEIDPlotLabels(canvas,r,target_lumi_fb);
            TLatex heading;heading.SetNDC();heading.SetTextFont(42);heading.SetTextSize(0.027);
            heading.DrawLatex(0.16,0.62,stage_names[stage]);
            const char* stage_tags[]={"negative_track_pid","selected_candidate_pid","all_particle_pid"};
            SaveStudyCanvas(canvas, png_dir, page, r,
                std::string("purity_") + stage_tags[stage] + (axis ? "_vs_pt" : "_vs_eta"));
        }
        h.candidate_pt[r]->SetTitle(";p_{T} (GeV);Expected candidates");
        h.candidate_eta[r]->SetTitle(";#eta;Expected candidates");
        h.n_candidates[r]->SetTitle(";Candidates per event;Expected events");
        h.truth_q2[r]->SetTitle(";Truth Q^{2} (GeV^{2});Expected events");
        h.pt_theta[r]->SetTitle(";#theta (degrees);p_{T} (GeV)");
        h.candidate_pt[r]->SetStats(false); h.candidate_eta[r]->SetStats(false);
        h.n_candidates[r]->SetStats(false); h.truth_q2[r]->SetStats(false);
        h.pt_theta[r]->SetStats(false);
        const char* one_d_names[]={"candidate_pt","candidate_eta","candidate_multiplicity","truth_q2"};
        for (int i = 0; i < 4; ++i) {
            canvas.Clear();
            TH1D* one_d = i == 0 ? h.candidate_pt[r].get() : i == 1 ? h.candidate_eta[r].get() : i == 2 ? h.n_candidates[r].get() : h.truth_q2[r].get();
            one_d->SetLineColor(kGray + 2);
            one_d->SetLineWidth(2);
            one_d->Draw("hist");
            AddEIDPlotLabels(canvas, r, target_lumi_fb);
            SaveStudyCanvas(canvas, png_dir, page, r, one_d_names[i]);
        }
        canvas.Clear();
        h.pt_theta[r]->Draw("colz");
        AddEIDPlotLabels(canvas, r, target_lumi_fb);
        SaveStudyCanvas(canvas, png_dir, page, r, "truth_electron_pt_vs_theta");

        // Signal efficiency is event-level; background survival is candidate-level.
        canvas.Clear(); canvas.SetLogy(false);
        canvas.SetTopMargin(0.34); canvas.SetBottomMargin(0.31);
        auto signal_flow = MakeSurvivalRatio(*cut_study.signal_event_flow[r],1,Form("signal_cutflow_eff_%d",r));
        auto background_flow = MakeSurvivalRatio(*cut_study.background_candidate_flow[r],4,Form("background_cutflow_survival_%d",r));
        signal_flow->SetTitle(";;Fraction of respective baseline");
        signal_flow->SetLineColor(kRed+1); signal_flow->SetMarkerColor(kRed+1); signal_flow->SetMarkerStyle(20); signal_flow->SetLineWidth(2);
        background_flow->SetLineColor(kGray+2); background_flow->SetMarkerColor(kGray+2); background_flow->SetMarkerStyle(21); background_flow->SetLineWidth(2);
        signal_flow->GetXaxis()->LabelsOption("v");
        signal_flow->GetXaxis()->SetLabelSize(0.035);
        signal_flow->Draw("hist e1"); background_flow->Draw("hist e1 same");
        TLegend flow_legend(0.57,0.72,0.91,0.86); flow_legend.SetBorderSize(0); flow_legend.SetFillStyle(0); flow_legend.SetTextSize(0.025);
        flow_legend.AddEntry(signal_flow.get(),"Truth-e event efficiency","lp");
        flow_legend.AddEntry(background_flow.get(),"Background-candidate survival","lp");
        flow_legend.Draw();
        AddEIDPlotLabels(canvas,r,target_lumi_fb);
        SaveStudyCanvas(canvas,png_dir,page,r,"selection_cutflow");

        const char* scan_labels[]={"Track-point threshold scan","Isolation threshold scan","E/(E+H) threshold scan"};
        for (int scan=0; scan<3; ++scan) {
            canvas.Clear(); canvas.SetLogy(false); canvas.SetTopMargin(0.34); canvas.SetBottomMargin(0.16);
            auto signal_scan = MakeEfficiencyRatio(*cut_study.signal_pass[scan][r],*cut_study.signal_total[scan][r],Form("signal_scan_ratio_%d_%d",r,scan));
            auto background_scan = MakeEfficiencyRatio(*cut_study.background_pass[scan][r],*cut_study.background_total[scan][r],Form("background_scan_ratio_%d_%d",r,scan));
            signal_scan->SetTitle(Form(";%s;Efficiency / acceptance",scan_labels[scan]));
            signal_scan->SetLineColor(kRed+1); signal_scan->SetMarkerColor(kRed+1); signal_scan->SetMarkerStyle(20); signal_scan->SetLineWidth(2);
            background_scan->SetLineColor(kGray+2); background_scan->SetMarkerColor(kGray+2); background_scan->SetMarkerStyle(21); background_scan->SetLineWidth(2);
            signal_scan->Draw("hist e1"); background_scan->Draw("hist e1 same");
            TLegend scan_legend(0.57,0.72,0.92,0.86); scan_legend.SetBorderSize(0); scan_legend.SetFillStyle(0); scan_legend.SetTextSize(0.025);
            scan_legend.AddEntry(signal_scan.get(),"Truth-e event efficiency","lp");
            scan_legend.AddEntry(background_scan.get(),"Background-candidate acceptance","lp");
            scan_legend.Draw();
            AddEIDPlotLabels(canvas,r,target_lumi_fb);
            const char* scan_tags[]={"track_point_threshold_scan","isolation_threshold_scan","e_over_eh_threshold_scan"};
            SaveStudyCanvas(canvas,png_dir,page,r,scan_tags[scan]);
        }
        canvas.SetTopMargin(0.1); canvas.SetBottomMargin(0.1);
    }
    std::cout << "Weighted eID study plots\n  input: " << input_name
              << "\n  output ROOT: " << root_name << "\n  output PNG directory: " << png_dir;
    std::cout
              << "\n  target luminosity: " << target_lumi_fb << " fb^-1"
              << "\n  selected Events rows: " << event_rows.size()
              << "\n  skipped rows with unsupported source Q2: " << skipped_regions
              << "\n  candidates with missing event key: " << missing_keys;
    for (int r = 0; r < kInputRegions; ++r)
        std::cout << "\n  minQ2=" << kRegionQ2Min[r] << ": processed="
                  << std::count_if(event_rows.begin(), event_rows.end(), [r](const EventRow& row){ return row.region == r; })
                  << ", N_generated=" << n_generated[r] << ", sigma=" << kSigmaFb[r]
                  << " fb, per-event weight=" << event_weights[r]
                  << (n_generated[r] == 0 ? " (no input supplied)" : "");
    std::cout << "\n  overall combines the processed contributions from all supplied regions.\n";
}

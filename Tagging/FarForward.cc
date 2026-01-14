#include "FarForward.hh"

#include <iostream>

FarForward::FarForward(std::string d_name, std::string b_name) {
    det_name = d_name;
    branch_name = b_name;

    if ( d_name == "rp" )
        det_title = "RP";
    else if ( d_name == "omd" ) 
        det_title = "OMD";
    else
        det_title = "ZDC";

    for ( int i = 0; i < 4; i ++ )
    {
        h_det_xy[i] = NULL;
    }

    h_cluster[0] = NULL;
    h_cluster[1] = NULL;
    h_cluster_sum = NULL;
}

FarForward::~FarForward() {

}

int FarForward::get_n_tracks() 
{
    int n_hit[4] = {0};
    for ( int i = 0; i < 4; i ++ )
        n_hit[i] = det[i]->size();

    int n_protons = *std::min_element(n_hit, n_hit + 4);

    return n_protons;
}

void FarForward::SetEvent(const podio::Frame* event) {
	mEvent = event;
}

void FarForward::set_rotation(double r) {
    rotate = r;
}

void FarForward::set_xShift(int p, double xs) {
    xShift[p] = xs;
}

void FarForward::set_zRange(int p, double min, double max) {
    zRange[p][0] = min;
    zRange[p][1] = max;
}

void FarForward::Process()
{
    if ( det_name == "zdc" ) // should be in a different class, but for now..
        GetEnergy();
    else
        GetHits();

    return;
}

void FarForward::GetHits() 
{
	auto& detHits = mEvent->get<edm4eic::TrackerHitCollection>(branch_name+"RecHits");

    for ( int p = 0; p < 4; p ++ )
    {
        det[p] = new edm4eic::TrackerHitCollection();
        det[p]->setSubsetCollection();
    }
        
	for(const auto& rp : detHits) 
    {
        for ( int pos = 0; pos < 4; pos ++ )
        {
            // hx_raw->Fill(rp.getPosition().x);
            // hz_raw->Fill(rp.getPosition().z);

            if (rp.getPosition().z > zRange[pos][0] && rp.getPosition().z < zRange[pos][1])
                det[pos]->push_back(rp);
        }
	}

    auto& rawAssoc = mEvent->get<edm4eic::MCRecoTrackerHitAssociationCollection>(branch_name + "RawHitAssociations");
    for (const auto& assoc : rawAssoc) 
    {
        std::cout << "rec: " << assoc.getRawHit() << std::endl;
        std::cout << "sim: " << assoc.getSimHit() << std::endl;
    }

	return;
}

double FarForward::GetEnergy() 
{
    is_ZDC_neutron = false;
    ZDC_E.clear();
    ZDC_PBG.clear();

	// auto& detCluster = mEvent->get<edm4eic::ClusterCollection>(branch_name);
    auto& detCluster = mEvent->get<edm4eic::ReconstructedParticleCollection>(branch_name);
 
    double e_sum = 0.0;
    ZDC_E.clear();
	for(const auto& e : detCluster) 
    {
        // cout << "ZDC particle ID " << e.getPDG() << endl;
        // cout << "ZDC energy " << e.getEnergy() << endl;

        ZDC_PBG.push_back(e.getPDG());
        ZDC_E.push_back(e.getEnergy());
    
        if ( e.getPDG() == ID_NEUTRON )
        {
            is_ZDC_neutron = true;
            e_sum += e.getEnergy();
        }
	}

    // std::cout << detCluster.size() << " clusters found in ZDC, total neutron energy: " << e_sum << " GeV" << std::endl;

	return e_sum;
}


void FarForward::define_histograms()
{
    if ( det_name == "zdc" )
    {
        h_cluster[0] = new TH1F(Form("h_%s_cluster_neutron", det_name.c_str()), Form("%s Energy;Cluster Energy[GeV];Counts", det_title.c_str()), 200, 0 , 200);
        h_cluster[1] = new TH1F(Form("h_%s_cluster_gamma", det_name.c_str()), Form("%s Energy;Cluster Energy[GeV];Counts", det_title.c_str()), 200, 0 , 200);
        h_cluster_sum = new TH1F(Form("h_%s_cluster_sum", det_name.c_str()), Form("%s Energy;Cluster Sum[GeV];Counts", det_title.c_str()), 200, 0 , 200);
    }
    else
    {
        for ( int i = 0; i < 4; i ++ )
            h_det_xy[i] = new TH2F(Form("h_%s_xy_%d", det_name.c_str(), i), Form("%s %d hits;x [mm];y [mm]", det_title.c_str(), i), 200, -200, 200, 200, -200, 200);

        // hx_raw = new TH1F(Form("h_%s_x_raw", det_name.c_str()), Form("%s Raw X;X [mm];Counts", det_title.c_str()), 2200, -1400, -800);

        // if ( det_name == "rp" )
        //     hz_raw = new TH1F(Form("h_%s_z_raw", det_name.c_str()), Form("%s Raw Z;Z [mm];Counts", det_title.c_str()), 2000, 32400, 34400);
        // else
        //     hz_raw = new TH1F(Form("h_%s_z_raw", det_name.c_str()), Form("%s Raw Z;Z [mm];Counts", det_title.c_str()), 6000, 22000, 28000);
    }

    return;
}

void FarForward::fill_energy_histograms()
{
    for ( int i = 0; i < ZDC_E.size(); i ++ )
    {
        // std::cout << "ZDC cluster energy: " << ZDC_E[i] << " GeV, PDG: " << ZDC_PBG[i] << std::endl;
        if ( abs(ZDC_PBG[i]) == ID_NEUTRON )
            h_cluster[0]->Fill(ZDC_E[i]);
        else if ( ZDC_PBG[i] == ID_GAMMA )
            h_cluster[1]->Fill(ZDC_E[i]);
    }

    // double sum = 0;
    // for ( auto e : ZDC_E )
    // {
    //     h_cluster->Fill(e);
    //     sum += e;
    // }

    // if ( ZDC_E.size() > 0 )
    //     h_cluster_sum->Fill(sum);

    return;
}

void FarForward::fill_hit_histograms()
{
    for ( int i = 0; i < 4; i ++ )
    {
        for ( const auto &hit : *det[i] )
        {
            double x = -(hit.getPosition().x-xShift[i])*rotate;
            double y = hit.getPosition().y;
            h_det_xy[i]->Fill(x, y);

            // if (det_name == "omd")
                // std::cout << "OMD hit z: " << hit.getPosition().z << " raw x: " << hit.getPosition().x << " x: " << x << " y: " << y << std::endl;
        }
    }

    return;
}

std::vector<TCanvas*> FarForward::draw_histograms()
{
    std::vector<TCanvas*> canvases;

    if ( det_name == "zdc" )
    {
        TCanvas* c_det_e = BookCanvas(Form("c_%s_e", det_name.c_str()), Form("c_%s_e", det_name.c_str()), 800, 500);

        gPad->SetGrid();
        gPad->SetLogy();
        h_cluster[0]->Draw("HIST");
        h_cluster[0]->SetLineColor(kRed);
        h_cluster[1]->Draw("HIST SAME");
        h_cluster[0]->GetYaxis()->SetRangeUser(1, std::max(h_cluster[0]->GetMaximum(), h_cluster[1]->GetMaximum())*1.5);

        TLegend* leg = new TLegend(0.70, 0.70, 0.88, 0.85);
        leg->SetTextSize(0.05);
    	leg->SetBorderSize(0);
    	leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_cluster[0], "n", "L");
        leg->AddEntry(h_cluster[1], "#gamma", "L");
        leg->Draw();

        canvases.push_back(c_det_e);
    }
       
    else
    {
        TCanvas* c_det_xy = BookCanvas(Form("c_%s_xy", det_name.c_str()), Form("c_%s_xy", det_name.c_str()), 1200, 1000);
        c_det_xy->Divide(2,2);

        for ( int i = 0; i < 4 ; i ++ )
        {
            c_det_xy->cd(i+1);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.13);
            gStyle->SetPalette(kLightTemperature);
            
            h_det_xy[i]->SetTitle(Form("%s%02d, %.0f hits", det_title.c_str(), i, h_det_xy[i]->GetEntries()));
            h_det_xy[i]->SetStats(0);
            // h_det_xy[i]->GetListOfFunctions()->Add(exUser);
            h_det_xy[i]->Draw("COLZ");
        }

        canvases.push_back(c_det_xy);

        // TCanvas* c_det_xz = BookCanvas(Form("c_%s_xz", det_name.c_str()), Form("c_%s_xz", det_name.c_str()), 1200, 500);
        // c_det_xz->Divide(2,1); 
         
        // c_det_xz->cd(1);
        // gPad->SetGrid();
        // hx_raw->Draw("HIST");

        // c_det_xz->cd(2);
        // gPad->SetGrid();
        // hz_raw->Draw("HIST");

        // canvases.push_back(c_det_xz);
    }

    return canvases;
}
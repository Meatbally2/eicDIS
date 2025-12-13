#include "DrawManager.hh"

DrawManager::DrawManager(std::string type_, std::string energy_, std::string campaign_) 
: type(type_), energy(energy_), campaign(campaign_) {
}

DrawManager::~DrawManager() {
}

void DrawManager::SetEPIC()
{
    set_ePIC_style();
    gStyle->SetOptTitle(0);
    return;
}

void DrawManager::SetLumi(double L)
{
    lumi = L;
    setLumi = true;
}

void DrawManager::SetQ2range(double Q2_min, double Q2_max)
{
    Q2min = Q2_min;
    Q2max = Q2_max;
    setQ2range = true;
}

void DrawManager::SetQ2min(double Q2_min)
{
    Q2min = Q2_min;
    setQ2min = true;
}

// Modified from ePIC example
void DrawManager::LableAndCollect(TCanvas* &c)
{
    // some loop to figure out number of pads because ROOT doesn't want to give it directly -.-
    int n_pad = 0;
    TList *prim = (TList*)c->GetListOfPrimitives();
    TIter next(prim);
    TObject *obj;
    while ((obj = next())) {
        if (obj->InheritsFrom(TPad::Class()))
            n_pad++;
    }
    int start = (n_pad == 0) ? 0 : 1;
    int end = (n_pad == 0) ? 0 : n_pad;
 
    for ( int i = start; i <= end; i++ )
    {
        if (i == 0)
            c->cd();  // Main canvas
        else
            c->cd(i);  // Sub-pad

        TLatex Text_com;
        Text_com.SetTextAlign(13);  //align at top
        Text_com.DrawLatexNDC(.20,.82,Form("%s, %s", type.c_str(), energy.c_str()));
        
        if ( setLumi )
            Text_com.DrawLatexNDC(.20,.79,Form("L = %.1f fb^{-1}", lumi));

        if ( setQ2range )
            Text_com.DrawLatexNDC(.20,.76,Form("%.0f #leq Q^{2} < %.0f GeV^{2}", Q2min, Q2max));

        if ( setQ2min )
            Text_com.DrawLatexNDC(.20,.76,Form("%.0f #leq Q^{2} GeV^{2}", Q2min));
            
        TLatex Text_ePIC;
        Text_ePIC.SetTextSize(0.05);
        Text_ePIC.SetTextFont(62);
        Text_ePIC.DrawLatexNDC(.29,.86,"Performance");  // performance plot
        //Text_ePIC.DrawLatexNDC(.35,.88,"ePIC Internal");  // for internal use only
        //Text_ePIC.DrawLatexNDC(.35,.88,"ePIC Preliminary"); // preliminary released version 
        //Text_ePIC.DrawLatexNDC(.35,.88,"ePIC Work in Progress"); // work in progress to be shown outside
        //Text_ePIC.DrawLatexNDC(.35,.88,"ePIC"); // final published version
        
        // Add dates: needed for performance plots
        TLatex Text_date;
        Text_date.SetTextSize(0.035);
        Text_date.SetTextFont(52);
        Text_date.DrawLatexNDC(.65,.96,Form("Simulation campaign: %s", campaign.c_str()));  // performance plot

        Text_com.Draw();
        Text_ePIC.Draw();
        Text_date.Draw();

        // ===== Add ePIC logo to the figure ======
        TImage *logo = TImage::Open("../GlobalUtil/EPIC-logo_black_transparent.png");
        // TPad *logo_pad = new TPad("logo_pad", "logo_pad", 0.8, 0.8, 0.93, 0.93); // Create a new pad and then draw the image in it
        TPad *logo_pad = new TPad("logo_pad", "logo_pad", 0.17, 0.83, 0.33, 0.93); // Create a new pad and then draw the image in it
        logo_pad->SetFillStyle(0);
        logo_pad->SetFillColor(0);
        logo_pad->SetFrameFillStyle(0);
        logo_pad->SetBorderMode(0);
        logo_pad->SetBorderSize(0);
        logo_pad->SetFrameBorderMode(0);
        logo_pad->SetLineWidth(0);
        logo_pad->SetLineColor(0);
        logo_pad->SetFrameLineWidth(0);
        logo_pad->SetFrameLineColor(0);
 
        TPad *pad_logo = (TPad*)logo_pad->Clone(Form("logo_pad_%d", i));
        logo->SetConstRatio(kTRUE);  // Maintain aspect ratio
        pad_logo->Draw("X");
        pad_logo->cd();
        logo->Draw();
    }
    
    c->Modified();
    c->Update();
    canvas_list.push_back(c);

    return;
}

void DrawManager::SaveToTree(TFile* &outFile)
{
    outFile->cd();
    for ( auto c : canvas_list )
        c->Write();

    return;
}
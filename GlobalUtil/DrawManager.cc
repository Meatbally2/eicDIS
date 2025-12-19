#include "DrawManager.hh"

DrawManager::DrawManager(std::string type_, std::string energy_, std::string campaign_) 
: type(type_), energy(energy_), campaign(campaign_) {
}

DrawManager::DrawManager(std::string type_, std::string campaign_) 
: type(type_), campaign(campaign_) {
    multiple_energies = true;
    std::cout << "** Multiple energies mode enabled in DrawManager **" << std::endl;
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

void DrawManager::SetDISrange(double Y_min, double Y_max, double W2_min, double Q2_min)
{
    Ymin = Y_min;
    Ymax = Y_max;
    W2min = W2_min;
    Q2min = Q2_min;
    setLumi = true;
    setDIS = true;
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

void DrawManager::SetPlotType(std::string plot_type_)
{
    plot_type = plot_type_;
}

void DrawManager::SaveToTree(TFile* &outFile)
{
    outFile->cd();
    for ( auto c : canvas_list )
        c->Write();

    return;
}

// Modified from ePIC example
void DrawManager::LableAndCollect(TCanvas* &c, int draw_position = 0)
{
    // Position code
    // 0: top-left (default)
    // 1: top-right 
    // 2: bottom-right
    // Not working very well for canvas with multiple pads yet

    std::string lumi_unit = type == "ep" ? "fb^{-1}" : "fb^{-1}/A";

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

       // Get canvas dimensions for scaling
        Double_t canvas_width = gPad->GetWw();
        Double_t canvas_height = gPad->GetWh();

        // Reference canvas size (your standard size)
        Double_t ref_width = 1398.0;
        Double_t ref_height = 575.0;

        // Calculate scale factor (use minimum to ensure logo fits)
        Double_t width_scale_factor = canvas_width / ref_width;
        Double_t height_scale_factor = canvas_height / ref_height;
        Double_t scale_factor = TMath::Min(width_scale_factor, height_scale_factor);

        // if (scale_factor < 1.0) {
            // Small canvas: use a less aggressive scaling (square root or custom formula)
            // scale_factor = TMath::Sqrt(scale_factor);  // Less shrinkage
            // scale_factor = 0.7 + 0.3 * scale_factor;  // Never go below 70%
        // }

        Double_t left_margin = 0.19;
        Double_t top_margin = 0.93;
        if ( draw_position == 1 ) {
            left_margin = 0.5;
        }
        if ( draw_position == 2 ) {
            left_margin = 0.6;
            top_margin = 0.42;
        }

        // ===== Add ePIC logo to the figure ======
        TImage *logo = TImage::Open("../GlobalUtil/EPIC-logo_black_transparent.png");
        UInt_t img_width = logo->GetWidth();    // Image width in pixels
        UInt_t img_height = logo->GetHeight();  // Image height in pixels
        Double_t img_aspect = (Double_t)img_width / (Double_t)img_height;

        // Base logo height in NDC (at reference canvas size)
        Double_t base_logo_height_ndc = 0.17;

        // Scale the logo height based on canvas size
        Double_t logo_height = base_logo_height_ndc * scale_factor;

        // Calculate logo width in NDC to maintain image aspect ratio
        // NDC is not square, so we need to account for canvas aspect ratio
        Double_t canvas_aspect = canvas_width / canvas_height;
        Double_t logo_width = logo_height * img_aspect / canvas_aspect;

        // Create and draw logo pad
        TPad *logo_pad = new TPad("logo_pad", "logo_pad", 
                                left_margin, top_margin - logo_height, 
                                left_margin + logo_width, top_margin);
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
        logo_pad->SetLeftMargin(0.0);
        logo_pad->SetBottomMargin(0.0);
        logo_pad->SetRightMargin(0.0);
        logo_pad->SetTopMargin(0.0);

        logo_pad->Draw();
        logo_pad->cd();
        logo->SetConstRatio(kTRUE);
        logo->Draw();

        // Back to main pad for text
        if (i == 0)
            c->cd();
        else
            c->cd(i);

        // Position "Internal" text right after the logo
        Double_t text_gap = 0.015 * scale_factor;  // Gap also scales
        Double_t internal_x = left_margin + logo_width + text_gap;
        Double_t internal_y = top_margin - 0.074 * scale_factor;

        if ( draw_position == 1 ) {
            internal_y = top_margin - 0.076 * scale_factor;
        }

        TLatex Text_ePIC;
        Text_ePIC.SetTextSize(0.065 * scale_factor);  // Text size also scales
        Text_ePIC.SetTextFont(62);
        Text_ePIC.SetTextAlign(13);  // Top-left alignment
        // Text_ePIC.DrawLatexNDC(internal_x, internal_y, "Internal");
        if ( draw_position == 1 ) {
            Text_ePIC.SetTextSize(0.075 * scale_factor);
        }
        Text_ePIC.DrawLatexNDC(internal_x, internal_y, "Performance");

        // Scale all other text elements too
        TLatex Text_com;
        Text_com.SetTextSize(0.055 * scale_factor);
        Text_com.SetTextAlign(13);

        // Base positions and spacing
        // Double_t text_x = 0.195;
        Double_t text_x = left_margin + 0.005;
        Double_t text_y_start = internal_y - 0.12*scale_factor;
        Double_t line_spacing = 0.07 * scale_factor;  // Scale the line spacing
        Double_t current_y = text_y_start;

        if ( draw_position == 1 ) {
            text_x = 0.605;
        }

        if ( !multiple_energies )
            Text_com.DrawLatexNDC(text_x, current_y, Form("%s, %s", type.c_str(), energy.c_str()));
        else
            Text_com.DrawLatexNDC(text_x, current_y, Form("%s, L = %.1f %s", type.c_str(), lumi, lumi_unit.c_str()));

        if ( setDIS )
        {
            current_y -= line_spacing;
            Text_com.DrawLatexNDC(text_x, current_y, Form("%.0f #leq y < %.0f", Ymin, Ymax));
            current_y -= line_spacing;
            Text_com.DrawLatexNDC(text_x, current_y, Form("Q^{2} #geq %.0f GeV^{2}", Q2min));
            current_y -= line_spacing;
            Text_com.DrawLatexNDC(text_x, current_y, Form("W^{2} #geq %.0f GeV^{2}", W2min));
        }

        if ( setLumi && !setDIS && !multiple_energies )
        {
            current_y -= line_spacing;
            Text_com.DrawLatexNDC(text_x, current_y, Form("L = %.1f %s", lumi, lumi_unit.c_str()));
        }

        if ( setQ2range )
        {
            current_y -= line_spacing;
            Text_com.DrawLatexNDC(text_x, current_y, Form("%.0f #leq Q^{2} < %.0f GeV^{2}", Q2min, Q2max));
        }

        if ( setQ2min )
        {
            current_y -= line_spacing;
            Text_com.DrawLatexNDC(text_x, current_y, Form("%.0f #leq Q^{2} GeV^{2}", Q2min));
        }

        TLatex Text_date;
        Text_date.SetTextSize(0.035 * scale_factor);
        Text_date.SetTextFont(52);
        Text_date.SetTextAlign(31);  // Right-aligned, bottom of text at y position
        if ( draw_position == 1 ) {
            Text_date.SetTextSize(0.045 * scale_factor);
        }

        // Get the right edge of the frame (accounting for right margin)
        Double_t right_margin = gPad->GetRightMargin();
        Double_t date_x = 1.0 - right_margin - 0.005;  // Right edge minus small padding

        Double_t date_y = 0.955;  // Near top
        Text_date.DrawLatexNDC(date_x, date_y, Form("Simulation campaign: %s", campaign.c_str()));  // performance plot

        Text_com.Draw();
        Text_ePIC.Draw();
        Text_date.Draw();
    }
    
    c->Modified();
    c->Update();
    canvas_list.push_back(c);

    return;
}
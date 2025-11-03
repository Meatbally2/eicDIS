
double x_gen_min = 0;
double q_gen_max = 0;

void BinLogX(TH1F* &h) // taken from https://root.cern.ch/root/roottalk/roottalk06/1213.html
{
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) 
    new_bins[i] = pow(10, from + i * width);

  axis->Set(bins, new_bins);
  delete[] new_bins;
}

void BinLogX(TH2F* &h, int a) // taken from https://root.cern.ch/root/roottalk/roottalk06/1213.html
{
  TAxis *axis = a == 0 ? h->GetXaxis() : h->GetYaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) 
    new_bins[i] = pow(10, from + i * width);

  axis->Set(bins, new_bins);
  delete[] new_bins;
}

void BinLogXY(TH2F* &h)
{
  for ( int j = 0; j < 2; j ++ )
  {
    BinLogX(h, j);
    // TAxis *axis = j == 0 ? h->GetXaxis() : h->GetYaxis();

    // int bins = axis->GetNbins();
    // Axis_t from = axis->GetXmin();
    // Axis_t to = axis->GetXmax();
    // Axis_t width = (to - from) / bins;
    // Axis_t *new_bins = new Axis_t[bins + 1];

    // for (int i = 0; i <= bins; i++) 
    //   new_bins[i] = pow(10, from + i * width);

    // axis->Set(bins, new_bins);

    // delete[] new_bins;
  }
}

TH1F* BookTH1(std::string name, std::string title, double n_bins, double min, double max)
{
  TH1F* h = new TH1F(name.c_str(), title.c_str(), n_bins, min, max);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.5);
  h->GetYaxis()->SetTitleOffset(1.5);

  BinLogX(h);

  return h;
}

TH2F* BookTH2(std::string name, std::string title, double n_xbins, double xmin, double xmax, double n_ybins, double ymin, double ymax, int colorMap)
{
  TH2F* h = new TH2F(name.c_str(), title.c_str(), n_xbins, xmin, xmax, n_ybins, ymin, ymax);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleOffset(1.2);

  TExec *exDefault = new TExec("ex1","gStyle->SetPalette(kBird)");
  TExec *exUser = new TExec("ex1",Form("gStyle->SetPalette(%d)",colorMap));

  if ( colorMap != 0 )
    h->GetListOfFunctions()->Add(exUser);
  else
    h->GetListOfFunctions()->Add(exDefault);

  BinLogXY(h);

  return h;
}

TCanvas* BookCanvas(std::string name, std::string title, double width, double height)
{
  TCanvas* c = new TCanvas(name.c_str(), title.c_str(), width, height);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.13);

  return c;
}

void plot_2d_xq(std::string c_name, std::string c_title, double c_min, double c_max, TH2F* h_xq2, double x_gen_min, double q_gen_max)
{
    TCanvas* c_xq = BookCanvas(c_name, c_title, c_min, c_max);
    c_xq->SetLogx();
    c_xq->SetLogy();
    // c_xq->SetLogz();
    h_xq2->SetStats(0);
    h_xq2->GetXaxis()->SetRangeUser(x_gen_min, 1);
    h_xq2->GetYaxis()->SetRangeUser(1,q_gen_max);
    h_xq2->Draw("COLZ TEXT");

    return;
}

void format_efficiency(TEfficiency* &g)
{
    auto graph_test = g->GetPaintedGraph();
    if ( g->GetPaintedGraph() == NULL )
        printf("cannot find \n");

    g->GetPaintedGraph()->GetXaxis()->CenterTitle();
    g->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);
   
    g->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.8);
    g->GetPaintedGraph()->GetYaxis()->CenterTitle();

    g->GetPaintedGraph()->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g->GetPaintedGraph()->GetXaxis()->SetLabelSize(20);
    g->GetPaintedGraph()->GetXaxis()->SetTitleSize(20);
    g->GetPaintedGraph()->GetXaxis()->SetTitleFont(43);

    g->GetPaintedGraph()->GetYaxis()->SetNdivisions(505);
    g->GetPaintedGraph()->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g->GetPaintedGraph()->GetYaxis()->SetLabelSize(20);
    g->GetPaintedGraph()->GetYaxis()->SetTitleSize(20);
    g->GetPaintedGraph()->GetYaxis()->SetTitleFont(43);

    return;
}

void format_histogram(TH1F* &h)
{
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.2);
   
    h->GetYaxis()->SetTitleOffset(1.8);
    h->GetYaxis()->CenterTitle();

    h->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetXaxis()->SetLabelSize(20);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleFont(43);

    h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetYaxis()->SetLabelSize(20);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleFont(43);

    return;
}

TCanvas* draw_2d_standard(TH2F* &h, std::string c_name, std::string c_title, double width, double height, bool is_log_z, bool is_std_scale)
{
    TCanvas* c = BookCanvas(c_name.c_str(), c_title.c_str(), width, height);
    c->SetLogx();
    c->SetLogy();
    h->SetStats(0);
    h->Draw("COLZ");

    if ( is_log_z )
        c->SetLogz();

    if (is_std_scale )
    {
        h->GetXaxis()->SetRangeUser(x_gen_min, 1);
        h->GetYaxis()->SetRangeUser(1,q_gen_max);
    }

    return c;
}

TCanvas* draw_2d_efficiency(TH2F* &h, std::string c_name, std::string c_title, double width, double height, bool is_log_z, bool is_std_scale)
{
    TCanvas* c = BookCanvas(c_name.c_str(), c_title.c_str(), width, height);
    c->SetLogx();
    c->SetLogy();
    h->SetStats(0);
    h->Draw("COLZ TEXT");

    if ( is_log_z )
        c->SetLogz();

    if (is_std_scale )
    {
        h->GetXaxis()->SetRangeUser(x_gen_min, 1);
        h->GetYaxis()->SetRangeUser(1,q_gen_max);
    }

    return c;
}

void set_2d_scale(TH2F* h)
{
    x_gen_min = h->GetXaxis()->GetBinCenter(h->FindFirstBinAbove(0, 1)) - h->GetXaxis()->GetBinWidth(h->FindFirstBinAbove(0, 1)/2.);
    q_gen_max = h->GetYaxis()->GetBinCenter(h->FindLastBinAbove(0, 2)) + h->GetYaxis()->GetBinWidth(h->FindLastBinAbove(0, 2)/2.);
    
    return;
}

void process_eff_hist(TH2F* &h_top, TH2F* &h_bottom)
{
    h_top->Divide(h_bottom);
    for ( int ix = 0; ix < h_top->GetXaxis()->GetNbins(); ix ++ )
	{
		for ( int iq = 0; iq < h_top->GetYaxis()->GetNbins(); iq ++ )
		{
            double eff = h_top->GetBinContent(ix+1, iq+1);
            double eff_round = std::round(eff * 1000) / 1000.;
            if ( eff_round > 1.5 )
                eff_round = 1.5;
            // h_top->SetBinContent(ix+1, iq+1, std::round(eff * 1000) / 1000.);
            // if ( eff_round < 2 )
                h_top->SetBinContent(ix+1, iq+1, eff_round);
            // else
            //     h_top->SetBinContent(ix+1, iq+1, 0);
        }
    }

    return;
}
#include <ctime>

void d_meson_ls_pt(TString filename = "results_all.root", int rebin = 10, int down_bin = 30, int up_bin = 120, int cond = 1) {
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	TFile *_file0 = TFile::Open(filename);
	TH2D *h_inv_mass_US_Kminus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_US_Kminus_vs_pt"));
	TH2D *h_inv_mass_US_Kplus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_US_Kplus_vs_pt"));
	TH2D *h_inv_mass_LS_Kminus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_LS_Kminus_vs_pt"));
	TH2D *h_inv_mass_LS_Kplus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_LS_Kplus_vs_pt"));
	
	h_Dstar_WS_vs_pt->ProjectionX("proj_ws", down_bin, up_bin);
	h_Dstar_RS_vs_pt->ProjectionX("proj_rs", down_bin, up_bin);
	h_Dstar_Sideband_vs_pt->ProjectionX("proj_sb", down_bin, up_bin);
	double int_rs = proj_rs->Integral(150,300);
	double int_sb = proj_sb->Integral(150,300);
	/*
	double int_rs = proj_rs->Integral();
	double int_sb = proj_sb->Integral();
	*/
	cout << "int_rs: " << int_rs << " int_sb: " << int_sb << endl;
	proj_ws->Rebin(rebin);
	proj_rs->Rebin(rebin);
	proj_sb->Rebin(rebin);
	const TH1D *proj_ws2 = proj_ws->Clone();
	proj_ws->Draw();
	proj_ws->SetTitle("");
	proj_ws->GetXaxis()->SetTitle("M_{K#pi#pi}-M_{K#pi}  [GeV/c^{2}]");
	proj_ws->GetYaxis()->SetTitle("Counts");
	proj_rs->SetLineColor(kBlue);
	proj_rs->GetXaxis()->SetTitle("M_{K#pi#pi}-M_{K#pi}  [GeV/c^{2}]");
	proj_rs->GetYaxis()->SetTitle("Counts");
	//proj_rs->SetTitle("Wrong sign");
	//proj_rs->SetTitle("RS - WS");
	proj_rs->SetTitle("");
	//proj_rs->Draw("sames");
	proj_rs->Draw();
	proj_ws->Draw("sames");
	proj_sb->Scale(int_rs/int_sb);
	proj_sb->SetLineColor(kGreen);
	proj_sb->Draw("sames");
	//proj_rs->Add(proj_ws2, -1);
	//proj_rs->Draw();
	TLegend *leg = new TLegend(0.6, 0.3, 0.8, 0.4);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(1);
	leg->SetShadowColor(0);
	leg->SetTextFont(22);
	leg->SetTextSize(0.04);
	//leg->SetHeader("Unlike sign K#pi pairs");
	leg->AddEntry(proj_rs,"RS","Pl");
	leg->AddEntry(proj_ws,"WS","Pl");
	leg->AddEntry(proj_sb,"SB","Pl");
	leg->Draw();


	TString output_file = "pics/"+filename;
	if (cond == 1) {
		output_file.ReplaceAll(".root","_dstar_Kpipi_Kpi.pdf");
	} else if (cond == 2) {
		output_file.ReplaceAll(".root","_dstar_Kpipi_Kpi_soft.pdf");
	}
	output_file.ReplaceAll(".pdf",Form("_%d_%d_%d.pdf",down_bin,up_bin,rebin));
	c2->Print(output_file);
	output_file.ReplaceAll(".pdf",".png");
	c2->Print(output_file);

	/*
	const TH1D *proj_rs2 = proj_rs->Clone();
	proj_rs2->Add(proj_sb, -1);
	//proj_rs->Rebin(rebin);
	   proj_rs2->SetTitle("RS - SB");
	   proj_rs2->Draw();
	   */
	/*
	float min_pt = down_bin/10.;
	float max_pt = up_bin/10.;
	   */
	//TCanvas *c4 = draw(signal, rebin, min_pt, max_pt, output_file_txt);
	//TCanvas *c4 = draw(signal, rebin, min_pt, max_pt, output_file_txt, inc_stat, Ntrig);
	//TCanvas *c_sb = draw(proj_rs2, rebin);
	const int n_pt = 3;
	//float bin_array[] = {0., 2., 3., 4., 5., 6.};
	//float bin_array[] = {1., 2., 3., 4., 5., 6.};
	float bin_array[] = {1., 2., 4., 6.};
	TLatex tx1;
	tx1.SetNDC();
	
	time_t  timev;
	time(&timev);
	std::stringstream ss;
	ss << timev;
	TString s_array = Form("%d", timev);
	
	TString sts = "_" + s_array + "_sb.txt";
	TString output_file_txt_sb = filename;
	output_file_txt_sb.ReplaceAll(".root", sts);
	output_file_txt_sb = "results_txt/" + output_file_txt_sb;
	ofstream results_file_sb (output_file_txt_sb,  ios::trunc);
	results_file_sb.close();
	cout << "output_file_txt_sb: " << output_file_txt_sb << endl;

	TString sts = "_" + s_array + "_ws.txt";
	TString output_file_txt_ws = filename;
	output_file_txt_ws.ReplaceAll(".root", sts);
	output_file_txt_ws = "results_txt/" + output_file_txt_ws;
	ofstream results_file_ws (output_file_txt_ws,  ios::trunc);
	results_file_ws.close();
	cout << "output_file_txt_ws: " << output_file_txt_ws << endl;
	
	/*
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	//c2->Divide(2,2);
	c2->Print(output_file+"[");
	*/
	for (int i_pt = 0; i_pt < n_pt; i_pt++) {
		float min_pt = bin_array[i_pt];
		float max_pt = bin_array[i_pt+1];
		/*
		float min_pt = bin_array[i_pt]*rebin/1000.;
		float max_pt = bin_array[i_pt+1]*rebin/1000.;
		float min_pt = down_bin/10.;
		float max_pt = up_bin/10.;
		*/
		cout << "pt: " << min_pt << " - " << max_pt << endl;

		down_bin = (int)min_pt*100;
		up_bin = (int)max_pt*100;

		TH1D *proj_us_Kminus = h_inv_mass_US_Kminus_vs_pt->ProjectionX("proj_us_Kminus", down_bin, up_bin);
		TH1D *proj_us_Kplus = h_inv_mass_US_Kplus_vs_pt->ProjectionX("proj_us_Kplus", down_bin, up_bin);
		TH1D *proj_ls_Kminus = h_inv_mass_LS_Kminus_vs_pt->ProjectionX("proj_ls_Kminus", down_bin, up_bin);
		TH1D *proj_ls_Kplus = h_inv_mass_LS_Kplus_vs_pt->ProjectionX("proj_ls_Kplus", down_bin, up_bin);
		
		TH1D *proj_us = proj_us_Kminus->Clone();
		proj_us->Add(proj_us_Kplus, 1);

		TH1D *proj_ls = proj_ls_Kminus->Clone();
		proj_ls->Add(proj_ls_Kplus, 1);

		int int_down = 1500;
		int int_up = 1700;
		double int_us = proj_us->Integral(int_down, int_up);

		cout << "int_rs: " << int_rs << " int_sb: " << int_sb << endl;
		proj_ls->Rebin(rebin);
		proj_us->Rebin(rebin);
		TH1D *proj_ls2a = proj_ls->Clone();
		proj_ls->Draw();
		proj_ls->SetTitle("");
		proj_ls->GetXaxis()->SetTitle("M_{K#pi}  [GeV/c^{2}]");
		proj_ls->GetYaxis()->SetTitle("Counts");
		proj_us->SetLineColor(kBlue);
		proj_us->GetXaxis()->SetTitle("M_{K#pi}  [GeV/c^{2}]");
		proj_us->GetYaxis()->SetTitle("Counts");
		//proj_us->SetTitle("Wrong sign");
		//proj_us->SetTitle("RS - WS");
		proj_us->SetTitle("");
		//proj_us->Draw("sames");
		proj_us->Draw();
		proj_ls->Draw("sames");
		const TH1D *proj_us2 = proj_us->Clone();
		proj_us->Add(proj_ls2a, -1);

		TCanvas *c_ws = draw(proj_us, rebin, min_pt, max_pt, output_file_txt_ws);
		c_ws->cd();
		tx1.DrawLatex(0.1,0.94,Form("D^0 p_{T}: %3.1f-%3.1f [GeV/c]", min_pt, max_pt));
		tx1.DrawLatex(0.15,0.8,"RS-WS");

		/*
		   proj_us->SetTitle("RS - WS");
		   proj_us->Draw();
		   */
		output_file = "pics/"+filename;
		if (cond == 1) {
			output_file.ReplaceAll(".root","_dstar_D0.pdf");
		} else if (cond == 2) {
			output_file.ReplaceAll(".root","_dstar_Kpipi_Kpi_RS_min_WS_soft.pdf");
		}
		output_file.ReplaceAll(".pdf",Form("_%d_%d_%d.pdf",down_bin,up_bin,rebin));
		c_ws->Print(output_file);
		output_file.ReplaceAll(".pdf",".png");
		c_ws->Print(output_file);

	}
	//c2->Print(output_file+"]");
	//cout << "output_file_txt: " << output_file_txt << endl;

}

//TCanvas *draw(TH1D *signal, int rebin_h = 1, double low_pt = 1.0, double high_pt = 2.0, TString output_file_txt = "def.txt", float inc_stat = 1, double Ntrig = 41E6) {
TCanvas *draw(TH1D *signal, int rebin_h = 1, double low_pt = 1.0, double high_pt = 2.0, TString output_file_txt = "def.txt", float inc_stat = 1, double Ntrig = 41E6) {
	const float fitRMin   = 1.72;
	//const float fitRMin   = 1.7;
	const float fitRMax   = 2.05;
	const float MKSize    = 1.;
	const float rotwthmin = 1.857;
	//const float rotwthmin = 1.85;
	const float rotwthmax = 1.89;

	 TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	TF1 *fun0 = new TF1("fun0","pol1(0)+gaus(2)",fitRMin,fitRMax);
	fun0->SetParameters(1.,1.,1.,1.865,0.015);
	fun0->SetLineColor(2);
	fun0->SetLineStyle(7);
	const int N = 400/rebin_h;
	const int binm1 = 1700/rebin_h;
	const int binm2 = 185;
	const int binm3 = 192;
	const int binm4 = 2100/rebin_h;

	/*
	const int N = 400/rebin_h;
	const int binm1 = 0;
	const int binm4 = 400/rebin_h;
	*/
	const double scalem = 1.;
	double mm[N],ym[N],yme[N],ym1[N];

	/*
	const int N = 400/rebin_h;
	//const int binm1 = 80/rebin_h;
	const int binm1 = 0;
	const double scalem = 1.;
	*/



	cout << signal-> GetNbinsX() << endl;
	int i=0;
	for(int ib=binm1; ib<binm4; ib++) {
		mm[i]  = signal->GetBinCenter(ib);
		ym[i]  = signal->GetBinContent(ib) + scalem;
		yme[i] = signal->GetBinError(ib);
		yme[i] = yme[i]/TMath::Sqrt(inc_stat);
		i++;
	}


	TGraphErrors *gm = new TGraphErrors(N,mm,ym,0,yme);
	gm->SetMarkerStyle(20);
	gm->SetMarkerSize(MKSize);
	gm->SetMarkerColor(4);
	gm->SetLineColor(4);
	//gm->Draw("psame");
	//gm->Draw("apsame");
	//gm->Draw("AP");
	fun0->SetParLimits(3,rotwthmin,rotwthmax);
	//fun0->SetParLimits(4,0.010,0.014);
	//fun0->SetParLimits(4,0.007,0.020);
	fun0->SetParLimits(4,0.007,0.015);
	//  signal->Fit(fun0,"NOR");
	gm->Fit(fun0,"NOR");

	TF1 *resfunm = new TF1("resfunm","pol1",fitRMin,fitRMax);
	resfunm->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1));
	resfunm->SetLineStyle(7);
	resfunm->SetLineWidth(1);
	//resfunm->Draw("same");
	//fun0->Draw("same");

	TF1 *resfunm1 = new TF1("resfunm1","pol1(0)+gaus(2)",fitRMin,fitRMax);
	resfunm1->SetParameters(0.,0.,fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
	resfunm1->SetLineColor(4);
	resfunm1->SetLineStyle(7);
	//resfunm1->Draw("same");

	for(int i=0; i<N; i++) ym1[i] = ym[i] - resfunm->Eval(mm[i]);

	TGraphErrors *gm1 = new TGraphErrors(N,mm,ym1,0,yme);
	gm1->SetMarkerStyle(24);
	gm1->SetMarkerSize(MKSize-0.1);
	gm1->SetMarkerColor(2);
	gm1->SetLineColor(2);
	//gm1->Draw("psame");

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(gm);
	mg->Add(gm1);
	mg->Draw("ap");
	//mg->GetYaxis()->SetTitle("Raw Counts (#times 10^{3})");
	mg->GetYaxis()->SetTitle("Raw Counts");
	//mg->GetYaxis()->SetTitleSize(titsize);
	mg->GetYaxis()->SetTitleOffset(1.5);
	//mg->GetYaxis()->SetLabelOffset(0.03);
	//mg->GetYaxis()->SetLabelSize(0.047);
	//mg->GetXaxis()->SetNdivisions(208);
	//  mg->GetXaxis()->CenterTitle();
	mg->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
	resfunm->Draw("same");
	fun0->Draw("same");
	resfunm1->Draw("same");

	TLatex tx1;
	tx1.DrawLatex(0.17,0.84,"SE-ME");

	float chi2m   = fun0->GetChisquare();
	float ndfm    = fun0->GetNDF();
	float meanm   = fun0->GetParameter(3);
	float meanme  = fun0->GetParError(3);
	float sigm    = fun0->GetParameter(4);
	float sigme   = fun0->GetParError(4);

	float width_bin = 0.001 * rebin_h;
	float yieldm  = sqrt(6.283)*fun0->GetParameter(2)*sigm/width_bin;
	float yieldme = yieldm*fun0->GetParError(2)/fun0->GetParameter(2);
	float bckg    = resfunm->Integral((fun0->GetParameter(3)-5*sigm), (5*sigm+fun0->GetParameter(3)))/width_bin;

	//float bckg    = resfunm->Integral((fun0->GetParameter(4)-3*sigm), (3*sigm+fun0->GetParameter(4)))/0.01;
	//float bckg    = resfunm->Integral((fun0->GetParameter(4)-sigm), (sigm+fun0->GetParameter(4)))/0.01;
	//float bckg    = resfunm->Integral((fun0->GetParameter(4)-3*sigm/0.01), (3*sigm/0.01+fun0->GetParameter(4)));
	//float bckg    = resfunm->Integral(1.8, 1.9);
	cout << " (fun0->GetParameter(3)-5*sigm): " << (fun0->GetParameter(3)-5*sigm) << " 5*sigm+fun0->GetParameter(3)): " << 5*sigm+fun0->GetParameter(3) << endl;
	float nsig    = yieldm / TMath::Sqrt(yieldm + bckg);
	cout << "sigm: " << sigm << " fun0->GetParameter(2): " << fun0->GetParameter(2) << endl;
	cout << " fun0->GetParameter(0): " << fun0->GetParameter(0) << endl;
	cout << " fun0->GetParameter(1): " << fun0->GetParameter(1) << endl;
	cout << " fun0->GetParameter(2): " << fun0->GetParameter(2) << endl;
	cout << " fun0->GetParameter(3): " << fun0->GetParameter(3) << endl;
	cout << " fun0->GetParameter(4): " << fun0->GetParameter(4) << endl;
	cout << " fun0->GetParameter(5): " << fun0->GetParameter(5) << endl;
	//cout << " fun0->GetParameter(6): " << fun0->GetParameter(6) << endl;
	cout << "bckg: " << bckg << " sig: " << yieldm << " N_sig: " << nsig << endl;

	char chh1[250],chh2[250],chh3[250],chh4[250], chh5[250], chh6[250];

	float yieldmBin = 0;
	float yieldmeBin = 0;
	int yieldmBinMin = 1850/rebin_h - binm1;
	int yieldmBinMax = 1900/rebin_h - binm1;
	for (int i = yieldmBinMin; i < yieldmBinMax; i++) {
		yieldmBin += ym1[i];
		yieldmeBin += yme[i]*yme[i];
		cout << ym1[i] << " " << yme[i] << endl;
	}
	yieldmeBin = TMath::Sqrt(yieldmeBin);
	cout << endl;
	cout << "Bin yield: " << yieldmBin << " +- " << yieldmeBin << endl;

	sprintf(chh1,"#chi^{2}/ndf = %4.1f / %d",chi2m,ndfm);
	sprintf(chh2,"yield = %7.0f #pm %7.0f",yieldm,yieldme);
	sprintf(chh3,"mean = %5.3f #pm %5.3f",meanm,meanme);
	sprintf(chh6,"yieldBin = %7.0f #pm %7.0f",yieldmBin,yieldmeBin);
	sprintf(chh4,"#sigma = %5.3f #pm %5.3f",sigm,sigme);
	sprintf(chh5,"N_sig = %2.1f", nsig);
	//sprintf(chh4,"N_sig = %2.1f #pm %2.1f", nsig, nsige);

	//TString output_file_txt = "def.txt";
	ofstream results_file (output_file_txt,  ios::app);
	if (results_file.is_open()) {
		//results_file << endl;
		//results_file << __DATE__ << "  " << __TIME__ << "    " << endl;
	} else cout << "Unable to open file.";
	results_file.precision(5);
	/*
	double Eff = 0.05;
	double pt = (high_pt + low_pt)/2.;
	double pt_bin_width = high_pt - low_pt;
	double inv_yield = yieldm / (4*TMath::Pi() * Ntrig * 0.0389 * pt_bin_width * pt * 2 * Eff);
	cout << "inv_yield = " << yieldm << " / (4 * TMath::Pi() * " << Ntrig << " * 0.0389 * " << pt_bin_width << " * " << pt << " * 2 * "<< Eff << ") = " << inv_yield << endl;
	double inv_yield_err = yieldme / (4*TMath::Pi() * Ntrig * 0.0389 * pt_bin_width * pt * 2 * Eff);
	cout << "Invariant Yield: " << low_pt << " " << high_pt << " " << inv_yield << " " << inv_yield_err << endl;
	//results_file << low_pt << " " << high_pt << " " << inv_yield << " " << inv_yield_err << endl;
	results_file << low_pt << "\t" << high_pt << "\t" << inv_yield << "\t" << inv_yield_err << endl;
	*/
	cout << "Yield: " << low_pt << " " << high_pt << " " << yieldm << " " << yieldme << endl;
	results_file << low_pt << " " << high_pt << " " << yieldm << " " << yieldme << endl;

	if (results_file.is_open()) {
		results_file.close();
	} else cout << "File is not open.";


	TPaveStats *pave1 = new TPaveStats(0.58,0.73,0.98,0.93,"brNDC");
	//TPaveStats *pave1 = new TPaveStats(0.48,0.78,0.98,0.95,"brNDC");
	pave1->SetName("stats");
	pave1->SetBorderSize(0);
	pave1->SetFillColor(10);
	//    pave1->SetTextAlign(12);
	pave1->SetTextFont(22);
	pave1->SetTextSize(0.04);
	TText *text1 = pave1->AddText(chh1);
	text1 = pave1->AddText(chh2);
	text1 = pave1->AddText(chh6);
	text1 = pave1->AddText(chh3);
	text1 = pave1->AddText(chh4);
	//text1 = pave1->AddText(chh5);
	pave1->Draw();


	c3->Update();
	return c3;


}

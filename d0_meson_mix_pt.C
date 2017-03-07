//for any backgroung estimations, currently mix = ls
#include <ctime>

void d0_meson_mix_pt(TString filename = "results_all.root", float inc_stat = 1, double Ntrig = 41E6) {
	TLatex tx1;
	tx1.SetNDC();

	const int n_pt = 3;
// 	int bin_array[] = {50,110,220,500};
	int bin_array[] = {20,120,220, 500};
	//int bin_array[] = {20,60,110,200};
	//int bin_array[] = {50,110,220,300};
	//int bin_array[] = {60,120,200,240};
	int rebin = 20;
	TFile *_file0 = TFile::Open(filename);
	TH2D *h_inv_mass_US_Kminus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_US_Kminus_vs_pt"));
	h_inv_mass_US_Kminus_vs_pt->Sumw2();
	//h_inv_mass_US_Kminus_vs_pt->Rebin(rebin);
	TH2D *h_inv_mass_US_Kplus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_US_Kplus_vs_pt"));
	h_inv_mass_US_Kplus_vs_pt->Sumw2();
	//h_inv_mass_US_Kplus_vs_pt->Rebin(rebin);
	TH2D *h_inv_mass_MIX_Kminus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_LS_Kminus_vs_pt"));
	h_inv_mass_MIX_Kminus_vs_pt->Sumw2();
	//h_inv_mass_MIX_Kminus_vs_pt->Rebin(rebin);
	TH2D *h_inv_mass_MIX_Kplus_vs_pt = (TH2D*)(_file0->Get("h_inv_mass_LS_Kplus_vs_pt"));
	h_inv_mass_MIX_Kplus_vs_pt->Sumw2();
	//h_inv_mass_MIX_Kplus_vs_pt->Rebin(rebin);
	TString output_file = "pics/"+filename;
	//output_file.ReplaceAll(".root","_d_meson_mix_pt_1700.pdf");
	output_file.ReplaceAll(".root","_d_meson_ls_pt.pdf");
	TString output_file_txt = filename;
	//output_file_txt.ReplaceAll(".root",".txt");
	time_t  timev;
	time(&timev);
	std::stringstream ss;
	ss << timev;
	//TString sts = "_" + ss.str() + ".txt";
	TString s_array = Form("%d_%d_%d_%d", timev, bin_array[0], bin_array[1], bin_array[2]);
	TString sts = "_" + s_array + ".txt";
	//output_file_txt = "results_txt/" + output_file_txt + ss.str();
	output_file_txt.ReplaceAll(".root", sts);
	output_file_txt = "results_txt/" + output_file_txt;
	//cout << output_file_txt << endl;
	ofstream results_file (output_file_txt,  ios::trunc);
	results_file.close();
	cout << "output_file_txt: " << output_file_txt << endl;
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	//c2->Divide(2,2);
	c2->Print(output_file+"[");
	for (int i_pt = 0; i_pt < n_pt; i_pt++) {
		c2->cd();
		//cout << "bin_array: " << bin_array[i_pt] << " " << bin_array[i_pt+1] << endl;
                TH1D *proj_us_Kminus = h_inv_mass_US_Kminus_vs_pt->ProjectionX("proj_us_Kminus", bin_array[i_pt], bin_array[i_pt+1]);
                TH1D *proj_us_Kplus = h_inv_mass_US_Kplus_vs_pt->ProjectionX("proj_us_Kplus", bin_array[i_pt], bin_array[i_pt+1]);
		TH1D *proj_mix_Kminus = h_inv_mass_MIX_Kminus_vs_pt->ProjectionX("proj_mix_Kminus", bin_array[i_pt], bin_array[i_pt+1]);
                TH1D *proj_mix_Kplus = h_inv_mass_MIX_Kplus_vs_pt->ProjectionX("proj_mix_Kplus", bin_array[i_pt], bin_array[i_pt+1]);
		/*
		//c2->cd(1);
		//proj_us_Kminus->DrawCopy();
		//proj_mix_Kminus->DrawCopy("sames");
		//proj_us_Kminus->Draw();
		//proj_mix_Kminus->Draw("sames");
		proj_mix_Kminus->Draw();
		proj_mix_Kminus->SetLineColor(kBlue);
		proj_us_Kminus->Draw("sames");
		gPad->Update();
		TPaveStats* sb=(TPaveStats*)(proj_mix_Kminus->GetListOfFunctions()->FindObject("stats"));
		if (sb != 0) {
		sb->SetY1NDC(0.6);
		sb->SetY2NDC(0.75);
		sb->SetLineColor(kBlue);
		sb->Draw();
		gPad->Modified();
		}

		//c2->cd(2);
		proj_mix_Kplus->SetLineColor(kBlue);
		proj_mix_Kplus->Draw();
		proj_us_Kplus->Draw("sames");
		gPad->Update();
		TPaveStats* sb2=(TPaveStats*)(proj_mix_Kplus->GetListOfFunctions()->FindObject("stats"));
		if (sb2 != 0) {
		sb2->SetY1NDC(0.6);
		sb2->SetY2NDC(0.75);
		sb2->SetLineColor(kBlue);
		sb2->Draw();
		gPad->Modified();
		}
		*/

		//c2->cd(3);
		TH1D *proj_us = proj_us_Kminus->Clone();
		proj_us->Add(proj_us_Kplus, 1);
                
                TH1D *proj_mix = proj_mix_Kminus->Clone();
		proj_mix->Add(proj_mix_Kplus, 1);
                
                int int_down = 1700;
		int int_up = 1800;
		
                double int_us = proj_us -> Integral(int_down, int_up);
		double int_mix = proj_mix -> Integral(int_down, int_up);
		cout << "integral us: " << int_us << " mix: " << int_mix << endl;
		proj_mix->Scale(int_us/int_mix);
		proj_mix->SetLineColor(kBlue);
// 		proj_mix->Draw("sames");
		
                TPaveStats* sb3=(TPaveStats*)(proj_mix->GetListOfFunctions()->FindObject("stats"));
		if (sb3 != 0) {
			sb3->SetY1NDC(0.6);
			sb3->SetY2NDC(0.75);
			sb3->SetLineColor(kBlue);
			sb3->Draw();
			gPad->Modified();
		}

		//c2->cd(4);
		TH1D *signal = proj_us->Clone();
		signal->Add(proj_mix, -1);
		signal->Rebin(rebin);
		//signal->GetXaxis()->SetRangeUser(1.7,2.0);
		signal->Draw();
		c2->Print(output_file);

		//TCanvas *c4 = draw(signal, rebin);
// 		float min_pt = bin_array[i_pt]*rebin/1000.;
                float min_pt = bin_array[i_pt]/100.;
                
// 		float max_pt = bin_array[i_pt+1]*rebin/1000.;
		float max_pt = bin_array[i_pt+1]/100.;
		                
		TCanvas *c4 = draw(signal, rebin, min_pt, max_pt, output_file_txt, inc_stat, Ntrig);
		c4->cd();
		tx1.DrawLatex(0.1,0.94,Form("p_{T}: %3.1f-%3.1f [GeV/c]", min_pt, max_pt));
		c4->Print(output_file);
		TString output_file2 = output_file;
		output_file2.ReplaceAll(".pdf",Form("_%d.png", i_pt));
		c4->Print(output_file2);
	}

	c2->Print(output_file+"]");
	cout << "output_file_txt: " << output_file_txt << endl;
        
}

TCanvas *draw(TH1D *signal, int rebin, double low_pt = 1.0, double high_pt = 2.0, TString output_file_txt = "def.txt", float inc_stat = 1, double Ntrig = 41E6) {
	//const float fitRMin   = 1.73;
	const float fitRMin   = 1.72;
	//const float fitRMin   = 1.7;
	const float fitRMax   = 2.05;
	const float MKSize    = 1.;
	const float rotwthmin = 1.85;
	//const float rotwthmin = 1.85;
	const float rotwthmax = 1.89;
	TLatex tx1;


	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	TF1 *fun0 = new TF1("fun0","pol1(0)+gaus(2)",fitRMin,fitRMax);
	fun0->SetParameters(1.,1.,1.,1.865,0.015);
	fun0->SetLineColor(4);
	fun0->SetLineStyle(7);

	const int N = 400/rebin;
	const int binm1 = 1700/rebin;
	const int binm2 = 185;
	const int binm3 = 192;
	const int binm4 = 2100/rebin;
	const double scalem = 1.;
	double mm[N],ym[N],yme[N],ym1[N];

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
        resfunm->SetLineColor(4);
	//resfunm->Draw("same");
	//fun0->Draw("same");

	TF1 *resfunm1 = new TF1("resfunm1","pol1(0)+gaus(2)",fitRMin,fitRMax);
	resfunm1->SetParameters(0.,0.,fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
	resfunm1->SetLineColor(2);
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

	tx1.DrawLatex(0.17,0.84,"SE-ME");

	float chi2m   = fun0->GetChisquare();
	float ndfm    = fun0->GetNDF();
	float meanm   = fun0->GetParameter(3);
	float meanme  = fun0->GetParError(3);
	float sigm    = fun0->GetParameter(4);
	float sigme   = fun0->GetParError(4);

	float width_bin = 0.001 * rebin;
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
	/*
	int yieldmBinMin = 185-binm1;
	int yieldmBinMax = 190-binm1;
	*/
	int yieldmBinMin = 1850/rebin - binm1;
	int yieldmBinMax = 1900/rebin - binm1;
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

	double Eff = 0.05;
	double pt = (high_pt + low_pt)/2.;
	double pt_bin_width = high_pt - low_pt;
	double inv_yield = yieldm / (4*TMath::Pi() * Ntrig * 0.0389 * pt_bin_width * pt * 2 * Eff);
	cout << "inv_yield = " << yieldm << " / (4 * TMath::Pi() * " << Ntrig << " * 0.0389 * " << pt_bin_width << " * " << pt << " * 2 * "<< Eff << ") = " << inv_yield << endl;
	double inv_yield_err = yieldme / (4*TMath::Pi() * Ntrig * 0.0389 * pt_bin_width * pt * 2 * Eff);
	cout << "Invariant Yield: " << low_pt << " " << high_pt << " " << inv_yield << " " << inv_yield_err << endl;
	//results_file << low_pt << " " << high_pt << " " << inv_yield << " " << inv_yield_err << endl;
	results_file << low_pt << "\t" << high_pt << "\t" << inv_yield << "\t" << inv_yield_err << endl;

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

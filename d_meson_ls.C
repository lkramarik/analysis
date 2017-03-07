void d_meson_ls(TString filename = "results_all.root", int rebin = 10) {

	const float txtsize   = 0.055;
	const float titsize   = 0.055;
	const float LFMargin  = 0.16;
	const float yToffset  = 1.4;
	const float xmin      = 1.72;
	const float xmax      = 2.07;
	const float y2min     = -1.4;
	const float y2max     = 7.6;
	const float MKSize    = 1.;
	const float fitRMin   = 1.73;
	const float fitRMax   = 2.05;
	const float rotwthmin = 1.857;
	const float rotwthmax = 1.89;
	const float lkswthmin = 1.859;
	const float lkswthmax = 1.89;
	const int   pTBinMin  = 3;
	const int   pTBinMax  = 22;

	//int rebin = 4;
	//int rebin = 10;
	//int rebin = 20;
	TFile *_file0 = TFile::Open(filename);
	TH1D *h_inv_mass_LS = (TH1D*)(_file0->Get("h_inv_mass_LS"));
	//h_inv_mass_LS->Sumw2();
	//h_inv_mass_LS->Rebin(rebin);
	TH1D *h_inv_mass_LS_neg = (TH1D*)(_file0->Get("h_inv_mass_LS_neg"));
	//h_inv_mass_LS_neg->Sumw2();
	//h_inv_mass_LS_neg->Rebin(rebin);
	TH1D *h_inv_mass_LS_pos = (TH1D*)(_file0->Get("h_inv_mass_LS_pos"));
	//h_inv_mass_LS_pos->Sumw2();
	//h_inv_mass_LS_pos->Rebin(rebin);
	TH1D *h_inv_mass_US = (TH1D*)(_file0->Get("h_inv_mass_US"));
	//double int_us = h_inv_mass_US->Integral(1100,1550) ;
	//double int_us = h_inv_mass_US->Integral(1200,1700) ;
	//double int_us = h_inv_mass_US->Integral(1500,1700) ;
	double int_us = h_inv_mass_US->Integral(1700,1800) ;
	//double int_us = 1 ;
	//h_inv_mass_US->Sumw2();
	//h_inv_mass_US->Rebin(rebin);
	//h_inv_mass_US_Kminus_vs_pt->ProjectionX("p1",0,130)
	//TH1D *h_inv_mass_LS_norm = new TH1D("h_inv_mass_LS_norm", "h_inv_mass_LS_norm", 3000, 0., 3.0);
	TH1D *h_inv_mass_LS_norm = new TH1D("h_inv_mass_LS_norm", "h_inv_mass_LS_norm", 3000, 0., 3.0);
	h_inv_mass_LS_norm->Sumw2();
	for (int i = 0; i < 3000; i++) {
		double norm = 2 * TMath::Sqrt(h_inv_mass_LS_neg->GetBinContent(i) * h_inv_mass_LS_pos->GetBinContent(i));
		h_inv_mass_LS_norm->SetBinContent(i, norm);
		/*
		double err;
		if (h_inv_mass_LS_pos->GetBinContent(i) > 0 && h_inv_mass_LS_neg->GetBinContent(i) > 0) {
			err = TMath::Sqrt( h_inv_mass_LS_neg->GetBinContent(i) * TMath::Power( h_inv_mass_LS_pos->GetBinError(i) / h_inv_mass_LS_pos->GetBinContent(i), 2) + 
				           h_inv_mass_LS_pos->GetBinContent(i) * TMath::Power( h_inv_mass_LS_neg->GetBinError(i) / h_inv_mass_LS_neg->GetBinContent(i), 2) );
		} else {
			err = 0;
		}
		*/
		double y1e = h_inv_mass_LS_pos->GetBinError(i);
		double y2e = h_inv_mass_LS_neg->GetBinError(i);
		//double err = 0.5*sqrt(y1e*y1e + y2e*y2e);
		double err = sqrt(y1e*y1e + y2e*y2e);
		h_inv_mass_LS_norm->SetBinError(i, err);
	}
	//double int_ls = h_inv_mass_LS_norm->Integral(1100,1550);
	//double int_ls = h_inv_mass_LS_norm->Integral(1200,1700);
	//double int_ls = h_inv_mass_LS_norm->Integral(1500,1700);
	double int_ls = h_inv_mass_LS_norm->Integral(1700,1800);
	//double int_ls = 1;
	h_inv_mass_LS_norm->Rebin(rebin);
	h_inv_mass_LS_neg->Rebin(rebin);
	h_inv_mass_LS_pos->Rebin(rebin);
	h_inv_mass_LS->Rebin(rebin);




	/*
	TFile *_file0 = TFile::Open(filename);
	TH1D *h_inv_mass_mix = (TH1D*)(_file0->Get("h_inv_mass_MIX"));
	h_inv_mass_mix->Sumw2();
	//double int_ls = h_inv_mass_mix->Integral(1100,1250);
	double int_ls = h_inv_mass_mix->Integral(1100,1550);
	h_inv_mass_mix->Rebin(rebin);
	TH1D *h_inv_mass_US = (TH1D*)(_file0->Get("h_inv_mass_US"));
	h_inv_mass_US->Sumw2();
	//double int_us = h_inv_mass_US->Integral(1100,1250) ;
	double int_us = h_inv_mass_US->Integral(1100,1550) ;
	*/
	h_inv_mass_US->Rebin(rebin);
	h_inv_mass_US->SetTitle("");
	//h_inv_mass_US->GetXaxis()->SetRangeUser(0.5,2.2);
	h_inv_mass_US->GetXaxis()->SetRangeUser(0.6,2.2);
	h_inv_mass_US->GetYaxis()->SetTitle("Raw Counts");
	h_inv_mass_US->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
	//TCanvas *c2 = new TCanvas("c2","c2",1800,700);

	cout << "bin US LS_norm LS " << endl;
	for(int ib=180; ib<200; ib++) {
		//cout << ib << " " << h_inv_mass_US->GetBinError(ib) << " " << h_inv_mass_LS_norm->GetBinError(ib) << " " << h_inv_mass_LS_pos->GetBinError(ib) << " " << h_inv_mass_LS_neg->GetBinError(ib) << endl;
		cout << ib << " " << h_inv_mass_US->GetBinError(ib) << " " << h_inv_mass_LS_norm->GetBinError(ib) << " " << h_inv_mass_LS->GetBinError(ib) << endl;
	}

	/*
	for(int ib=40; ib<200; ib++) {
		cout << ib << " " << h_inv_mass_US->GetBinError(ib) << " " << h_inv_mass_US->GetBinContent(i) << h_inv_mass_US->GetBinError(ib)/ h_inv_mass_US->GetBinContent(i) << endl;
	}
	*/




	TCanvas *c2 = new TCanvas("c2","c2",1200,500);
	//c2->Divide(2,2);
	c2->Divide(3,1);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	//gStyle->SetEndErrorSize(0.01);
	/*
	c2->SetFillColor(10);
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameFillColor(0);
	c2->SetFrameBorderMode(0);
	c2->SetBorderSize(2);
	//  c2->SetLogy();
	//  c2->SetGridx();
	//  c2->SetGridy();
	c2->SetLeftMargin(0.12);
	c2->SetBottomMargin(0.15);
	c2->SetTopMargin(0.025);
	c2->SetRightMargin(0.025);
	*/

	TLatex tx1;
	tx1.SetNDC();


	c2->cd(1);
	//h_inv_mass_mix->Scale(h_inv_mass_US->Integral(1100,1250)/h_inv_mass_mix->Integral(1100,1250));
	cout << "integral us: " << int_us << " ls: " << int_ls << endl;
	//h_inv_mass_LS_norm->Scale(int_us/int_ls);
	h_inv_mass_LS_norm->SetLineColor(kBlue);
	h_inv_mass_LS_norm->SetTitle("");
	h_inv_mass_LS_norm->GetXaxis()->SetRangeUser(0.5,2.2);
	//h_inv_mass_LS_norm->GetYaxis()->SetNdivisions(208);
	h_inv_mass_LS_norm->GetYaxis()->SetTitleOffset(1.5);
	h_inv_mass_LS_norm->GetYaxis()->SetTitle("Raw Counts");
	h_inv_mass_LS_norm->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
	h_inv_mass_LS_norm->SetLineColor(kRed);
	h_inv_mass_LS_norm->DrawCopy();
	//h_inv_mass_LS_norm->Draw("sames");
	h_inv_mass_US->DrawCopy("sames");
	//h_inv_mass_LS_norm->DrawCopy("sames");
	
	TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetShadowColor(0);
	leg->SetTextFont(22);
	leg->SetTextSize(0.04);
	//leg->SetHeader("Unlike sign K#pi pairs");
	leg->SetHeader("K#pi pairs");
	leg->AddEntry(h_inv_mass_US,"Unlike sign (US)","P");
	//leg->AddEntry(h_inv_mass_US,"Same event");
	leg->AddEntry(h_inv_mass_LS_norm,"Like sign (LS)","P");
	leg->Draw();


	/*
	c2->cd(2);
	h_inv_mass_LS_norm->SetLineColor(kBlue);
	h_inv_mass_LS_norm->GetXaxis()->SetRangeUser(1.7,2.0);
	h_inv_mass_LS_norm->DrawCopy();
	//h_inv_mass_LS_norm->Draw("sames");
	h_inv_mass_US->GetXaxis()->SetRangeUser(1.7,2.0);
	h_inv_mass_US->DrawCopy("sames");
	h_inv_mass_LS_norm->SetLineColor(kRed);
	h_inv_mass_LS_norm->GetXaxis()->SetRangeUser(1.7,2.0);
	h_inv_mass_LS_norm->DrawCopy("sames");
	*/

	//c2->cd(3);
	c2->cd(2);
	h_inv_mass_US->Add(h_inv_mass_LS_norm, -1);
	h_inv_mass_US->GetYaxis()->SetTitleOffset(1.5);
	h_inv_mass_US->DrawCopy();
	tx1.DrawLatex(0.37,0.84,"US-LS");

	/*
	c2->cd(3);
	h_inv_mass_US->SetMarkerStyle(kOpenCircle);
	h_inv_mass_US->SetMarkerColor(kBlack);
	h_inv_mass_US->Draw("pe");
	h_inv_mass_US->GetXaxis()->SetRangeUser(1.7,2.1);
	//h_inv_mass_US->SetStats(kFALSE);
	//h_inv_mass_US->Fit(fun0,"NOR");
	   TF1 *resfunm = new TF1("resfunm","pol2",fitRMin,fitRMax);
	   resfunm->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1),fun0->GetParameter(2));
	   resfunm->SetLineStyle(7);
	   resfunm->SetLineWidth(1);
	   resfunm->Draw("same");
	   fun0->Draw("same");
	   */

	c2->cd(3);

	/*
	TH1D *htmp1 = new TH1D("htmp1","",1,xmin,xmax);
	htmp1->GetYaxis()->SetTitle("Raw Counts (#times 10^{3})");
	htmp1->GetYaxis()->SetTitleSize(titsize);
	htmp1->GetYaxis()->SetTitleOffset(yToffset);
	htmp1->GetYaxis()->SetLabelOffset(0.03);
	htmp1->GetYaxis()->SetLabelSize(0.047);
	htmp1->SetMinimum(y2min);
	htmp1->SetMaximum(y2max);
	htmp1->GetXaxis()->SetNdivisions(208);
	//  htmp1->GetXaxis()->CenterTitle();
	htmp1->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
	htmp1->GetXaxis()->SetTitleOffset(1.1);
	htmp1->GetXaxis()->SetTitleSize(titsize);
	htmp1->GetXaxis()->SetLabelOffset(0.01);
	htmp1->GetXaxis()->SetLabelSize(0.047);
	htmp1->GetXaxis()->SetLabelFont(42);
	htmp1->GetXaxis()->SetTitleFont(42);
	htmp1->GetYaxis()->SetNdivisions(208);
	htmp1->Draw();
	*/

	TF1 *fun0 = new TF1("fun0","pol1(0)+gaus(2)",fitRMin,fitRMax);
	fun0->SetParameters(1.,1.,1.,1.865,0.015);
	fun0->SetLineColor(2);
	fun0->SetLineStyle(7);

	//const int N = 22;
	//const int N = 35;
	//const int N = 88;
	const int N = 400/rebin;
	const int binm1 = 1700/rebin;
	//const int binm1 = 432;
	//const int binm1 = 173;
	//const int binm1 = 1730;
	const int binm2 = 185;
	const int binm3 = 192;
	const int binm4 = 2100/rebin;
	//const int binm4 = 520;
	//const int binm4 = 208;
	//const int binm4 = 2080;
	const double scalem = 1.;
	double mm[N],ym[N],yme[N],ym1[N];

	cout << h_inv_mass_US-> GetNbinsX() << endl;
	int i=0;
	//for(int ib=binm2; ib<binm3; ib++) {
	for(int ib=binm1; ib<binm4; ib++) {
		mm[i]  = h_inv_mass_US->GetBinCenter(ib);
		ym[i]  = h_inv_mass_US->GetBinContent(ib) + scalem;
		yme[i] = h_inv_mass_US->GetBinError(ib);
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
	fun0->SetParLimits(4,0.010,0.014);
	//  h_inv_mass_US->Fit(fun0,"NOR");
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
	//mg->GetXaxis()->SetTitleOffset(1.5);
	//mg->GetXaxis()->SetTitleSize(titsize);
	//mg->GetXaxis()->SetLabelOffset(0.01);
	//mg->GetXaxis()->SetLabelSize(0.047);
	//mg->GetXaxis()->SetLabelFont(42);
	//mg->GetXaxis()->SetTitleFont(42);
	//mg->GetYaxis()->SetNdivisions(208);

	//mg->Draw("psame");
	resfunm->Draw("same");
	fun0->Draw("same");
	resfunm1->Draw("same");

	tx1.DrawLatex(0.17,0.84,"US-LS");

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
	cout << " (fun0->GetParameter(4)-5*sigm): " << (fun0->GetParameter(4)-5*sigm) << " 5*sigm+fun0->GetParameter(4)): " << 5*sigm+fun0->GetParameter(4) << endl;
	float nsig    = yieldm / TMath::Sqrt(yieldm + bckg);
	cout << "sigm: " << sigm << " fun0->GetParameter(3): " << fun0->GetParameter(3) << endl;
	cout << " fun0->GetParameter(1): " << fun0->GetParameter(1) << endl;
	cout << " fun0->GetParameter(2): " << fun0->GetParameter(2) << endl;
	cout << " fun0->GetParameter(3): " << fun0->GetParameter(3) << endl;
	cout << " fun0->GetParameter(4): " << fun0->GetParameter(4) << endl;
	cout << " fun0->GetParameter(5): " << fun0->GetParameter(5) << endl;
	cout << " fun0->GetParameter(6): " << fun0->GetParameter(6) << endl;
	cout << "bckg: " << bckg << " sig: " << yieldm << " N_sig: " << nsig << endl;

	char chh1[250],chh2[250],chh3[250],chh4[250], chh5[250], chh6[250];

	float yieldmBin = 0;
	float yieldmeBin = 0;
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
	sprintf(chh6,"yieldBin = %7.0f #pm %7.0f",yieldmBin,yieldmeBin);
	sprintf(chh3,"mean = %5.3f #pm %5.3f",meanm,meanme);
	sprintf(chh4,"#sigma = %5.3f #pm %5.3f",sigm,sigme);
	sprintf(chh5,"N_sig = %2.1f", nsig);
	//sprintf(chh4,"N_sig = %2.1f #pm %2.1f", nsig, nsige);

	//TPaveStats *pave1 = new TPaveStats(0.58,0.73,0.98,0.93,"brNDC");
	TPaveStats *pave1 = new TPaveStats(0.48,0.78,0.98,0.95,"brNDC");
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

	
	c2->Update();



	//TString output_file = filename;
	TString output_file = "pics/"+filename;
	output_file.ReplaceAll(".root","_d_meson_ls.pdf");
	c2->Print(output_file);
	output_file.ReplaceAll(".pdf",".png");
	c2->Print(output_file);
}

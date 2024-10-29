void testReferenceGammaPeak()
{
	//TFile f("/home/terawatt/server_/ZFSRAID/Ing27-HPGe/root/test147_HPGe_2_Labr_4_Ing27_87prc_Generator_100hz_Ba_center_Pb_1mm_Box_2cm.root");
	TFile f("/home/terawatt/server_/ZFSRAID/Ing27-HPGe/root/test161_HPGe_2_Labr_4_Ing27_84prc_Generator_100hz_Ba_center_Pb_1mm_CH2Cl2_4cm.root");
	TH2F *h_sp=(TH2F*)f.Get("h_3_0");
	
	h_sp->RebinX(5);
	h_sp->RebinY(2);
	
	ATOFProcess p;
	p.FullSpectrum=*h_sp;
	p.GenerateAntiCoincedence(-450,-50,-20,50);
	
	cout<<h_sp->GetYaxis()->GetBinCenter(1)<<"\n";
	cout<<p.PureCoincedence.GetYaxis()->GetBinCenter(1)<<"\n";
	ReferenceGammaPeak RGP;
	RGP.XMin=4100;
	RGP.XMax=4800;
	RGP.PeakMin=4400;
	RGP.PeakMax=4500;
	RGP.Energy=4450;
	RGP.fProcess=&p;
	RGP.Use2Dhist=true;
	RGP.GenerateSubstrateHistogram();
	TFile f_out("out.root","recreate");
	f_out.WriteTObject(h_sp);
	f_out.WriteTObject(&RGP.FullHist);
	f_out.WriteTObject(&RGP.SubstrateHist);
	f_out.WriteTObject(&RGP.PeakHist);
	f_out.WriteTObject(&RGP.FullX);
	f_out.WriteTObject(&RGP.FullY);
	f_out.WriteTObject(&RGP.SubstrateX);
	f_out.WriteTObject(&RGP.SubstrateY);
	f_out.WriteTObject(&RGP.PeakX);
	f_out.WriteTObject(&RGP.PeakY);
	f_out.Close();
}

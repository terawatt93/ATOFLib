#include "ATOFLib.hh"

void Transpose()
{
	TFile f("Spectra2D.root");
	TH2F *h=(TH2F*)f.Get("h_8_5");
	h->Draw();
	TH2F h2=TransposeTH2(h);
	TFile f2("f2.root","recreate");
	f2.WriteTObject(&h2);
	f2.Close();
}

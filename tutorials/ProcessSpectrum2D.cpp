#include "ATOFLib.hh"

void FindPeakFromSample(vector<TF1> &Components,double Pos, double &PeakPos, double &Err, int &PeakNumber)
{
	double Penalty=1e9;
	for(unsigned int i=0;i<Components.size();i++)
	{
		cout<<"pos: "<<Pos<<" "<<Pos-Components[i].GetParameter(1)<<" "<<Components[i].GetParameter(2)<<" "<<abs(Pos-Components[i].GetParameter(1))+pow(Components[i].GetParameter(2),2)<<"\n";
		if(Penalty>(abs(Pos-Components[i].GetParameter(1))+pow(Components[i].GetParameter(2),2)))
		{
			Penalty=abs(Pos-Components[i].GetParameter(1))+pow(Components[i].GetParameter(2),2);
			PeakPos=Components[i].GetParameter(1);
			Err=Components[i].GetParError(1);
			PeakNumber=i;
		}
	}
}


void EstimateBorders(vector<TF1> &Components,TF1 &Function,double &Left, double &Right,double Fraction=0.96)
{
	Left=1e9; Right=-1e9;
	for(unsigned int i=0;i<Components.size();i++)
	{
		if((Left>Components[i].GetParameter(1)-3*Components[i].GetParameter(2)))
		{
			Left=Components[i].GetParameter(1)-3*Components[i].GetParameter(2);
		}
		if((Right<Components[i].GetParameter(1)+3*Components[i].GetParameter(2)))
		{
			Right=Components[i].GetParameter(1)+3*Components[i].GetParameter(2);
		}
	}
	double InitIntegral=Function.Integral(Left,Right);
	double Step=0.1;
	while(Function.Integral(Left,Right)/InitIntegral>Fraction)
	{
		double Left_r=Left+Step;
		double Right_r=Right-Step;
		double Integral1=Function.Integral(Left_r,Right);
		double Integral2=Function.Integral(Left,Right_r);
		if(Integral1<Integral2)
		{
			Right=Right-Step;
		}
		else
		{
			Left=Left+Step;
		}
	}
}

void ProcessHPGe(TOFWindow *w,double ReferencePosition)// у нас по-разному обрабатываются спектры для HPGe и LaBr
{
	//У HPGe одна компонента
	TOFComponent *comp=w->fATOF->GetOrCreateTOFComponent(0);
	comp->PredefAnalysisType="BordersGraph";
	int N=comp->LeftBorderGraph.GetN();
	//за окно мы принимаем 2сигма относительно самых удаленных друг от друга пиков
	double Left=1e9,Right=-1e9,LeftErr=0,RightErr=0;
	EstimateBorders(w->Components,w->FitFunction,Left,Right);
	
	/*for(unsigned int i=0;i<w->Components.size();i++)
	{
		if((Left>w->Components[i].GetParameter(1)-2*w->Components[i].GetParameter(2)) && (w->Components[i].GetParameter(2)<5))
		{
			Left=w->Components[i].GetParameter(1)-2*w->Components[i].GetParameter(2);
			LeftErr=sqrt(pow(w->Components[i].GetParError(1),2)+pow(2*w->Components[i].GetParError(2),2));
		}
		if((Right<w->Components[i].GetParameter(1)+2*w->Components[i].GetParameter(2)) && (w->Components[i].GetParameter(2)<5))
		{
			Right=w->Components[i].GetParameter(1)+2*w->Components[i].GetParameter(2);
			RightErr=sqrt(pow(w->Components[i].GetParError(1),2)+pow(2*w->Components[i].GetParError(2),2));
		}
	}*/
	if(N==0)
	{
		w->fATOF->TOFDependenceLeft=w->Centroid;
	}
	else
	{
		w->fATOF->TOFDependenceRight=w->Centroid;
	}
	cout<<"Left,Right="<<Left<<" "<<Right<<"\n";
	comp->LeftBorderGraph.SetPoint(N,w->Centroid,Left);
	comp->LeftBorderGraph.SetPointError(N,0,LeftErr);
	comp->RightBorderGraph.SetPoint(N,w->Centroid,Right);
	comp->RightBorderGraph.SetPointError(N,0,RightErr);
}

void ProcessLaBr(TOFWindow *w,double ReferencePosition)
{
	if(w->fATOF->TOFComponents.size()<3)
	{
		w->fATOF->TOFComponents.resize(3);
	}
	TOFComponent *Sample=w->fATOF->GetOrCreateTOFComponent(0);
	TOFComponent *Generator=w->fATOF->GetOrCreateTOFComponent(1);
	TOFComponent *Neutrons=w->fATOF->GetOrCreateTOFComponent(2);
	
	Sample->PredefAnalysisType="BordersGraph";
	Generator->PredefAnalysisType="BordersGraph";
	Neutrons->PredefAnalysisType="BordersGraph";
	
	vector<TOFComponent*> Components={Sample,Generator,Neutrons};
	double SamplePos,SamplePosError,GeneratorPos,GeneratorError,NeutronsPos,NeutronsErr;
	int SamplePeakNumber=-1,GeneratorPeakNumber=-1,NeutronPeakNumber=1;
	
	FindPeakFromSample(w->Components,ReferencePosition,SamplePos,SamplePosError,SamplePeakNumber);
	FindPeakFromSample(w->Components,ReferencePosition-6,GeneratorPos,GeneratorError,GeneratorPeakNumber);
	FindPeakFromSample(w->Components,ReferencePosition+8,NeutronsPos,NeutronsErr,NeutronPeakNumber);
	
	int N=Sample->LeftBorderGraph.GetN();
	cout<<"N="<<N<<"\n";
	cout<<"SamplePeakNumber:"<<SamplePeakNumber<<"\n";
	cout<<"LeftBorderGraph:"<<N<<" "<<w->Centroid<<" "<<SamplePos-2*w->Components[SamplePeakNumber].GetParameter(2)<<"\n";
	Sample->LeftBorderGraph.SetPoint(N,w->Centroid,SamplePos-2*w->Components[SamplePeakNumber].GetParameter(2));
	Sample->RightBorderGraph.SetPoint(N,w->Centroid,SamplePos+2*w->Components[SamplePeakNumber].GetParameter(2));
	
	Generator->LeftBorderGraph.SetPoint(N,w->Centroid,GeneratorPos-2*w->Components[GeneratorPeakNumber].GetParameter(2));
	Generator->RightBorderGraph.SetPoint(N,w->Centroid,GeneratorPos+w->Components[GeneratorPeakNumber].GetParameter(2));
	
	Neutrons->LeftBorderGraph.SetPoint(N,w->Centroid,SamplePos+2*w->Components[SamplePeakNumber].GetParameter(2));
	Neutrons->RightBorderGraph.SetPoint(N,w->Centroid,NeutronsPos+5*w->Components[NeutronPeakNumber].GetParameter(2));
}

void ProcessSpectrum2D()
{
	TH1::AddDirectory(false);
	TF1 f_HPGe("f_HPGe","gaus(0)+gaus(3)+gaus(6)+gaus(9)+pol1(12)",5.518050,81.534700);
	TF1 f_LaBr("f_HPGe","gaus(0)+gaus(3)+gaus(6)+gaus(9)+pol1(12)",33.956700,104.399000);
	vector<double> ParHPGe={3392.527507,41.790317,3.627782,3946.630890,43.808792,1.852620,1466.782824,50.321583,1.887797,1542.721396,47.341043,9.096406,2318.567243,-0.083779};
	vector<double> ParLaBr={11008.018270,59.432300,1.460086,11520.801565,65.955535,1.725138,4405.836594,72.620964,2.144907,1352.663122,79.522543,2.391178,15520.623933,1.888519};
	
	vector<vector<double> > TWindows={{100,200},{150,300},{250,500},{400,800},{700,1000},{800,1500},{1200,3000},{2000,4000},{3000,5000},{4000,7000},{6000,10000}};
	for(int i=0;i<f_HPGe.GetNpar();i++)
	{
		f_HPGe.SetParameter(i,ParHPGe[i]);
		f_LaBr.SetParameter(i,ParLaBr[i]);
	}
	
	
	TFile f("Spectra2D.root");
	ATOFProcess p_HPGe, p_LaBr;
	
	p_HPGe.GenerateTH2=true;
	p_LaBr.GenerateTH2=true;
	
	p_HPGe.ReadFromTFile_Full(&f,"h_0_0");
	p_LaBr.ReadFromTFile_Full(&f,"h_8_5");
	
	
	p_HPGe.GenerateAntiCoincedence(-300,0,0,100);
	p_LaBr.GenerateAntiCoincedence(-300,0,0,100);
	
	p_HPGe.GenerateTOFWindows(&TWindows);
	p_LaBr.GenerateTOFWindows(&TWindows);
	
	p_HPGe.AttachFitFunction(&f_HPGe);
	p_LaBr.AttachFitFunction(&f_LaBr);
	
	for(unsigned int i=0;i<TWindows.size();i++)
	{
		p_HPGe.TOFWindows[i].FitWindow();
		p_LaBr.TOFWindows[i].FitWindow();
	}
	
	p_HPGe.SelectPeaksFCN=ProcessHPGe;
	p_LaBr.SelectPeaksFCN=ProcessLaBr;
	
	p_HPGe.GenerateComponents(43.8);
	//p_HPGe.ProcessComponents();
	p_LaBr.GenerateComponents(65.95);
	
	/*p_HPGe.ProcessComponents();
	p_LaBr.ProcessComponents();
	*/
	TCanvas c;
	
	c.Print("HPGe.pdf[","pdf");
	p_HPGe.PureCoincedence.ProjectionX()->Draw("hist");
	c.Print("HPGe.pdf","pdf");
	
	p_HPGe.PureCoincedence.ProjectionY()->Draw("hist");
	c.Print("HPGe.pdf","pdf");

	for(unsigned int i=0;i<p_HPGe.TOFWindows.size();i++)
	{
		p_HPGe.TOFWindows[i].Draw();
		p_HPGe.TOFWindows[i].FitFunction.Draw("same");
		c.Print("HPGe.pdf","pdf");
	}
	
	for(unsigned int i=0;i<p_HPGe.TOFComponents.size();i++)
	{
		p_HPGe.TOFComponents[i].SpectrumHist2D.RebinX(10);
		p_HPGe.TOFComponents[i].SpectrumHist2D.RebinY(2);
		
		p_HPGe.TOFComponents[i].SpectrumHist2D.Draw("colz");
		p_HPGe.TOFComponents[i].LeftBorderGraph.Draw("lp");
		p_HPGe.TOFComponents[i].RightBorderGraph.Draw("lp");
		c.Print("HPGe.pdf","pdf");
		
		p_HPGe.TOFComponents[i].SpectrumHist2D.ProjectionY()->Draw();
		c.Print("HPGe.pdf","pdf");
		p_HPGe.TOFComponents[i].SpectrumHist.Draw("hist");
		c.Print("HPGe.pdf","pdf");
		
	}
	c.Print("HPGe.pdf]","pdf");
	
	c.Print("LaBr.pdf[","pdf");
	p_LaBr.PureCoincedence.ProjectionX()->Draw("hist");
	c.Print("LaBr.pdf","pdf");
	
	p_LaBr.PureCoincedence.ProjectionY()->Draw("hist");
	c.Print("LaBr.pdf","pdf");

	for(unsigned int i=0;i<p_LaBr.TOFWindows.size();i++)
	{
		p_LaBr.TOFWindows[i].Draw();
		c.Print("LaBr.pdf","pdf");
	}
	
	for(unsigned int i=0;i<p_LaBr.TOFComponents.size();i++)
	{
		p_LaBr.TOFComponents[i].SpectrumHist2D.Draw("colz");
		p_LaBr.TOFComponents[i].LeftBorderGraph.Draw("lp");
		p_LaBr.TOFComponents[i].RightBorderGraph.Draw("lp");
		c.Print("LaBr.pdf","pdf");
	}
	p_LaBr.TOFComponents[0].SpectrumHist.Draw("hist");
	for(unsigned int i=1;i<p_LaBr.TOFComponents.size();i++)
	{
		p_LaBr.TOFComponents[i].SpectrumHist.SetLineColor(i+1);
		p_LaBr.TOFComponents[i].SpectrumHist.Draw("hist same");
		
	}
	c.Print("LaBr.pdf","pdf");
	c.Print("LaBr.pdf]","pdf");
}

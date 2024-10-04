
#include "ATOFLib.hh"
#include <future>
#include <TROOT.h>   
#include <TSpectrum.h>
#pragma once

void LinearRegression(vector<Float_t> *xData, vector<Float_t> *sData,int indexMin, int indexMax,double &a, double &b)
{
	//cout<<"indexMin indexMax:"<<indexMin<<" "<<indexMax<<"\n";
	double Sx=0,Sy=0,Sxy=0,Sxx=0;
	a=0;
	b=0;
	if(indexMin<0)
	{
		indexMin=0;
	}
	if(indexMax>=sData->size()-1)
	{
		indexMax=sData->size()-1;
	}
	int n=0;
	//cout<<"indexMin indexMax:"<<indexMin<<" "<<indexMax<<"\n";
	if((indexMin<sData->size())&&(indexMax<sData->size())&&(indexMin<indexMax))
	{
		for(int i=indexMin; i < indexMax; i++)
		{
		 Sx+=xData->at(i);//x
		 Sy+=sData->at(i);
		 Sxy+=xData->at(i)*sData->at(i);
		 Sxx+=xData->at(i)*xData->at(i);
		 n++;
		}
	}
	
	Sx/=n;
	Sy/=n;
	Sxy/=n;
	Sxx/=n;
	a=(Sx*Sy-Sxy)/(Sx*Sx-Sxx);
	b=(Sxy-a*Sxx)/Sx;
}

int Type(int det)
{
	if(det<2)
	{
		return 1;
	}
	if(det<6)
	{
		return 2;
	}
	return 0;
}
void AddPointToTGraphErrors(TGraphErrors *gr,double x,double y,double ex, double ey)
{
	if(gr->GetN()==0)
	{
		gr->SetPoint(0,x,y);
		gr->SetPointError(0,ex,ey);
	}
	else
	{
		double xMax,yMax;
		gr->GetPoint(gr->GetN()-1,xMax,yMax);
		cout<<"xMax:"<<xMax<<" "<<x<<"\n";
		if(xMax<x)
		{
			int N=gr->GetN();
			gr->SetPoint(N,x,y);
			gr->SetPointError(N,ex,ey);
		}
		else
		{
			bool Found=false;
			for(unsigned int i=0;i<gr->GetN()-1;i++)
			{
				double _x,_y;
				gr->GetPoint(i,_x,_y);
				if(abs(_x-x)<0.5)//замена точки
				{
					gr->SetPoint(i,x,y);
					gr->SetPointError(i,ex,ey);
					Found=true;
				}
			}
			if(!Found)
			{
				for(unsigned int i=0;i<gr->GetN()-1;i++)
				{
					double _x,_y;
					double __x,__y;
					gr->GetPoint(i,_x,_y);
					gr->GetPoint(i+1,__x,__y);
					if((x>_x)&&(x<__x-0.1))
					{
						gr->InsertPointBefore(i+1,x,y);
						gr->SetPointError(i+1,ex,ey);
					}
				}
			}
		}
		cout<<"x,y:"<<x<<" "<<y<<"\n";
		
	}
	
}

ATOFProcess *FindBestMatch(ATOFProcess* current, vector<ATOFProcess> &ATOFP, bool Fitted=false)
{
	
	if(!current->FittedTOF)
	{
		//сначала ищем с ручным фитом и тем же детектором
		//сначала пытаемся найти с тем же детектором
		bool HasManualFits=false;
		for(unsigned int i=0;i<ATOFP.size();i++)
		{
			if((current->DetNumber==ATOFP[i].DetNumber)&&(ATOFP[i].GeneratedAnti))
			{
				if((!Fitted))
				{
					cout<<"current: "<<current->FullSpectrum.GetName()<<" Found: "<<ATOFP[i].FullSpectrum.GetName()<<"\n"; return &ATOFP[i];
				}
				else if((ATOFP[i].FittedTOF)&&(ATOFP[i].HasManualFits()>0))
				{
					HasManualFits=true;
					cout<<"current: "<<current->FullSpectrum.GetName()<<" Found: "<<ATOFP[i].FullSpectrum.GetName()<<"\n"; return &ATOFP[i];
				}
				
			}
		}
		//пытаемся найти с тем же типом детектора
		for(unsigned int i=0;i<ATOFP.size();i++)
		{
			if((Type(current->DetNumber)==Type(ATOFP[i].DetNumber))&&(ATOFP[i].GeneratedAnti))
			{
				if((!Fitted))
				{
					cout<<"current: "<<current->FullSpectrum.GetName()<<" Found: "<<ATOFP[i].FullSpectrum.GetName()<<"\n"; return &ATOFP[i];
				}
				else if((ATOFP[i].FittedTOF)&&(ATOFP[i].HasManualFits()>0))
				{
					cout<<"current: "<<current->FullSpectrum.GetName()<<" Found: "<<ATOFP[i].FullSpectrum.GetName()<<"\n"; return &ATOFP[i];
				}
			}
		}
		if(!HasManualFits)
		{
			for(unsigned int i=0;i<ATOFP.size();i++)
			{
				if((current->DetNumber==ATOFP[i].DetNumber)&&(ATOFP[i].GeneratedAnti))
				{
					if(ATOFP[i].FittedTOF)
					{
						cout<<"current: "<<current->FullSpectrum.GetName()<<" Found: "<<ATOFP[i].FullSpectrum.GetName()<<"\n"; return &ATOFP[i];
					}
					
				}
			}
		}
		
	}
	return 0;
	
}

void TOFComponent::FitTOFComponent(string Type)
{
	bool UseFit=false;
	if(Type=="HPGe")
	{
		PeakPositions=TF1("PeakPositions","pol1(0)+exp([2]*log(x)+[3]*log(x)^2)",0,10000);
		PeakPositions.SetParameters(86,1e-4,1,-1e-1);
		PeakPositions.FixParameter(2,0);
		PeakPositions.FixParameter(3,0);
		PeakPositionGraph.Fit(&PeakPositions);
		PeakPositions.SetParLimits(2,0.1,1);
		PeakPositions.SetParLimits(3,-1,0);
		UseFit=true;
	}
	else if(Type=="LaBr")
	{
		PeakPositions=TF1("PeakPositions","pol2(0)",0,10000);
		UseFit=true;
	}
	if(UseFit)
	{
		Fitted=true;
		PeakPositionGraph.Fit(&PeakPositions);
	}
}

void TOFWindow::GetParametersFromTOFWindow(TOFWindow *w, bool UseTSpectrum)
{
	//FitFunction=w->FitFunction;
	LeftBorder=w->LeftBorder; RightBorder=w->RightBorder;
	//MinE=w->MinE; MaxE=w->MaxE;
	NPeaks=w->NPeaks;
	PosValues=w->PosValues;
	if(UseTSpectrum)
	{
		FindPeaksWithTSpectrum();
	}
	if(PosValues.size()>2)
	{
		double x1,x2;
		FitFunction.GetRange(x1,x2);
		PosValues[0]=x1;
		PosValues[PosValues.size()-1]=x2;
		/*for(int i=0;i<NPeaks;i++)
		{
			PosValues[i+1]=FitFunction.GetParameter(TString::Format("Pos_%d",i));
		}*/
		CreateFitFunction();
		//BuildFitFunction();
	}
	
}

void TOFWindow::GenerateNames()
{
	TOFSpectrum.SetName(fATOF->FullSpectrum.GetName()+TString::Format("_TW_%d_%d",(int)MinE,(int)MaxE));
	TOFSpectrum.SetTitle(fATOF->FullSpectrum.GetName()+TString::Format(" Time window %.1f_%.1f",MinE,MaxE));
}

int TOFWindow::FindTOFComponent(double PosValue, double &_diff_value)
{
	int Index=-1;
	double Diff=0;
	int iterator=0;
	for(unsigned int i=0;i<fATOF->TOFComponents.size();i++)
	{
		if(fATOF->TOFComponents[i].PeakPositionGraph.GetN()>0)
		{
			int NPoints=fATOF->TOFComponents[i].PeakPositionGraph.GetN();
			double x,y;
			fATOF->TOFComponents[i].PeakPositionGraph.GetPoint(NPoints-1,x,y);
			if(iterator==0)
			{
				Diff=abs(PosValue-y);
				Index=0;
			}
			if(abs(PosValue-y)<Diff)
			{
				Index=i;
				Diff=abs(PosValue-y);
			}
			iterator++;
		}
	}
	return Index;
}



void TOFWindow::SaveToRoot(TFile *f)
{
	f->WriteTObject(&TOFSpectrum);
	if(Fitted)
	{
		f->WriteTObject(&FitFunction);
	}

}

void TOFWindow::Draw()
{
	TOFSpectrum.Draw("hist");
	TLine *l=new TLine();
	if(Fitted)
	{
		for(int i=0;i<PosValues.size();i++)
		{
			l->SetLineColor(i+1);
			l->DrawLine(PosValues[i],TOFSpectrum.GetMinimum(),PosValues[i],TOFSpectrum.GetMaximum());
		}
		gStyle->SetOptFit(1111);
		//FitFunction.SetLineWidth(2);
		FitFunction.SetLineColor(2);
		FitFunction.Draw("same");
		gPad->GetCanvas()->Modified();
		gPad->GetCanvas()->Update();
	}
}

void TOFWindow::RemoveFunctionFromHist()
{
	TF1 *fun = TOFSpectrum.GetFunction(fATOF->FullSpectrum.GetName()+TString::Format("_TOFFit_%d_%d",(int)MinE,(int)MaxE));
	if(fun)
	TOFSpectrum.GetListOfFunctions()->Remove(fun);
}
void TOFWindow::BuildFitFunction()
{
	cout<<"PosValues.size():"<<PosValues.size()<<"\n";
	cout<<"Sigma!!!!!!!!!: "<<FitFunction.GetParameter(2)<<"\n";

	TOFSpectrum.Fit(&FitFunction,"RB","",PosValues[0],PosValues[PosValues.size()-1]);
	//TOFSpectrum.Fit(&FitFunction);
	//FitFunction.FixParameter(2,1000);
	//TOFSpectrum.Fit(&FitFunction);

	if((int)(fATOF->TOFComponents.size())<NPeaks)
	{
		fATOF->TOFComponents.resize(NPeaks);
		for(int i=0;i<NPeaks;i++)
		{
			if(fATOF->TOFComponents[i].PeakPositionGraph.GetN()==0)
			{
				int ParNumberPos=FitFunction.GetParNumber(TString::Format("Pos_%d",i));
				int ParNumberHeigth=FitFunction.GetParNumber(TString::Format("Heigth_%d",i));
				fATOF->TOFComponents[i].fATOF=fATOF;
				fATOF->TOFComponents[i].CompNumber=i;
				
				double E=MinE+(MaxE-MinE)/2;

				fATOF->TOFComponents[i].PeakPositionGraph.SetName(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d_Pos",i));
				fATOF->TOFComponents[i].SigmaGraph.SetName(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d_Sig",i));
				fATOF->TOFComponents[i].AmplitudeGraph.SetName(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d_Heigth",i));

				fATOF->TOFComponents[i].LeftBorderGraph.SetName(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d_LB",i));
				fATOF->TOFComponents[i].RightBorderGraph.SetName(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d_RB",i));



				fATOF->TOFComponents[i].PeakPositionGraph.SetLineColor(i+1);
				fATOF->TOFComponents[i].SigmaGraph.SetLineColor(i+1);
				fATOF->TOFComponents[i].PeakPositionGraph.SetMarkerColor(i+1);
				fATOF->TOFComponents[i].SigmaGraph.SetMarkerColor(i+1);
				fATOF->TOFComponents[i].LeftBorderGraph.SetLineColor(i+1);
				fATOF->TOFComponents[i].RightBorderGraph.SetMarkerColor(i+1);


				fATOF->TOFComponents[i].PeakPositionGraph.SetPoint(0,E,FitFunction.GetParameter(ParNumberPos));
				fATOF->TOFComponents[i].PeakPositionGraph.SetPointError(0,0,FitFunction.GetParError(ParNumberPos));
				int ParNumberSig=FitFunction.GetParNumber(TString::Format("Sigma_%d",i));
				fATOF->TOFComponents[i].SigmaGraph.SetPoint(0,E,FitFunction.GetParameter(ParNumberSig));
				fATOF->TOFComponents[i].SigmaGraph.SetPointError(0,0,FitFunction.GetParError(ParNumberSig));
				
				fATOF->TOFComponents[i].AmplitudeGraph.SetPoint(0,0,FitFunction.GetParameter(ParNumberHeigth));
				fATOF->TOFComponents[i].AmplitudeGraph.SetPointError(0,0,FitFunction.GetParError(ParNumberHeigth));
			}
		}
	}
	else
	{
		for(int i=0;i<NPeaks;i++)
		{
			double Diff=0;
			double E=MinE+(MaxE-MinE)/2;
			int ParNumberPos=FitFunction.GetParNumber(TString::Format("Pos_%d",i));
			int ParNumberSig=FitFunction.GetParNumber(TString::Format("Sigma_%d",i));
			int ParNumberHeigth=FitFunction.GetParNumber(TString::Format("Heigth_%d",i));
			int CompNumber=FindTOFComponent(FitFunction.GetParameter(ParNumberPos),Diff);
			cout<<"CompNumber:"<<CompNumber<<"\n";
			if(CompNumber>-1)
			{
				TGraphErrors *pos=&(fATOF->TOFComponents[CompNumber].PeakPositionGraph);
				TGraphErrors *sig=&(fATOF->TOFComponents[CompNumber].SigmaGraph);
				TGraphErrors *amp=&(fATOF->TOFComponents[CompNumber].AmplitudeGraph);
				
				cout<<"N="<<pos->GetN()<<"\n";
				if(ManualFit)
				{
					AddPointToTGraphErrors(pos,E,FitFunction.GetParameter(ParNumberPos),0,FitFunction.GetParError(ParNumberPos));
					AddPointToTGraphErrors(sig,E,FitFunction.GetParameter(ParNumberSig),0,FitFunction.GetParError(ParNumberSig));
					AddPointToTGraphErrors(amp,E,FitFunction.GetParameter(ParNumberHeigth),0,FitFunction.GetParError(ParNumberHeigth));
				}
			}
			
			
		}
	}
	Fitted=true;
	fATOF->FittedTOF=true;
	for(unsigned int i=0;i<fATOF->TOFWindows.size();i++)
	{
		if(!fATOF->TOFWindows[i].Fitted)
		{
			fATOF->FittedTOF=false;
		}
	}
	
}

void TOFWindow::AddPointToCompGraphAuto()
{
	bool ComFound=false;
	for(int i=0;i<NPeaks;i++)
	{
		double Diff=0;
		double E=MinE+(MaxE-MinE)/2;
		int ParNumberPos=FitFunction.GetParNumber(TString::Format("Pos_%d",i));
		int ParNumberSig=FitFunction.GetParNumber(TString::Format("Sigma_%d",i));
		int ParNumberHeigth=FitFunction.GetParNumber(TString::Format("Heigth_%d",i));
		int CompNumber=FindTOFComponent(FitFunction.GetParameter(ParNumberPos),Diff);
		cout<<"CompNumber:"<<CompNumber<<"\n";
		if(CompNumber>-1)
		{
			TGraphErrors *pos=&(fATOF->TOFComponents[CompNumber].PeakPositionGraph);
			TGraphErrors *sig=&(fATOF->TOFComponents[CompNumber].SigmaGraph);
			TGraphErrors *amp=&(fATOF->TOFComponents[CompNumber].AmplitudeGraph);
			
			cout<<"N="<<pos->GetN()<<"\n";
			AddPointToTGraphErrors(pos,E,FitFunction.GetParameter(ParNumberPos),0,FitFunction.GetParError(ParNumberPos));
			AddPointToTGraphErrors(sig,E,FitFunction.GetParameter(ParNumberSig),0,FitFunction.GetParError(ParNumberSig));
			AddPointToTGraphErrors(amp,E,FitFunction.GetParameter(ParNumberHeigth),0,FitFunction.GetParError(ParNumberHeigth));
			ComFound=true;
		}

	}
	if(!ComFound)
	{
		for(int i=0;i<NPeaks;i++)
		{
			double Diff=0;
			double E=MinE+(MaxE-MinE)/2;
			int ParNumberPos=FitFunction.GetParNumber(TString::Format("Pos_%d",i));
			int ParNumberSig=FitFunction.GetParNumber(TString::Format("Sigma_%d",i));
			int ParNumberHeigth=FitFunction.GetParNumber(TString::Format("Heigth_%d",i));
			int CompNumber=i;
			cout<<"CompNumber:"<<CompNumber<<"\n";
			if(CompNumber>-1)
			{
				TGraphErrors *pos=&(fATOF->TOFComponents[CompNumber].PeakPositionGraph);
				TGraphErrors *sig=&(fATOF->TOFComponents[CompNumber].SigmaGraph);
				TGraphErrors *amp=&(fATOF->TOFComponents[CompNumber].AmplitudeGraph);
				
				cout<<"N="<<pos->GetN()<<"\n";
				AddPointToTGraphErrors(pos,E,FitFunction.GetParameter(ParNumberPos),0,FitFunction.GetParError(ParNumberPos));
				AddPointToTGraphErrors(sig,E,FitFunction.GetParameter(ParNumberSig),0,FitFunction.GetParError(ParNumberSig));
				AddPointToTGraphErrors(amp,E,FitFunction.GetParameter(ParNumberHeigth),0,FitFunction.GetParError(ParNumberHeigth));
				ComFound=true;
			}

		}
	}
}

void TOFWindow::FindPeaksWithTSpectrum()
{
	if(PosValues.size()>3)
	{
		cout<<"FindPeaksWithTSpectrum!\n";
		int NPeaks=PosValues.size()-2;
		TSpectrum *s = new TSpectrum(PosValues.size()-2,4);
		int nfound = s->Search(&TOFSpectrum,1,"",0.01);
		Double_t *xpeaks;
		xpeaks = s->GetPositionX();
		vector<double> Found(nfound);
		for(int i=0;i<nfound;i++)
		{
			Found[i]=xpeaks[i];
		}
		sort(Found.begin(),Found.end());
		for(int i=1;i<PosValues.size()-1;i++)
		{
			if(nfound>i-1)
			{
				PosValues[i]=Found[i-1];
			}
		}
	}
	
}

void TOFWindow::CreateFitFunction()
{
	if(PosValues.size()<3)
	{
		return;
	}
	double MaxSigma=(PosValues[PosValues.size()-1]-PosValues[0])/((PosValues.size()-2)*2);
	double MinSigma=TOFSpectrum.GetBinWidth(1);
	TString FunctionTS;
	int ParIterator=0;
	NPeaks=0;
	for(unsigned int i=1;i<PosValues.size()-1;i++)
	{
		FunctionTS+=TString::Format("gaus(%d)",ParIterator);
		if(i<PosValues.size()-2)
		{
			FunctionTS+="+";
		}
		ParIterator+=3;
	}
	FunctionTS+=TString::Format("+pol0(%d)",ParIterator);
	ParIterator=0;
	cout<<"Fit:"<<FunctionTS<<"\n";
	FitFunction=TF1(fATOF->FullSpectrum.GetName()+TString::Format("_TOFFit_%d_%d",(int)MinE,(int)MaxE),FunctionTS,PosValues[0],PosValues[PosValues.size()-1]);
	if(!ManualFit)
	{
		if(PosValues.size()==3)
		{
			TOFSpectrum.GetXaxis()->SetRangeUser(PosValues[0],PosValues[PosValues.size()-1]);
			PosValues[1]=TOFSpectrum.GetXaxis()->GetBinCenter(TOFSpectrum.GetMaximumBin());
		}
	}
	for(unsigned int i=1;i<PosValues.size()-1;i++)
	{
		FitFunction.SetParName(ParIterator,TString::Format("Heigth_%d",i-1));
		FitFunction.SetParName(ParIterator+1,TString::Format("Pos_%d",i-1));
		FitFunction.SetParName(ParIterator+2,TString::Format("Sigma_%d",i-1));
		FitFunction.SetParLimits(ParIterator,0.5*TOFSpectrum.Interpolate(PosValues[i]),2*TOFSpectrum.GetMaximum());
		FitFunction.SetParameter(ParIterator,TOFSpectrum.Interpolate(PosValues[i]));
		FitFunction.SetParameter(ParIterator+1,PosValues[i]);
		//FitFunction.FixParameter(ParIterator+1,(PosValues[i]));
		double diff;
		int CompNumber=FindTOFComponent(PosValues[i],diff);
		TGraphErrors *sigGr=0;

		if(CompNumber<0)
		{
			if(fATOF->ReferenceATOF)
			{
				CompNumber=fATOF->ReferenceATOF->TOFWindows[0].FindTOFComponent(PosValues[i],diff);
				sigGr=&(fATOF->ReferenceATOF->TOFComponents[CompNumber].SigmaGraph);
			}
		}///!
		bool UseFixedLimitsForSigma=false;
		if(fATOF)
		{
			if(fATOF->UseFixedLimitsForSigma)
			{
				UseFixedLimitsForSigma=true;
				MaxSigma=fATOF->MaxTOFSigma;
				MinSigma=fATOF->MinTOFSigma;

			}
		}///!
		if((CompNumber>-1)&&(!ManualFit)&&(!UseFixedLimitsForSigma))
		{
			if(!sigGr)
			{
				sigGr=&(fATOF->TOFComponents[CompNumber].SigmaGraph);
			}
			double x,y;

			int N=sigGr->GetN();
			sigGr->GetPoint(N-1,x,y);
			double E=MinE+(MaxE-MinE)/2;
			if(x<E)
			{
				MaxSigma=2*sigGr->Eval(E);
				//MinSigma=0.5*fATOF->TOFComponents[CompNumber].Eval(E);
			}
			else
			{
				MaxSigma=2*y;
				//MinSigma=0.5*y;
			}
		}
		cout<<"MinSigma,MaxSigma: "<<ParIterator+2<<" "<<MinSigma<<" "<<MaxSigma<<"\n";
		FitFunction.SetParLimits(ParIterator+1,PosValues[i]-MaxSigma,PosValues[i]+MaxSigma);
		cout<<"MinPos,MaxPos: "<<PosValues[i]<<" "<<PosValues[i]-MaxSigma<<" "<<PosValues[i]+MaxSigma<<"\n";
		FitFunction.SetParameter(ParIterator+2,(MaxSigma-MinSigma)/2);
		FitFunction.SetParLimits(ParIterator+2,MinSigma,MaxSigma);
		NPeaks++;
		ParIterator+=3;
	}
	FitFunction.SetParameter(ParIterator,0);
	
	/*TOFSpectrum.Fit(&FitFunction,"R","",PosValues[0],PosValues[PosValues.size()-1]);
	cout<<"FitFunction:"<<FunctionTS<<"\n";
	ParIterator=0;
	double SigmaValue=TOFSpectrum.GetBinWidth(1);
	{
		for(unsigned int i=1;i<PosValues.size()-1;i++)
		{
			if(FitFunction.GetParameter(ParIterator+2)>TOFSpectrum.GetBinWidth(1))
			{
				SigmaValue=FitFunction.GetParameter(ParIterator+2);
			}
			ParIterator+=3;
		}
	}
	ParIterator=0;
	for(unsigned int i=1;i<PosValues.size()-1;i++)
	{
		if(FitFunction.GetParameter(ParIterator+2)<TOFSpectrum.GetBinWidth(1))
		{
			FitFunction.SetParLimits(ParIterator+1,(PosValues[i]-SigmaValue),(PosValues[i]+SigmaValue));
			FitFunction.SetParLimits(ParIterator+2,0.5*SigmaValue,2*SigmaValue);
		}
		else
		{
			FitFunction.SetParLimits(ParIterator+1,(PosValues[i]-FitFunction.GetParameter(ParIterator+2)),(PosValues[i]+FitFunction.GetParameter(ParIterator+2)));
		}
		ParIterator+=3;
	}*/
}

void TOFComponent::FillComponent()
{
	SpectrumHist=TH1D(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d",CompNumber),fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d; E,keV; Count",CompNumber),fATOF->FullSpectrum.GetNbinsX(),fATOF->FullSpectrum.GetXaxis()->GetXmin(),fATOF->FullSpectrum.GetXaxis()->GetXmax());
	bool Generate2d=false;
	
	if(fATOF)
	{
		if(fATOF->GenerateTH2)
		{
			Generate2d=true;
			SpectrumHist2D=TH2F(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d",CompNumber),fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d; E,keV; Count",CompNumber),fATOF->FullSpectrum.GetNbinsX(),fATOF->FullSpectrum.GetXaxis()->GetXmin(),fATOF->FullSpectrum.GetXaxis()->GetXmax(),100,-50,50);
		}
	}
	//SpectrumHist2D

	for(int i=1;i<fATOF->FullSpectrum.GetNbinsX()+1;i++)
	{
		double Pos=0;
		double E=fATOF->FullSpectrum.GetXaxis()->GetBinCenter(i);
		if(Fitted)
		{
			Pos=PeakPositions.Eval(E);
		}
		else
		{
			Pos=PeakPositionGraph.Eval(E);
		}
		double Sig=SigmaGraph.Eval(E);

		if(fATOF)
		{
			if(E<fATOF->TOFDependenceLeft)
			{
				Pos=PeakPositionGraph.Eval(fATOF->TOFDependenceLeft);
				Sig=SigmaGraph.Eval(fATOF->TOFDependenceLeft);
			}
			else if(E>fATOF->TOFDependenceRight)
			{
				Pos=PeakPositionGraph.Eval(fATOF->TOFDependenceRight);
				Sig=SigmaGraph.Eval(fATOF->TOFDependenceRight);
			}
		}


		LeftBorderGraph.SetPoint(LeftBorderGraph.GetN(),fATOF->FullSpectrum.GetXaxis()->GetBinCenter(i),Pos-NSigmaLeft*Sig);
		RightBorderGraph.SetPoint(RightBorderGraph.GetN(),fATOF->FullSpectrum.GetXaxis()->GetBinCenter(i),Pos+NSigmaRight*Sig);
		int YMin=fATOF->PureCoincedence.GetYaxis()->FindBin(Pos-NSigmaLeft*Sig);
		int YMax=fATOF->PureCoincedence.GetYaxis()->FindBin(Pos+NSigmaRight*Sig);
		double Value=0,Error=0;
		for(int j=YMin;j<=YMax;j++)
		{
			Value+=fATOF->PureCoincedence.GetBinContent(i,j);
			Error+=pow(fATOF->PureCoincedence.GetBinError(i,j),2);
			if(Generate2d)
			{
				int BinY=SpectrumHist2D.GetYaxis()->FindBin(fATOF->PureCoincedence.GetYaxis()->GetBinCenter(j));
				SpectrumHist2D.SetBinContent(i,BinY,fATOF->PureCoincedence.GetBinContent(i,j));
				SpectrumHist2D.SetBinError(i,BinY,fATOF->PureCoincedence.GetBinError(i,j));
			}
		}
		Error=sqrt(Error);
		SpectrumHist.SetBinContent(i,Value);
		SpectrumHist.SetBinError(i,Error);
	}
}

void ATOFProcess::ParseName()
{
	TString name=FullSpectrum.GetName();
	name.ReplaceAll("_"," ");
	stringstream sstr(name.Data());
	sstr>>name>>BeamNumber>>DetNumber;
}



void ATOFProcess::AutomaticTOFFit(ATOFProcess *proc)
{
	TF1::DefaultAddToGlobalList(false); 	
	TOFWindow* PrevWindow=0;
	ReferenceATOF=proc;
	bool ThereAreInitFits=false;
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		if(TOFWindows[i].ManualFit==true)
		{
			ThereAreInitFits=true;
			break;
		}
	}

	//что делаем:
	//1)есть вручную профитированные спектры-берем параметры из графика
	//2) нет-берем параметры из начального приближения для другой комбинации пучок-детектор
	if(ThereAreInitFits)
	{
		TOFWindow *refWindow=0;
		for(unsigned int i=0;i<TOFWindows.size();i++)
		{
			if(TOFWindows[i].ManualFit)
			{
				refWindow=&TOFWindows[i];
			}
			else if(refWindow)
			{
				TOFWindows[i].PosValues=refWindow->PosValues;
				double E=TOFWindows[i].MinE+(TOFWindows[i].MaxE-TOFWindows[i].MinE)/2;
				for(unsigned int j=0;j<TOFComponents.size();j++)
				{
					double x,y;
					int N=TOFComponents[j].PeakPositionGraph.GetN();
					TOFComponents[j].PeakPositionGraph.GetPoint(N-1,x,y);
					
					if(x<E)
					{
						TOFWindows[i].PosValues[j+1]=TOFComponents[j].PeakPositionGraph.Eval(E);
					}
				}
			}
		}
	}
	//else if(TOFComponents.size()>0)
	else if(proc)
	{
		//берем параметры из начального приближения для другой комбинации пучок-детектор
		for(unsigned int i=0;i<TOFWindows.size();i++)
		{
			if(proc->TOFComponents.size()>1)
			{
				if(i==0)
				{
					TOFWindows[i].GetParametersFromTOFWindow(&(proc->TOFWindows[i]),true);
					//TOFWindows[i].FindPeaksWithTSpectrum();
				}
				else
				{
					TOFWindows[i].GetParametersFromTOFWindow(&(TOFWindows[i-1]),false);//если компонент >1, то зависимость от энергии слабая
					double E=TOFWindows[i].MinE+(TOFWindows[i].MaxE-TOFWindows[i].MinE)/2;
					for(unsigned int j=0;j<TOFComponents.size();j++)
					{
						TOFWindows[i].PosValues[j+1]=TOFComponents[j].PeakPositionGraph.Eval(E);
					}
				} 
			}
			else
			{
				TOFWindows[i].GetParametersFromTOFWindow(&(proc->TOFWindows[i]),true);
			}
			//if(proc->TOFWindows.size()>i)
		}

	}

	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		if(TOFWindows[i].PosValues.size()>2)
		{
			TOFWindows[i].CreateFitFunction();
			TOFWindows[i].BuildFitFunction();
			TOFWindows[i].AddPointToCompGraphAuto();
		}
		
	}
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].FitTOFComponent(DetType);
	}
	/*for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		for(unsigned int j=0;j<TOFComponents.size();j++)
		{
			if((TOFWindows[i].PosValues.size()<j+1)&&(TOFComponents[j].Fitted==true))
			{
				double E=TOFWindows[i].MinE+(TOFWindows[i].MaxE-TOFWindows[i].MinE)/2;
				TOFWindows[i].PosValues[j+1]=TOFComponents[j].PeakPositions.Eval(E);
			}
		}
		TOFWindows[i].CreateFitFunction();
		TOFWindows[i].BuildFitFunction();
		TOFWindows[i].AddPointToCompGraphAuto();
	}
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].FitTOFComponent(DetType);
	}*/
	
}

void ATOFProcess::SaveToRoot(TFile *f)
{
	f->WriteTObject(&FullSpectrum);
	f->WriteTObject(&Coincedence);
	f->WriteTObject(&PureCoincedence);
	f->WriteTObject(&Anticoincedence);
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		TOFWindows[i].SaveToRoot(f);
	}
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].SaveToRoot(f);
	}
}

void ATOFProcess::GetParametersFromATOF(ATOFProcess *p)
{
	if(!GeneratedAnti)
	{
		GenerateAntiCoincedence(p);
	}
	for(unsigned int i=0;i<p->TOFComponents.size();i++)
	{
		TOFComponents.push_back(p->TOFComponents[i]);
		TOFComponents[i].SpectrumHist.Reset();
	}
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		if(p->TOFWindows.size()>i)
		{
			TOFWindows[i].GetParametersFromTOFWindow(&(p->TOFWindows[i]));
		}
	}
	AssignPointers(fGUI);
	CurrentWindowNumber=p->CurrentWindowNumber;
}


void ATOFProcess::GenerateTOFWindows(double WidthValue)
{
	double Min=FullSpectrum.GetXaxis()->GetXmin(), Max=FullSpectrum.GetXaxis()->GetXmax();
	double Start=Min;
	while(Start<Max)
	{
		PureCoincedence.GetXaxis()->SetRangeUser(Start,Start+WidthValue);
		TOFWindow TW;
		TW.fATOF=this;
		TW.TOFSpectrum=*(PureCoincedence.ProjectionY());
		TW.MinE=Start; TW.MaxE=Start+WidthValue;
		TW.GenerateNames();
		TW.WindowNumber=TOFWindows.size();
		TOFWindows.push_back(TW);
		Start+=WidthValue;
	}
	if(TOFWindows.size()>0)
	{
		CurrentWindow=&TOFWindows[0];
	}
	PureCoincedence.GetXaxis()->SetRangeUser(Min,Max);
}
void ATOFProcess::GenerateAntiCoincedence(double LeftBorder,double RightBorder, double CoinLeftBorder,double CoinRightBorder,bool UseLinearRegression)
{
	int YMin=FullSpectrum.GetYaxis()->FindBin(LeftBorder);
	int YMax=FullSpectrum.GetYaxis()->FindBin(RightBorder);
	
	int YMinC=FullSpectrum.GetYaxis()->FindBin(CoinLeftBorder);
	int YMaxC=FullSpectrum.GetYaxis()->FindBin(CoinRightBorder);
	
	vector<Float_t> BkgBinContents,x_values;
	
	BkgBinContents.resize(YMax-YMin);
	x_values.resize(YMax-YMin);
	
	Anticoincedence=TH2F(TString(FullSpectrum.GetName())+"_anti",TString(FullSpectrum.GetName())+"_anti",FullSpectrum.GetNbinsX(),FullSpectrum.GetXaxis()->GetXmin(),FullSpectrum.GetXaxis()->GetXmax(),YMax-YMin,LeftBorder,RightBorder);
	Coincedence=TH2F(TString(FullSpectrum.GetName())+"_coin",TString(FullSpectrum.GetName())+"_coin",FullSpectrum.GetNbinsX(),FullSpectrum.GetXaxis()->GetXmin(),FullSpectrum.GetXaxis()->GetXmax(),YMaxC-YMinC,CoinLeftBorder,CoinRightBorder);
	PureCoincedence=TH2F(TString(FullSpectrum.GetName())+"_pure",TString(FullSpectrum.GetName())+"_pure",FullSpectrum.GetNbinsX(),FullSpectrum.GetXaxis()->GetXmin(),FullSpectrum.GetXaxis()->GetXmax(),YMaxC-YMinC,CoinLeftBorder,CoinRightBorder);
	MaxPosition=Coincedence.GetYaxis()->GetBinCenter(Coincedence.ProjectionY()->GetMaximumBin());
	for(int i=1;i<=FullSpectrum.GetNbinsX();i++)
	{
		double BkgValue=0;
		int iterator=0;
		for(int j=YMin+1;j<=YMax;j++)
		{
			Anticoincedence.SetBinContent(i,j-YMin,FullSpectrum.GetBinContent(i,j));
			BkgValue+=FullSpectrum.GetBinContent(i,j);
			Anticoincedence.SetBinError(i,j-YMin,sqrt(FullSpectrum.GetBinContent(i,j)));
			BkgBinContents[j-YMin-1]=(FullSpectrum.GetBinContent(i,j));
			x_values[j-YMin-1]=(j);
			iterator++;
		}
		//cout<<"YMax-YMin="<<YMax-YMin<<" "<<iterator<<" "<<x_values.size()<<"\n";
		BkgValue=BkgValue/(YMax-YMin);
		double k=0;
		if(UseLinearRegression)
		{
			double bkg0=BkgValue;
			LinearRegression(&x_values,&BkgBinContents,0,BkgBinContents.size(),k,BkgValue);
			//cout<<"k,a="<<k<<" "<<BkgValue<<" "<<bkg0<<"\n";
		}
		//x_values.resize(0);
		//BkgBinContents.resize(0);
		
		for(int j=YMinC+1;j<=YMaxC;j++)
		{
			double Corr=BkgValue+k*j;
			Coincedence.SetBinContent(i,j-YMinC,FullSpectrum.GetBinContent(i,j));
			Coincedence.SetBinError(i,j-YMinC,sqrt(FullSpectrum.GetBinContent(i,j)));
			PureCoincedence.SetBinContent(i,j-YMinC,FullSpectrum.GetBinContent(i,j)-Corr);
			PureCoincedence.SetBinError(i,j-YMinC,sqrt(FullSpectrum.GetBinContent(i,j)+Corr));
			PeakSubstrateRatio+=BkgValue;
		}
		
	}
	GeneratedAnti=true;
	PeakSubstrateRatio=PureCoincedence.Integral()/PeakSubstrateRatio;
	Anticoincedence.Scale((CoinRightBorder-CoinLeftBorder)/(RightBorder-LeftBorder));
	GenerateTOFWindows(500);
}

void ATOFProcess::ProcessComponents()
{
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		if(FittedTOF)
		TOFComponents[i].FillComponent();
	}
}

void ATOFProcess::GenerateAntiCoincedence(ATOFProcess *proc)
{
	if(GeneratedAnti)
	{
		return;
	}
	
	
	/*double YMinAnti=proc->MaxPosition-proc->Anticoincedence.GetYaxis()->GetXmin();
	double YMaxAnti=proc->MaxPosition-proc->Anticoincedence.GetYaxis()->GetXmax();
	
	double YMin=proc->MaxPosition-proc->Coincedence.GetYaxis()->GetXmin();
	double YMax=proc->MaxPosition-proc->Coincedence.GetYaxis()->GetXmax();
	
	double MaxPosition=FullSpectrum.GetYaxis()->GetBinCenter(FullSpectrum.ProjectionY()->GetMaximumBin());*/
	
	double YMinAnti=proc->Anticoincedence.GetYaxis()->GetXmin();
	double YMaxAnti=proc->Anticoincedence.GetYaxis()->GetXmax();
	
	double YMin=proc->Coincedence.GetYaxis()->GetXmin();
	double YMax=proc->Coincedence.GetYaxis()->GetXmax();
	
	double MaxPos=proc->FullSpectrum.GetYaxis()->GetBinCenter(proc->FullSpectrum.ProjectionY()->GetMaximumBin());

	double DistMinAnti=MaxPos-YMinAnti;
	double DistMaxAnti=MaxPos-YMaxAnti;

	double DistMinCoin=MaxPos-YMin;
	double DistMaxCoin=MaxPos-YMax;

	double MaxPosition=FullSpectrum.GetYaxis()->GetBinCenter(FullSpectrum.ProjectionY()->GetMaximumBin());
	
	//GenerateAntiCoincedence(MaxPosition-YMinAnti,MaxPosition-YMaxAnti,MaxPosition-YMin,MaxPosition-YMax);
	GenerateAntiCoincedence(MaxPosition-DistMinAnti,MaxPosition-DistMaxAnti,MaxPosition-DistMinCoin,MaxPosition-DistMaxCoin);

	/*cout<<" Dist:"<<proc->FullSpectrum.GetName()<<" "<<FullSpectrum.GetName()<<" YMinAnti:"<<YMinAnti<<" YMaxAnti:"<<YMaxAnti<<" YMin:"<<YMin<<" YMax:"<<YMax<<" MaxPos:"<<MaxPos<<" DistMinAnti:"<<DistMinAnti<<" DistMaxAnti:"<<DistMaxAnti<<" "<<MaxPosition-DistMinAnti<<" "<<MaxPosition-DistMaxAnti<<" MaxPosition:"<<MaxPosition<<"\n";*/
}

void ATOFProcess::DrawInGUI()
{
	if(fGUI->DrawType==0)
	{
		TCanvas* c=fGUI->fCanvas->GetCanvas();
		c->Clear();
		c->cd();
		//c->Divide(2,1);
		//c->cd(1);
		//FullSpectrum.Draw("colz");
		//c->cd(2);
		fGUI->CurrentHist=FullSpectrum.ProjectionY();
		fGUI->CurrentHist->Draw("hist");
		Coincedence.ProjectionY()->Draw("hist same");
		PureCoincedence.ProjectionY()->Draw("hist same");
		TLine *l1=new TLine(),*l2=new TLine(),*l3=new TLine(),*l4=new TLine();
		l1->SetLineColor(1);
		l2->SetLineColor(3);
		l3->SetLineColor(4);
		l4->SetLineColor(6);
		l1->DrawLine(Anticoincedence.GetYaxis()->GetXmin(),0,Anticoincedence.GetYaxis()->GetXmin(),Anticoincedence.GetMaximum());
		l2->DrawLine(Anticoincedence.GetYaxis()->GetXmax(),0,Anticoincedence.GetYaxis()->GetXmax(),Anticoincedence.GetMaximum());
		l3->DrawLine(Coincedence.GetYaxis()->GetXmin(),0,Coincedence.GetYaxis()->GetXmin(),Coincedence.GetMaximum());
		l4->DrawLine(Coincedence.GetYaxis()->GetXmax(),0,Coincedence.GetYaxis()->GetXmax(),Coincedence.GetMaximum());
		fGUI->fCanvas->GetCanvas()->Modified();
		fGUI->fCanvas->GetCanvas()->Update();
		if((GeneratedAnti)&&(!fGUI->fDrawGraphs->GetState()))
		{
			
			fGUI->fGraphCanvas->GetCanvas()->cd();
			FullSpectrum.ProjectionX()->SetMinimum(1);
			FullSpectrum.SetLineColor(1);
			FullSpectrum.ProjectionX()->Draw("hist");
			Anticoincedence.SetLineColor(3);
			Anticoincedence.ProjectionX()->Draw("hist same");
			Coincedence.ProjectionX()->Draw("hist same");
			PureCoincedence.SetLineColor(6);
			PureCoincedence.ProjectionX()->Draw("hist same");
			gPad->BuildLegend(0.7,0.7,1,1);
			fGUI->fGraphCanvas->GetCanvas()->Modified();
			fGUI->fGraphCanvas->GetCanvas()->Update();
		}
		else if(fGUI->fDrawGraphs->GetState())
		{
			c=fGUI->fGraphCanvas->GetCanvas();
			c->Clear();
			c->cd();
			c->Divide(1,2);
			c->cd(1);
			TMultiGraph *mgrPos=new TMultiGraph();
			TMultiGraph *mgrSig=new TMultiGraph();
			for(unsigned int i=0;i<TOFComponents.size();i++)
			{
				mgrPos->Add(&(TOFComponents[i].PeakPositionGraph));
				mgrSig->Add(&(TOFComponents[i].SigmaGraph));
			}
			mgrPos->Draw("alp");
			c->cd(2);
			mgrSig->Draw("alp");
			c->cd();
			c->Modified();
			c->Update();
			fGUI->fCanvas->GetCanvas()->Modified();
			fGUI->fCanvas->GetCanvas()->Update();
			fGUI->fCanvas->GetCanvas()->cd();
		}
		c=fGUI->fCanvas->GetCanvas();
		c->cd();
	}
	else
	{
		TCanvas* c=fGUI->fCanvas->GetCanvas();
		c->cd();
		if(!CurrentWindow)
		{
			CurrentWindow=&(TOFWindows[CurrentWindowNumber]);
		}
		cout<<"WindowNumber="<<fGUI->CurrentATOF->CurrentWindowNumber<<" "<<TOFWindows.size()<<"\n";
		fGUI->CurrentHist=&(CurrentWindow->TOFSpectrum);
		CurrentWindow->Draw();
		fGUI->fCanvas->GetCanvas()->Modified();
		fGUI->fCanvas->GetCanvas()->Update();

		c=fGUI->fGraphCanvas->GetCanvas();
		c->Clear();
		c->cd();
		c->Divide(1,2);
		c->cd(1);
		TMultiGraph *mgrPos=new TMultiGraph();
		TMultiGraph *mgrSig=new TMultiGraph();
		for(unsigned int i=0;i<TOFComponents.size();i++)
		{
			mgrPos->Add(&(TOFComponents[i].PeakPositionGraph));
			mgrSig->Add(&(TOFComponents[i].SigmaGraph));
			/*if(TOFComponents[i].Fitted)
			{
				mgrPos->Add(&(TOFComponents[i].PeakPositions));
			}*/
		}
		mgrPos->Draw("alp");
		for(unsigned int i=0;i<TOFComponents.size();i++)
		{
			if(TOFComponents[i].Fitted)
			{
				TOFComponents[i].PeakPositions.Draw("same");
			}
		}
		c->cd(2);
		mgrSig->Draw("alp");
		c->cd();
		c->Modified();
		c->Update();
		fGUI->fCanvas->GetCanvas()->Modified();
		fGUI->fCanvas->GetCanvas()->Update();
		fGUI->fCanvas->GetCanvas()->cd();
	}
}

int ATOFProcess::HasManualFits()
{
	int result=0;
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		if(TOFWindows[i].ManualFit)
		{
			result++;
		}
	}
	return result;
}

void ATOFProcess::Draw()
{
	TCanvas *c=0;
	if(gPad)
	{
		c=gPad->GetCanvas();
	}
	if(!c)
	{
		c=new TCanvas();
	}
	c->Clear();
	c->Divide(2,2);
	FullSpectrum.SetLineColor(1);
	Anticoincedence.SetLineColor(3);
	PureCoincedence.SetLineColor(6);
	c->cd(1);
	FullSpectrum.ProjectionY()->Draw("hist");
	Coincedence.ProjectionY()->Draw("hist same");
	PureCoincedence.ProjectionY()->Draw("hist same");

	c->cd(2);
	FullSpectrum.ProjectionX()->Draw("hist");
	Anticoincedence.ProjectionX()->Draw("hist same");
	Coincedence.ProjectionX()->Draw("hist same");
	PureCoincedence.ProjectionX()->Draw("hist same");

	c->cd(3);
	FullSpectrum.Draw("col");
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].PeakPositionGraph.SetLineColor(i+1);
		TOFComponents[i].SigmaGraph.SetLineColor(i+1);
		TOFComponents[i].PeakPositionGraph.SetMarkerColor(i+1);
		TOFComponents[i].SigmaGraph.SetMarkerColor(i+1);
		TOFComponents[i].LeftBorderGraph.SetLineColor(i+1);
		TOFComponents[i].RightBorderGraph.SetMarkerColor(i+1);
		TOFComponents[i].LeftBorderGraph.SetLineColor(i+1);
		TOFComponents[i].RightBorderGraph.SetLineColor(i+1);
		TOFComponents[i].LeftBorderGraph.Draw("l");
		TOFComponents[i].RightBorderGraph.Draw("l");
	}

	c->cd(4);
	double max=0;
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].SpectrumHist.SetLineColor(i+1);
		if(TOFComponents[i].SpectrumHist.GetMaximum()>max)
		{
			max=TOFComponents[i].SpectrumHist.GetMaximum();
		}
	}
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].SpectrumHist.SetMaximum(max);
		TOFComponents[i].SpectrumHist.SetMinimum(0);
		if(i==0)
		{
			TOFComponents[i].SpectrumHist.Draw("hist");
		}
		else
		{
			TOFComponents[i].SpectrumHist.Draw("hist same");
		}
	}


}
void ATOFProcess::AssignPointers(GUIclass *GUI)
{
	fGUI=GUI;
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].fATOF=this;
	}
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		TOFWindows[i].fATOF=this;
	}
}

void ATOFProcess::RefitAutoFitted()
{
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].PeakPositionGraph=TGraphErrors();
		TOFComponents[i].SigmaGraph=TGraphErrors();
		TOFComponents[i].AmplitudeGraph=TGraphErrors();
		TOFComponents[i].LeftBorderGraph=TGraphErrors();
		TOFComponents[i].RightBorderGraph=TGraphErrors();

		TOFComponents[i].PeakPositionGraph.SetName(FullSpectrum.GetName()+TString::Format("_comp_%d_Pos",i));
		TOFComponents[i].SigmaGraph.SetName(FullSpectrum.GetName()+TString::Format("_comp_%d_Sig",i));
		TOFComponents[i].AmplitudeGraph.SetName(FullSpectrum.GetName()+TString::Format("_comp_%d_Heigth",i));
		TOFComponents[i].LeftBorderGraph.SetName(FullSpectrum.GetName()+TString::Format("_comp_%d_LB",i));
		TOFComponents[i].RightBorderGraph.SetName(FullSpectrum.GetName()+TString::Format("_comp_%d_RB",i));
	}

	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
	//	double E=TOFWindows[i].MinE+(TOFWindows[i].MaxE-TOFWindows[i].MinE)/2;
		if(TOFWindows[i].ManualFit)
		{
			TOFWindows[i].AddPointToCompGraphAuto();
		}
		else
		{
			TOFWindows[i].Fitted=false;
		}

	}

}

void MoveTH2F(TH2F *f1,double Mv)
{
	vector<vector<double> > Content;
	vector<vector<double> > Error;
	Content.resize(f1->GetNbinsX());
	Error.resize(f1->GetNbinsX());
	for(int i=0;i<f1->GetNbinsX();i++)
	{
		Content[i].resize(f1->GetNbinsY());
		Error[i].resize(f1->GetNbinsY());
		for(int j=0;j<f1->GetNbinsY();j++)
		{
			double BinCenterY=f1->GetYaxis()->GetBinCenter(j+1);
			double NewBinCenter=BinCenterY-Mv;
			int BinNumber=f1->GetYaxis()->FindBin(NewBinCenter);
			if(BinNumber==0)
			{
				BinNumber=1;
			}
			if(BinNumber==f1->GetNbinsY()+1)
			{
				BinNumber=f1->GetNbinsY();
			}
			Content[i][j]=f1->GetBinContent(i+1,BinNumber);
			Error[i][j]=f1->GetBinError(i+1,BinNumber);
		}
	}
	for(int i=0;i<f1->GetNbinsX();i++)
	{
		for(int j=0;j<f1->GetNbinsY();j++)
		{
			f1->SetBinContent(i,j,Content[i][j]);
			f1->SetBinError(i,j,Error[i][j]);
		}
	}
}

void AddMV(TH2F *f1, TH2F *f2, double k ,double Mv)
{
	for(int i=0;i<f1->GetNbinsX();i++)
	{
		for(int j=0;j<f1->GetNbinsY();j++)
		{
			double BinCenterY=f1->GetYaxis()->GetBinCenter(j+1);
			double NewBinCenter=BinCenterY-Mv;
			int BinNumber=f2->GetYaxis()->FindBin(NewBinCenter);
			if(BinNumber==0)
			{
				BinNumber=1;
			}
			if(BinNumber==f2->GetNbinsY()+1)
			{
				BinNumber=f2->GetNbinsY();
			}
			double BinValue=f1->GetBinContent(i+1,j+1)+k*f2->GetBinContent(i+1,BinNumber);
			double BinError=sqrt(pow(f1->GetBinError(i+1,j+1),2)+pow(k*f2->GetBinError(i+1,BinNumber),2));
			f1->SetBinContent(i+1,j+1,BinValue);
			f1->SetBinError(i+1,j+1,BinError);
		}
	}
}

void ATOFProcess::Add(ATOFProcess &p,double k,double MV)
{
	TH2F _FullSpectrum, _Anticoincedence, _Coincedence, _PureCoincedence;
	
	/*if(MV!=0)
	{
		Move(&(p.FullSpectrum),MV);
		Move(&(p.Anticoincedence),MV);
		Move(&(p.Coincedence),MV);
		Move(&(p.PureCoincedence),MV);
	}*/
	
	_FullSpectrum=p.FullSpectrum;
	_Anticoincedence=p.Anticoincedence;
	_Coincedence=p.Coincedence;
	_PureCoincedence=p.PureCoincedence;
	
	if(MV==0)
	{
		FullSpectrum.Add(&_FullSpectrum,k);
		Anticoincedence.Add(&_Anticoincedence,k);
		Coincedence.Add(&_Coincedence,k);
		PureCoincedence.Add(&_PureCoincedence,k);
	}
	else
	{
		AddMV(&FullSpectrum,&_FullSpectrum,k,MV);
		AddMV(&Coincedence,&_Coincedence,k,MV);
		AddMV(&Anticoincedence,&_Anticoincedence,k,MV);
		AddMV(&PureCoincedence,&_PureCoincedence,k,MV);
	}

	//7 vfz 2024 продолжить здесь!!!
	if(TOFComponents.size()==p.TOFComponents.size())
	{
		for(unsigned int i=0;i<p.TOFComponents.size();i++)
		{
			TOFComponents[i].SpectrumHist.Add(&(p.TOFComponents[i].SpectrumHist),k);
			TOFComponents[i].SpectrumHist2D.Add(&(p.TOFComponents[i].SpectrumHist),k);
		}
	}
	TOFComponents.resize(0);
	TOFWindows.resize(0);

}

bool ATOFSignature(string str)
{
	if(str.find("_TW_")!=string::npos)
	{
		return false;
	}
	if(str.find("_TOFFit_")!=string::npos)
	{
		return false;
	}
	
	if(str.size()>0)
	{
		int _iterator=0;
		bool digits=true;
		for(int i=str.size()-1;i>=0;i--)
		{
			if(str[i]=='_')
			{
				_iterator++;
			}
			if(_iterator<2)
			{
				if((str[i]!='_')&&((str[i]<'0')||(str[i]>'9')))
				return false;
			}
			if(_iterator==2)
			{
				return true;
			}	
		}
	}
	return false;
}

vector<string> GetVectorOfATOFNames(vector<string> &Names)//гистограмма может иметь любое имя вида %s_%d_%d
{
	vector<string> result;
	for(unsigned int i=0;i<Names.size();i++)
	{
		if(ATOFSignature(Names[i]))
		{
			result.push_back(Names[i]);
		}
	}
	return result;
}


vector<TOFWindow> ReadTOFWindows(TFile *f,string HName, vector<string> Names)
{
	vector<TOFWindow> res;
	vector<string> ObjNames;
	string Name="_TW_";
	for(unsigned int i=0;i<Names.size();i++)
	{
		int NMInd=Names[i].find(Name);
		if(NMInd==0)
		{
			ObjNames.push_back(Names[i].substr(Name.size()));
		}
	}
	for(unsigned int i=0;i<ObjNames.size();i++)
	{
		TH1D *h=(TH1D*)f->Get((HName+Name+ObjNames[i]).c_str());
		if(h)
		{
			TOFWindow result;
			h->SetDirectory(0);
			result.TOFSpectrum=*h;
			//delete h;
			TString ts(ObjNames[i].c_str());
			ts.ReplaceAll("_"," ");
			stringstream sstr(ts.Data());
			sstr>>result.MinE>>result.MaxE;
			TF1 *fft=(TF1*)f->Get((HName+"_TOFFit_"+ObjNames[i]).c_str());
			if(fft)
			{
				result.FitFunction=*fft;
				result.Fitted=true;
				fft->GetRange(result.LeftBorder,result.RightBorder);
			}
			res.push_back(result);
		}
	}
	return res;
}

ATOFProcess ReadATOF(TFile *f,string Name, vector<string> &Names)
{
	TH1::AddDirectory(false);
	ATOFProcess result;
	if(!ATOFSignature(Name))
	{
		cout<<"This is ReadATOF: name "<<Name<<" does not corresponds ATOF object!\n";
	}
	vector<string> ObjNames;
	for(unsigned int i=0;i<Names.size();i++)
	{
		int NMInd=Names[i].find(Name);
		if(NMInd==0)
		{
			ObjNames.push_back(Names[i].substr(Name.size()));
		}
	}
	TH2F *hFull=((TH2F*)f->Get(Name.c_str()));
	if(!hFull)
	{
		cout<<"Cannot read base histogram "<<Name<<"\n";
		return result;
	}

	result.FullSpectrum=*hFull;
	//delete hFull;
	for(unsigned int i=0;i<ObjNames.size();i++)
	{
		if(ObjNames[i].find("_anti")!=string::npos)
		{
			TH2F *hA=((TH2F*)f->Get((Name+ObjNames[i]).c_str()));
			cout<<"Name+ObjNames[i]:"<<Name+ObjNames[i]<<"\n";
			if(hA)
			{
				hA->SetDirectory(0);
				result.Anticoincedence=*hA;
				//delete hA;
			}
		}
		if(ObjNames[i].find("_pure")!=string::npos)
		{
			TH2F *hP=((TH2F*)f->Get((Name+ObjNames[i]).c_str()));
			if(hP)
			{
				hP->SetDirectory(0);
				result.PureCoincedence=*hP;
				//delete hP;
			}
		}
		if(ObjNames[i].find("_coin")!=string::npos)
		{
			TH2F *hP=((TH2F*)f->Get((Name+ObjNames[i]).c_str()));
			if(hP)
			{
				hP->SetDirectory(0);
				result.Coincedence=*hP;
				//delete hP;
			}
		}
		result.TOFWindows=ReadTOFWindows(f,Name,ObjNames);

	}
	return result;
}


vector<ATOFProcess> ReadATOF(TFile *f)
{
	vector<ATOFProcess> result;
	vector<string> Names;
	TList *KeyList=f->GetListOfKeys();
	TIter next(KeyList);
	while (TObject *obj = next())
	{
		Names.push_back(obj->GetName());
	}
	for(unsigned int i=0;i<Names.size();i++)
	{
		if(ATOFSignature(Names[i]))
		{
			result.push_back(ReadATOF(f,Names[i],Names));
		}
	}
	return result;
}

GUIclass::~GUIclass()
{
   // Clean up

   Cleanup();
}

//______________________________________________________________________________
void GUIclass::CloseWindow()
{
   // Called when window is closed via the window manager.

   delete this;
}

//GUIclass::GUIclass(vector<TH2F> histograms, TString _fName): TGMainFrame(gClient->GetRoot(), 100, 100)
GUIclass::GUIclass(vector<ATOFProcess> &_ATOFP, TString _fName, TString _f_time_calibr)
{
	ATOFP=_ATOFP;
	for(unsigned int i=0;i<_ATOFP.size();i++)
	{
		ATOFP[i].ParseName();
		ATOFP[i].AssignPointers(this);
	}
	fName=_fName;
	/*Spectra=_Spectra;
	for(unsigned int i=0;i<Spectra.size();i++)
	{
		Spectra[i].fGUI=this;
	}
	fName=_fName;
	//hist=_hist;

	/*LeftBorder=new TLine(); RightBorder=new TLine(); Centroid=new TLine();
	LeftBorder->SetLineColor(2); RightBorder->SetLineColor(6); Centroid->SetLineColor(3);*/

	SetCleanup(kDeepCleanup);
	fMainFrame = new TGHorizontalFrame(this, 0, 0, 0);
	AddFrame(fMainFrame, new TGLayoutHints(kLHintsRight | kLHintsExpandX |
                                   kLHintsExpandY));
	fVFrame1= new TGVerticalFrame(fMainFrame, 0, 0, 0);
	fVFrame2= new TGVerticalFrame(fMainFrame, 0, 0, 0);

	fCanvas = new TRootEmbeddedCanvas("Canvas", fVFrame1, 600, 400);
	fGraphCanvas = new TRootEmbeddedCanvas("Canvas2", fVFrame2, 600, 400);
	fLcan = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
	fMainFrame->AddFrame(fVFrame1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX | kLHintsExpandY));
	fMainFrame->AddFrame(fVFrame2, new TGLayoutHints(kLHintsRight|kLHintsExpandX | kLHintsExpandY));
	fVFrame1->AddFrame(fCanvas, fLcan);
	fVFrame2->AddFrame(fGraphCanvas, fLcan);



	fHframe0 = new TGHorizontalFrame(fVFrame1, 0, 0, 0);
	fHframe0->Resize(200, 50);
	fHframe1 = new TGHorizontalFrame(fVFrame1, 0, 0, 0);
	fHframe1->Resize(200, 25);

	fHframe2 = new TGHorizontalFrame(fVFrame1, 0, 0, 0);
	fHframe2->Resize(200, 25);
	fHframe3 = new TGHorizontalFrame(fVFrame1, 0, 0, 0);
	fHframe3->Resize(200, 25);

	fHframe4 = new TGHorizontalFrame(fVFrame1, 0, 0, 0);
	fHframe4->Resize(200, 25);

	Slider = new TGDoubleHSlider(fHframe1, 10, kDoubleScaleBoth,0);
  	Slider->SetRange(0,1);
 	Slider->SetPosition(0,1);

 	Select_CACWin=new TGTextButton(fHframe2, "Select AC/C windows");
 	PrevTimeWindow=new TGTextButton(fHframe2, "Prev Time window");
 	NextTimeWindow=new TGTextButton(fHframe2, "Next Time window");
 	AutoTOFFit=new TGTextButton(fHframe2, "Auto TOF fit");
 	FitTOF=new TGTextButton(fHframe2, "Start TOF fit");
 	RefitButton=new TGTextButton(fHframe2, "Refit");
 	EnergyForm=new TGNumberEntry(fHframe3, 0, 5);
 	LabelFitResult = new TGLabel(fHframe3, "No input.");
 	ProcessComponents=new TGTextButton(fHframe3, "Process Components");

 	PrevSpectrum=new TGTextButton(fHframe4, "<< Previous spectrum");
 	NextSpectrum=new TGTextButton(fHframe4, "Next spectrum >>");
 	Save=new TGTextButton(fHframe4, "Save");
 	fDrawGraphs=new TGCheckButton(fHframe2, "Draw graphs", 10);

 	Select_CACWin->Connect("Clicked()", "GUIclass",
                      this, "SelectCoincedenceAnticoincedence()");

 	PrevTimeWindow->Connect("Clicked()", "GUIclass",
                      this, "PrevTimeWindowFunction()");
 	NextTimeWindow->Connect("Clicked()", "GUIclass",
                      this, "NextTimeWindowFunction()");
 	PrevSpectrum->Connect("Clicked()", "GUIclass",
                      this, "PrevSpectrumFunction()");
 	NextSpectrum->Connect("Clicked()", "GUIclass",
                      this, "NextSpectrumFunction()");
    FitTOF->Connect("Clicked()", "GUIclass",
                      this, "FitTOFFunction()");
    AutoTOFFit->Connect("Clicked()", "GUIclass",
                      this, "AutoFitFunction()");
    ProcessComponents->Connect("Clicked()", "GUIclass",
                      this, "ProcessComponentsFunction()");
 	Save->Connect("Clicked()", "GUIclass",
                      this, "SaveFunction()");
 	RefitButton->Connect("Clicked()", "GUIclass",
                      this, "RefitButtonFunction()");


 	fMainLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5);
	SliderLayout = new TGLayoutHints(kLHintsLeft|kLHintsExpandX);
	
	Slider->Connect("PositionChanged()", "GUIclass",
                      this, "DoSlider()");
	fCanvas->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "GUIclass",this,"DoCanvas(Int_t,Int_t,Int_t,TObject*)");

	fHframe1->AddFrame(Slider,SliderLayout);
	fHframe2->AddFrame(Select_CACWin,SliderLayout);
	fHframe2->AddFrame(PrevTimeWindow,SliderLayout);
	fHframe2->AddFrame(NextTimeWindow,SliderLayout);
	fHframe2->AddFrame(FitTOF,SliderLayout);
	fHframe2->AddFrame(AutoTOFFit,SliderLayout);
	fHframe2->AddFrame(RefitButton,SliderLayout);
	fHframe2->AddFrame(fDrawGraphs,SliderLayout);
	//fHframe2->AddFrame(ResetCalibration,SliderLayout);
	fHframe3->AddFrame(EnergyForm,SliderLayout);
	fHframe3->AddFrame(LabelFitResult,SliderLayout);

	fHframe4->AddFrame(PrevSpectrum,SliderLayout);
	fHframe4->AddFrame(NextSpectrum,SliderLayout);
	fHframe4->AddFrame(Save,SliderLayout);

	fVFrame1->AddFrame(fHframe0, fMainLayout);
	fVFrame1->AddFrame(fHframe1, fMainLayout);
	fVFrame1->AddFrame(fHframe2, fMainLayout);
	fVFrame1->AddFrame(fHframe3, fMainLayout);
	fVFrame1->AddFrame(fHframe4, fMainLayout);
	SetWindowName("Calibration");
	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();
	if(ATOFP.size()>0)
	{
		CurrentATOF=&ATOFP[0];
		CurrentATOF->DrawInGUI();
		fCanvas->GetCanvas()->Modified();
		fCanvas->GetCanvas()->Update();
	}
}

void GUIclass::DoSlider()
{
	double Min, Max;
	Slider->GetPosition(Min,Max);
	//cout<<"Slider:"<<Min<<" "<<Max<<"\n";
	//if(hist)
	{
		/*double XMin, XMax;
		fCanvas->GetCanvas()->cd();
		XMin=hist->GetXaxis()->GetXmin();
		XMax=hist->GetXaxis()->GetXmax();
		//cout<<XMin<<" "<<XMax<<"\n";
		double Range=XMax-XMin;
		hist->GetXaxis()->SetRangeUser(XMin+Min*Range,XMin+Max*Range);
		hist->Draw();
		fCanvas->GetCanvas()->Modified();
  		fCanvas->GetCanvas()->Update();*/
	}
}
void GUIclass::DoCanvas(Int_t event, Int_t x, Int_t y, TObject *selected)
{
	TCanvas *c = (TCanvas *) gTQSender;
	double xrange = gPad->GetX2()-gPad->GetX1();
	double yrange = gPad->GetY2()-gPad->GetY1();
	
	if(ActionType==1)
	{
		TCanvas* c2=fCanvas->GetCanvas();
		c2->cd();
		double MinY=CurrentATOF->FullSpectrum.GetYaxis()->GetXmin(), MaxY=CurrentATOF->FullSpectrum.GetYaxis()->GetXmax();
		if(event==11)
		{
			cout<<"CurrentATOF->ClickIterator:"<<CurrentATOF->ClickIterator<<"\n";
			int Bin=CurrentHist->GetXaxis()->FindBin(c->AbsPixeltoX(x));
			double BinContent=CurrentHist->GetBinContent(Bin);
			cout<<"c->AbsPixeltoX(x)>MinY)&&(c->AbsPixeltoY(y)<MaxY:"<<c->AbsPixeltoX(x)<<" "<<MinY<<" "<<c->AbsPixeltoY(y)<<" "<<MaxY<<"\n";
			if((c->AbsPixeltoX(x)>MinY)&&(c->AbsPixeltoX(x)<MaxY))
			{
				if(CurrentATOF->ClickIterator==0)
				{
					TLine *l=new TLine();
					l->SetLineColor(2);
					l->DrawLine(c->AbsPixeltoX(x),CurrentHist->GetMinimum(), c->AbsPixeltoX(x),BinContent);
					CurrentATOF->AntiLeft=c->AbsPixeltoX(x);
					CurrentATOF->ClickIterator++;
				}
				else if(CurrentATOF->ClickIterator==1)
				{
					TLine *l=new TLine();
					l->SetLineColor(3);
					l->DrawLine(c->AbsPixeltoX(x),CurrentHist->GetMinimum(), c->AbsPixeltoX(x),BinContent);
					CurrentATOF->AntiRight=c->AbsPixeltoX(x);
					CurrentATOF->ClickIterator++;
				}
				else if(CurrentATOF->ClickIterator==2)
				{
					TLine *l=new TLine();
					l->SetLineColor(4);
					l->DrawLine(c->AbsPixeltoX(x),CurrentHist->GetMinimum(), c->AbsPixeltoX(x),BinContent);
					CurrentATOF->CoinLeft=c->AbsPixeltoX(x);
					CurrentATOF->ClickIterator++;
				}
				else if(CurrentATOF->ClickIterator==3)
				{
					TLine *l=new TLine();
					l->SetLineColor(6);
					l->DrawLine(c->AbsPixeltoX(x),CurrentHist->GetMinimum(), c->AbsPixeltoX(x),BinContent);
					CurrentATOF->CoinRight=c->AbsPixeltoX(x);
					ActionType=0;

					CurrentATOF->TOFWindows.resize(0);
					CurrentATOF->GenerateAntiCoincedence(CurrentATOF->AntiLeft,CurrentATOF->AntiRight,CurrentATOF->CoinLeft,CurrentATOF->CoinRight);
					CurrentATOF->ClickIterator=0;
					CurrentATOF->DrawInGUI();
				}
			}
			fCanvas->GetCanvas()->Modified();
			fCanvas->GetCanvas()->Update();
		}
		
	}
	if(ActionType==2)
	{
		TCanvas* c2=fCanvas->GetCanvas();
		c2->cd();
		double MinX=CurrentHist->GetXaxis()->GetXmin(), MaxX=CurrentHist->GetXaxis()->GetXmax();
		if(event==11)
		{
			c2->cd();
			TLine *l=new TLine();
			int Bin=CurrentHist->GetXaxis()->FindBin(c->AbsPixeltoX(x));
			double BinContent=CurrentHist->GetBinContent(Bin);
			l->SetLineColor(CurrentATOF->CurrentWindow->PosValues.size()+1);
			l->DrawLine(c->AbsPixeltoX(x),CurrentHist->GetMinimum(), c->AbsPixeltoX(x),CurrentHist->GetMaximum());
			CurrentATOF->CurrentWindow->PosValues.push_back(c->AbsPixeltoX(x));
			//CurrentATOF->DrawInGUI();
			fCanvas->GetCanvas()->Modified();
			fCanvas->GetCanvas()->Update();
		}
	}
	cout<<"w,h "<<c->AbsPixeltoX(x) <<" "<<c->AbsPixeltoY(y)<<" "<<xrange<<" "<<yrange<<"\n";
  // printf("Canvas %s: event=%d, x=%d, y=%d, selected=%s\n", c->GetName(),
    //      event, x, y, selected->IsA()->GetName());
}
/*void GUIclass::ResetCurrentFitFunction()
{

}*/

void GUIclass::SelectCoincedenceAnticoincedence()
{
	cout<<"SelectCoincedenceAnticoincedence()\n";
	if(CurrentATOF)
	{
		CurrentATOF->GeneratedAnti=false;
		DrawType=0;
		ActionType=1;
	}
	
}

void GUIclass::PrevTimeWindowFunction()
{
	DrawType=1;
	if((CurrentATOF->TOFWindows.size()>0))
	{
		if(CurrentATOF->CurrentWindowNumber>0)
		CurrentATOF->CurrentWindowNumber--;
		CurrentATOF->CurrentWindow=&(CurrentATOF->TOFWindows[CurrentATOF->CurrentWindowNumber]);
	}
	else
	{
		return;
	}
	cout<<"PrevTimeWindowFunction\n";
	CurrentATOF->DrawInGUI();
}
void GUIclass::NextTimeWindowFunction()
{
	DrawType=1;
	if(((CurrentATOF->CurrentWindowNumber+1<CurrentATOF->TOFWindows.size())&&(CurrentATOF->TOFWindows.size()>0)))
	{
		CurrentATOF->CurrentWindowNumber++;
		CurrentATOF->CurrentWindow=&(CurrentATOF->TOFWindows[CurrentATOF->CurrentWindowNumber]);
	}
	if(CurrentATOF->TOFWindows.size()==0)
	{
		return;
	}
	cout<<"NextTimeWindowFunction\n";
	CurrentATOF->DrawInGUI();
}

void GUIclass::FitTOFFunction()
{
	if(!CurrentATOF)
	{
		return;
	}
	if(!(CurrentATOF->CurrentWindow))
	{
		return;
	}
	CurrentATOF->CurrentWindow->FitNow=!CurrentATOF->CurrentWindow->FitNow;
	if(!CurrentATOF->CurrentWindow->FitNow)
	{
		FitTOF->SetText("Start TOF fit");
		CurrentATOF->CurrentWindow->ManualFit=true;
		CurrentATOF->CurrentWindow->CreateFitFunction();
		CurrentATOF->CurrentWindow->BuildFitFunction();
		CurrentATOF->DrawInGUI();
		ActionType=0;
		
	}
	else
	{
		CurrentATOF->CurrentWindow->PosValues.resize(0);
		FitTOF->SetText("Stop TOF fit");
		//CurrentATOF->CurrentWindow->BuildFitFunction();
		CurrentATOF->DrawInGUI();
		ActionType=2;
	}
}

void GUIclass::PrevSpectrumFunction()
{
	if(CurrentSpectrumNumber>0)
	{
		CurrentSpectrumNumber--;
	}
	DrawType=0;
	CurrentATOF=&ATOFP[CurrentSpectrumNumber];
	ATOFP[CurrentSpectrumNumber].DrawInGUI();
}
void GUIclass::NextSpectrumFunction()
{
	if(CurrentSpectrumNumber<(int)ATOFP.size()-1)
	{
		CurrentSpectrumNumber++;
	}
	DrawType=0;
	CurrentATOF=&ATOFP[CurrentSpectrumNumber];
	ATOFP[CurrentSpectrumNumber].DrawInGUI();
}

void SaveInThread(GUIclass *GUI)
{
	if(GUI->SavingNow)
	{
		return;
	}
	GUI->SavingNow=true;
	TFile f(GUI->fName,"Recreate");
	{
		for(unsigned int i=0;i<GUI->ATOFP.size();i++)
		{
			//ATOFP[i].SaveToRoot(&f);
			GUI->ATOFP[i].SetName(GUI->ATOFP[i].FullSpectrum.GetName());
			for(unsigned int j=0;j<GUI->ATOFP[i].TOFWindows.size();j++)
			{
				GUI->ATOFP[i].TOFWindows[j].RemoveFunctionFromHist();
			}
			f.WriteTObject(&GUI->ATOFP[i]);
		}
	}
	f.Close();
	if(GUI->f_time_calibr.Length()>0)
	{
		TFile f_time(GUI->f_time_calibr,"recreate");
		for(unsigned int i=0;i<GUI->ATOFP.size();i++)
		{
			for(unsigned int j=0;j<GUI->ATOFP[i].TOFComponents.size();j++)
			{
				GUI->ATOFP[i].TOFComponents[j].PeakPositions.SetName(TString::Format("PeakPositions_%d_%d_%d",GUI->ATOFP[i].ChannelAlpha,GUI->ATOFP[i].ChannelGamma,j));//номер х, номер гамма и номер компоненты
				GUI->ATOFP[i].TOFComponents[j].SigmaValues.SetName(TString::Format("SigmaValues_%d_%d_%d",GUI->ATOFP[i].ChannelAlpha,GUI->ATOFP[i].ChannelGamma,j));//номер х, номер гамма и номер компоненты
				f_time.WriteTObject(&(GUI->ATOFP[i].TOFComponents[j].PeakPositions));
				f_time.WriteTObject(&(GUI->ATOFP[i].TOFComponents[j].SigmaValues));
			}
		}
		f_time.Close();
	}
	GUI->SavingNow=false;
}

void GUIclass::SaveFunction()
{
	//ROOT::EnableThreadSafety();
	//async(SaveInThread,this);
	SaveInThread(this);
}
void GUIclass::ProcessComponentsFunction()
{
	for(unsigned int i=0;i<ATOFP.size();i++)
	{
		ATOFP[i].ProcessComponents();
	}
}

void GUIclass::AutoFitFunction()
{

	for(unsigned int i=0;i<ATOFP.size();i++)
	{
		CurrentATOF->AutomaticTOFFit(0);
		if(!ATOFP[i].GeneratedAnti)
		{
			ATOFProcess *p=FindBestMatch(&ATOFP[i],ATOFP);
			if(p)
			{
				ATOFP[i].GenerateAntiCoincedence(p);
				p=FindBestMatch(&ATOFP[i],ATOFP,true);
				if(p)
				ATOFP[i].AutomaticTOFFit(p);
			}
		}
		ATOFProcess *p=FindBestMatch(&ATOFP[i],ATOFP,true);
		if(p)
		{
			ATOFP[i].AutomaticTOFFit(p);
		}
	}
	
}
void GUIclass::RefitButtonFunction()
{
	if(CurrentATOF)
	{
		CurrentATOF->RefitAutoFitted();
	}
}


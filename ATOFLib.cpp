
#include "ATOFLib.hh"
#include <future>
#include <TROOT.h>   
#include <TSpectrum.h>
#pragma once


int GetIndexOfMinimum(vector<double> &v)
{
	if(v.size()==0)
	{
		return -1;
	}
	double Min=v[0];
	int MinIndex=-1;
	for(unsigned int i=0;i<v.size();i++)
	{
		if(v[i]<Min)
		{
			MinIndex=i;
			Min=v[i];
		}
	}
	return MinIndex;
}

int GetIndexOfMaximum(vector<double> &v)
{
	if(v.size()==0)
	{
		return -1;
	}
	double Max=v[0];
	int MaxIndex=-1;
	for(unsigned int i=0;i<v.size();i++)
	{
		if(v[i]>Max)
		{
			MaxIndex=i;
			Max=v[i];
		}
	}
	return MaxIndex;
}


/*std::vector<std::string> SplitString(const std::string &s, char delim, bool MergeDelimeters=true) 
{
	std::vector<std::string> result;
	std::stringstream ss (s);
	std::string item;
	while (getline (ss, item, delim)) 
	{
		if(!MergeDelimeters)
		{
			result.push_back (item);
		}
		else if (item.size()>0)
		{
			result.push_back (item);
		}
	}
	return result;
}
*/
void LinearRegression(vector<double> *x, vector<double> *y, vector<double> *x_err, vector<double> *y_err,vector<double> &result)
{
	//Митин, формула над формулой 44
	double sum1=0;  // \sum\limits_{i=0}^n x_i^2/\sigma_i^2
	double sum5=0; // \sum\limits_{i=0}^n x_i/\sigma_i^2
	double sum2=0; // \sum\limits_{i=0}^n 1/\sigma_i^2
	double sum3=0; // \sum\limits_{i=0}^n x_i*y_i/\sigma_i^2
	double sum4=0; // \sum\limits_{i=0}^n y_i/\sigma_i^2
	for(unsigned int i=0;i<x->size();i++)
	{
		double Error=0;
		if(x_err!=0)
		{
			Error+=pow(x_err->at(i),2);
		}
		if(y_err!=0)
		{
			Error+=pow(y_err->at(i),2);
		}
		if(Error==0)
		{
			Error=1;
		}
		sum2+=1/Error;
		sum1+=pow(x->at(i),2)/Error;
		sum5+=x->at(i)/Error;
		sum3+=x->at(i)*y->at(i)/Error;
		sum4+=y->at(i)/Error;
	}
	double Delta=sum1*sum2-pow(sum5,2);
	double Delta_a=sum3*sum2-sum5*sum4;
	double Delta_b=sum1*sum4-sum5*sum3;
	//cout<<"Delta="<<Delta<<"\n";
	result.resize(4);
	result[1]=Delta_a/Delta;//a
	result[0]=Delta_b/Delta;//b
	result[3]=sum2/Delta;//a
	result[2]=sum1/Delta;//b

}

TH2F TransposeTH2(TH2* h)
{
	int NBinsX=h->GetNbinsX();
	int NBinsY=h->GetNbinsY();
	
	double XMin=h->GetXaxis()->GetXmin();
	double XMax=h->GetXaxis()->GetXmax();
	
	double YMin=h->GetYaxis()->GetXmin();
	double YMax=h->GetYaxis()->GetXmax();
	
	TH2F result(TString(h->GetName())+"_T",h->GetTitle(),NBinsY,YMin,YMax,NBinsX,XMin,XMax);
	for(int i=0;i<NBinsX;i++)
	{
		for(int j=0;j<NBinsY;j++)
		{
			result.SetBinContent(j+1,i+1,h->GetBinContent(i+1,j+1));
			result.SetBinError(j+1,i+1,h->GetBinError(i+1,j+1));
		}
	}
	return result;
}

ReferenceGammaPeak::ReferenceGammaPeak(double _XMin,double _XMax,double _PeakMin,double _PeakMax,double _Energy,TH2F *_FullHist)
{
	XMin=_XMin;
	XMax=_XMax;
	PeakMin=_PeakMin;
	PeakMax=_PeakMax;
	Energy=_Energy;
	if(_FullHist)
	FullHist=*_FullHist;
}

void ReferenceGammaPeak::GenerateSubstrateHistogram()
{
	if(!fProcess)
	{
		cout<<"This is ReferenceGammaPeak::GenerateSubstrateHistogram(): cannot generate histograms because of invalid pointer to ATOFProcess* fProcess. Returned!\n";
		return;
	}
	TH2F *h_full=&(fProcess->PureCoincedence);
	//cout<<"1:"<<h_full->GetYaxis()->GetBinCenter(1)<<"\n";
	
	int x1Bin=h_full->GetXaxis()->FindBin(XMin);
	int x2Bin=h_full->GetXaxis()->FindBin(XMax);
	
	int x1BinPeak=h_full->GetXaxis()->FindBin(PeakMin);
	int x2BinPeak=h_full->GetXaxis()->FindBin(PeakMax);
	
	//int NBinsX=h_full->GetNbinsX();
	int NBinsY=h_full->GetNbinsY();
	double BinWidthY=h_full->GetYaxis()->GetBinWidth(1);
	double BinWidthX=h_full->GetXaxis()->GetBinWidth(1);
	
	//cout<<"2:"<<h_full->GetYaxis()->GetBinCenter(1)<<"\n";
	
	FullY=*(h_full->ProjectionY());  SubstrateY=*(h_full->ProjectionY());  PeakY=*(h_full->ProjectionY());
	TString OneDimName(h_full->GetName());
	
	
	FullX=TH1D(OneDimName+"_full_px",OneDimName+"_full_px",x2Bin-x1Bin+1,h_full->GetXaxis()->GetBinCenter(x1Bin)-0.5*BinWidthX,h_full->GetXaxis()->GetBinCenter(x2Bin)+0.5*BinWidthX);
	
	//cout<<"3:"<<h_full->GetYaxis()->GetBinCenter(1)<<"\n";
	SubstrateX=FullX; PeakX=FullX;
	
	OneDimName.ReplaceAll("_pure","_peak_%d");
	OneDimName=TString::Format(OneDimName,(int)(Energy*10));
	
	FullX.SetName(OneDimName+"_full_px");
	FullY.SetName(OneDimName+"_full_py");
	SubstrateX.SetName(OneDimName+"_substrate_px");
	SubstrateY.SetName(OneDimName+"_substrate_py");
	PeakX.SetName(OneDimName+"_peak_px");
	PeakY.SetName(OneDimName+"_peak_py");
	
	PeakX.Reset();
	PeakY.Reset();
	
	SubstrateX.Reset();
	SubstrateY.Reset();
	//cout<<"4:"<<h_full->GetYaxis()->GetBinCenter(1)<<"\n";
	
	if(Use2Dhist)
	{
		//cout<<"5:"<<h_full->GetYaxis()->GetBinCenter(1)<<"\n";
		//cout<<"FullHist: "<<XMin<<" "<<XMax<<" "<<h_full->GetYaxis()->GetBinCenter(1)-BinWidthY<<" "<<h_full->GetYaxis()->GetBinCenter(NBinsY+1)+BinWidthY<<"\n";
		FullHist=CutTH2(h_full,XMin,XMax,h_full->GetYaxis()->GetBinCenter(1)-BinWidthY,h_full->GetYaxis()->GetBinCenter(NBinsY+1)+BinWidthY);
		TString HistName(h_full->GetName());
		HistName.ReplaceAll("_pure","_peak_%d");
		HistName=TString::Format(HistName,(int)(Energy*10));
		FullHist.SetName(HistName);
		
		SubstrateHist=FullHist;
		SubstrateHist.SetName(HistName+"_substrate");
		PeakHist=FullHist;
		PeakHist.SetName(HistName+"_peak");
	}
	for(int i=1;i<NBinsY+1;i++)
	{
		vector<double> x_val,y_val,y_err;
		double BinCenterY=h_full->GetYaxis()->GetBinCenter(i);
		for(int j=x1Bin;j<x2Bin;j++)
		{
			if(j>=x1BinPeak && j<=x2BinPeak)
			{
				continue;
			}
			double BinCenterX=h_full->GetXaxis()->GetBinCenter(j);
			double Y=0,Y_err=0;
			if(Averaging>0)
			{
				int left=i-Averaging/2;
				int right=i+Averaging/2;
				int Count=0;
				if(left<1)
				{
					left=1;
				}
				if(right>NBinsY)
				{
					right=NBinsY;
				}
				for(int p=left;p<=right;p++)
				{
					Y+=h_full->GetBinContent(j,p);
					Y_err+=pow(h_full->GetBinError(j,p),2);
					Count++;
				}
				Y=Y/Count;
				Y_err=sqrt(Y_err)/sqrt(Count);
			}
			else
			{
				Y=h_full->GetBinContent(j,i);
				Y_err=h_full->GetBinError(j,i);
			}
			
			x_val.push_back(h_full->GetXaxis()->GetBinCenter(j));
			y_val.push_back(Y);
			y_err.push_back(Y_err);
			if(Use2Dhist)
			{
				PeakHist.SetBinContent(j-x1Bin+1,i,0);
				PeakHist.SetBinError(j-x1Bin+1,i,0);
				SubstrateHist.SetBinContent(j-x1Bin+1,i,h_full->GetBinContent(j,i));
				SubstrateHist.SetBinError(j-x1Bin+1,i,h_full->GetBinError(j,i));
			}
			int XBin=FullX.GetXaxis()->FindBin(BinCenterX);
			int YBin=FullY.GetXaxis()->FindBin(BinCenterY);
			FullX.SetBinContent(XBin,FullX.GetBinContent(XBin)+h_full->GetBinContent(j,i));
			SubstrateX.SetBinContent(XBin,SubstrateX.GetBinContent(XBin)+h_full->GetBinContent(j,i));
			FullY.SetBinContent(YBin,FullY.GetBinContent(YBin)+h_full->GetBinContent(j,i));
			SubstrateY.SetBinContent(YBin,SubstrateY.GetBinContent(YBin)+h_full->GetBinContent(j,i));
			
			FullX.SetBinError(XBin,FullX.GetBinError(XBin)+pow(h_full->GetBinError(j,i),2));
			SubstrateX.SetBinError(XBin,SubstrateX.GetBinError(XBin)+pow(h_full->GetBinError(j,i),2));
			FullY.SetBinError(YBin,FullY.GetBinError(YBin)+pow(h_full->GetBinError(j,i),2));
			SubstrateY.SetBinError(YBin,SubstrateY.GetBinError(YBin)+pow(h_full->GetBinError(j,i),2));
			
			//PeakX.SetBinContent(j-x1Bin+1,h_full->GetBinContent(j,i)+PeakX.GetBinContent(j-x1Bin+1));
			//PeakY
		}
		//cout<<"x1BinPeak: "<<x1BinPeak<<" x2BinPeak: "<<x2BinPeak<<"\n";
		vector<double> result;
		LinearRegression(&x_val,&y_val,0,&y_err,result);
		for(int j=x1BinPeak;j<=x2BinPeak;j++)
		{
			double x=h_full->GetXaxis()->GetBinCenter(j);
			double y_substrate=result[0]+h_full->GetXaxis()->GetBinCenter(j)*result[1];
			double y=h_full->GetBinContent(j,i)-y_substrate;
			
			
			int XBin=FullX.GetXaxis()->FindBin(x);
			int YBin=FullY.GetXaxis()->FindBin(BinCenterY);
			
			double Error=h_full->GetBinError(j,i);
			if(Use2Dhist)
			{
				
				PeakHist.SetBinContent(XBin,i,y);
				PeakHist.SetBinError(XBin,i,h_full->GetBinError(j,i));
				SubstrateHist.SetBinContent(XBin,i,y_substrate);
				SubstrateHist.SetBinError(XBin,i,Error);
			}
			
			FullX.SetBinContent(XBin,FullX.GetBinContent(XBin)+h_full->GetBinContent(j,i));
			SubstrateX.SetBinContent(XBin,SubstrateX.GetBinContent(XBin)+y_substrate);
			FullY.SetBinContent(YBin,FullY.GetBinContent(YBin)+h_full->GetBinContent(j,i));
			SubstrateY.SetBinContent(YBin,SubstrateY.GetBinContent(YBin)+y_substrate);
			PeakX.SetBinContent(XBin,PeakX.GetBinContent(XBin)+y);
			PeakY.SetBinContent(YBin,PeakY.GetBinContent(YBin)+y);
			
			
			FullX.SetBinError(XBin,FullX.GetBinError(XBin)+pow(h_full->GetBinError(j,i),2));
			SubstrateX.SetBinError(XBin,SubstrateX.GetBinError(XBin)+pow(Error,2));
			FullY.SetBinError(YBin,FullY.GetBinError(YBin)+pow(Error,2));
			SubstrateY.SetBinError(YBin,SubstrateY.GetBinError(YBin)+pow(Error,2));
			PeakX.SetBinError(XBin,PeakX.GetBinError(XBin)+pow(Error,2));
			PeakY.SetBinError(YBin,PeakY.GetBinError(YBin)+pow(Error,2));
		}
	}
	for(int i=1;i<=FullX.GetNbinsX();i++)
	{
		FullX.SetBinError(i,sqrt(FullX.GetBinError(i)));
		SubstrateX.SetBinError(i,sqrt(SubstrateX.GetBinError(i)));
		PeakX.SetBinError(i,sqrt(PeakX.GetBinError(i)));
	}
	for(int i=1;i<=FullY.GetNbinsX();i++)
	{
		FullY.SetBinError(i,sqrt(FullY.GetBinError(i)));
		SubstrateY.SetBinError(i,sqrt(SubstrateY.GetBinError(i)));
		PeakY.SetBinError(i,sqrt(PeakY.GetBinError(i)));
	}
}

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
	if(indexMax>=int(sData->size())-1)
	{
		indexMax=sData->size()-1;
	}
	int n=0;
	//cout<<"indexMin indexMax:"<<indexMin<<" "<<indexMax<<"\n";
	if((indexMin<int(sData->size()))&&(indexMax<int(sData->size()))&&(indexMin<indexMax))
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
			for(int i=0;i<gr->GetN()-1;i++)
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
				for(int i=0;i<gr->GetN()-1;i++)
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
void TOFComponent::FitReferencePeaks(TFitFunction *fit)
{
	
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
void TOFWindow::GenerateFunctionComponents()
{
	vector<int> Colors={1,2,3,4,6,7,8,9,25,28,30,32,41,52};
	string FuncStr(FitFunction.GetTitle());
	vector<string> Elements=SplitString(FuncStr,'+');
	int ParIter=0;
	int NComponents=0;
	for(unsigned int i=0;i<Elements.size();i++)
	{
		if(Elements[i].find("gaus(")!=string::npos)
		{
			NComponents++;
		}
	}
	if(Components.size()<NComponents)
	{
		Components.resize(NComponents);
		for(int i=0;i<NComponents;i++)
		{
			Components[i]=TF1(TString::Format("Peak_%d",i),"gaus(0)",FitFunction.GetXmin(),FitFunction.GetXmax());
		}
	}
	for(unsigned int i=0;i<Elements.size();i++)
	{
		if(Elements[i].find("gaus(")!=string::npos)
		{
			Components[i].SetParameters(FitFunction.GetParameter(ParIter),FitFunction.GetParameter(ParIter+1),FitFunction.GetParameter(ParIter+2));
			for(int pp=0;pp<3;pp++)
			{
				Components[i].SetParError(pp,FitFunction.GetParError(ParIter+pp));
			}
			Components[i].SetLineColor(Colors[i]);
			ParIter+=3;
		}
	}
}

/*class PeakParameters
{
	public:
	double H,H_err,Pos,Pos_err,Sig,Sig_err;
	
};*/

void TOFWindow::AttachFitFunction(TF1 *PrevFit,double PosDelta,double kSigma_min,double kSigma_max)
{
	string FuncStr(PrevFit->GetTitle());
	vector<string> Elements=SplitString(FuncStr,'+');
	vector<double> Parameters;
	vector<double> PeakHeigths;
	vector<double> PeakPositions;
	vector<double> PeakSigmas;
	
	vector<double> Coefficients;//показывает долю в конкретном пике от максимума
	TF1 gaus("gaus","gaus",-100,100);
	TF1 Substrate;
	for(int i=0;i<PrevFit->GetNpar();i++)
	{
		Parameters.push_back(PrevFit->GetParameter(i));
	}
	int ParIter=0;
	Components.resize(0);
	for(unsigned int i=0;i<Elements.size();i++)
	{
		if(Elements[i].find("gaus(")!=string::npos)
		{
			//cout<<"Elements[i]:"<<Elements[i]<<"\n";
			PeakHeigths.push_back(Parameters[ParIter]);
			PeakPositions.push_back(Parameters[ParIter+1]);
			PeakSigmas.push_back(Parameters[ParIter+2]);
			ParIter+=3;
			Components.push_back(TF1(TString::Format("Peak_%d",i),"gaus(0)",PrevFit->GetXmin(),PrevFit->GetXmax()));
		}
		else if(Elements[i].find("pol")!=string::npos)
		{
			vector<string> tmp=SplitString(Elements[i],'(');
			Substrate=TF1("Substrate",tmp[0].c_str(),PrevFit->GetXmin(),PrevFit->GetXmax());
			//cout<<"Elements[i]:"<<Elements[i]<<"\n";
			for(int p=0;p<Substrate.GetNpar();p++)
			{
				Substrate.SetParameter(p,Parameters[ParIter]);
				//cout<<"p:"<<p<<" "<<Parameters[ParIter]<<"\n";
				ParIter++;
			}
		}
	}
	
	
	for(unsigned int i=0;i<PeakHeigths.size();i++)
	{
		double Sum=0;
		for(unsigned int j=0;j<PeakHeigths.size();j++)
		{
			gaus.SetParameters(PeakHeigths[j],PeakPositions[j],PeakSigmas[j]);
			Sum+=gaus.Eval(PeakPositions[i]);
		}
		Coefficients.push_back(PeakHeigths[i]/Sum);
	}
	
	double MinSubstrate=PrevFit->GetXmax(),MaxSubstrate=PrevFit->GetXmin();
	for(unsigned int i=0;i<PeakPositions.size();i++)
	{
		if(MinSubstrate>PeakPositions[i]-4*PeakSigmas[i])
		{
			MinSubstrate=PeakPositions[i]-4*PeakSigmas[i];
		}
		if(MaxSubstrate<PeakPositions[i]+4*PeakSigmas[i])
		{
			MaxSubstrate=PeakPositions[i]+4*PeakSigmas[i];
		}
	}
	
	//int MaxIndex=GetIndexOfMaximum(PeakHeigths);
	
	TGraphErrors gr;
	
	int MinBin=TOFSpectrum.FindBin(PrevFit->GetXmin());
	int MaxBin=TOFSpectrum.FindBin(PrevFit->GetXmax());
	TH1D Peaks("Peaks","Peaks",MaxBin-MinBin,PrevFit->GetXmin(),PrevFit->GetXmax());
	
	for(int i=MinBin;i<=MaxBin;i++)
	{
		int N=gr.GetN();
		double Center=TOFSpectrum.GetBinCenter(i);
		if(Center<MinSubstrate||Center>MaxSubstrate)
		{
			gr.SetPoint(N,TOFSpectrum.GetBinCenter(i),TOFSpectrum.GetBinContent(i));
			gr.SetPointError(N,0,TOFSpectrum.GetBinError(i));
		}
	}
	gr.Fit(&Substrate,"R","",PrevFit->GetXmin(),PrevFit->GetXmax());

	
	for(int i=MinBin;i<=MaxBin;i++)
	{
		Peaks.SetBinContent(i+1-MinBin,TOFSpectrum.GetBinContent(i)-Substrate.Eval(TOFSpectrum.GetBinCenter(i)));
		Peaks.SetBinError(i+1-MinBin,TOFSpectrum.GetBinError(i));
	}

	TF1 NewFit(TString(PrevFit->GetName())+"_1",FuncStr.c_str(),PrevFit->GetXmin(),PrevFit->GetXmax());
	ParIter=0;
	for(unsigned int i=0;i<PeakHeigths.size();i++)
	{
		NewFit.SetParameter(ParIter,Coefficients[i]*PeakHeigths[i]);
		NewFit.SetParameter(ParIter+1,PeakPositions[i]);
		NewFit.SetParameter(ParIter+2,PeakSigmas[i]);
		NewFit.SetParLimits(ParIter,0,5*Coefficients[i]*PeakHeigths[i]);
		NewFit.SetParLimits(ParIter+1,PeakPositions[i]-PosDelta,PeakPositions[i]+PosDelta);
		NewFit.SetParLimits(ParIter+2,kSigma_min*PeakSigmas[i],kSigma_max*PeakSigmas[i]);
		
		ParIter+=3;
	}
	for(int i=ParIter;i<NewFit.GetNpar();i++)
	{
		NewFit.SetParameter(i,Substrate.GetParameter(i-ParIter));
	}
	FitFunction=NewFit;
}

void TOFWindow::FitWindow()
{
	TOFSpectrum.Fit(&(FitFunction),"R","",FitFunction.GetXmin(),FitFunction.GetXmax());
	Fitted=true;
	int ParIndex=0;
	for(unsigned int i=0;i<Components.size();i++)
	{
		Components[i].SetParameter(0,FitFunction.GetParameter(ParIndex));
		Components[i].SetParameter(1,FitFunction.GetParameter(ParIndex+1));
		Components[i].SetParameter(2,FitFunction.GetParameter(ParIndex+2));
		
		Components[i].SetParError(0,FitFunction.GetParError(ParIndex));
		Components[i].SetParError(1,FitFunction.GetParError(ParIndex+1));
		Components[i].SetParError(2,FitFunction.GetParError(ParIndex+2));
		ParIndex+=3;
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

void TOFWindow::Draw(Option_t * 	option)
{
	GenerateFunctionComponents();
	TOFSpectrum.Draw(option);
	
	TLine *l=new TLine();
	if(Fitted)
	{
		for(unsigned int i=0;i<PosValues.size();i++)
		{
			l->SetLineColor(i+1);
			l->DrawLine(PosValues[i],TOFSpectrum.GetMinimum(),PosValues[i],TOFSpectrum.GetMaximum());
		}
		gStyle->SetOptFit(1111);
		//FitFunction.SetLineWidth(2);
		FitFunction.SetLineColor(2);
		FitFunction.Draw("same"+TString(option));
		
		for(unsigned int i=0;i<Components.size();i++)
		{
			Components[i].Draw("same"+TString(option));
		}
		
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
			//double Diff=0;
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
		//int NPeaks=PosValues.size()-2;
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
		for(unsigned int i=1;i<PosValues.size()-1;i++)
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

void TOFComponent::FillComponent(string AnalysisType)
{
	if(AnalysisType=="")
	{
		AnalysisType=PredefAnalysisType;
	}
	SpectrumHist=TH1D(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d",CompNumber),fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d; E,keV; Count",CompNumber),fATOF->FullSpectrum.GetNbinsX(),fATOF->FullSpectrum.GetXaxis()->GetXmin(),fATOF->FullSpectrum.GetXaxis()->GetXmax());
	bool Generate2d=false;
	//cout<<"AnalysisType="<<AnalysisType<<"\n";
	if(fATOF)
	{
		if(fATOF->GenerateTH2)
		{
			Generate2d=true;
			SpectrumHist2D=TH2F(fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d2D",CompNumber),fATOF->FullSpectrum.GetName()+TString::Format("_comp_%d; E,keV; Count",CompNumber),fATOF->FullSpectrum.GetNbinsX(),fATOF->FullSpectrum.GetXaxis()->GetXmin(),fATOF->FullSpectrum.GetXaxis()->GetXmax(),fATOF->PureCoincedence.GetNbinsY(),fATOF->PureCoincedence.GetYaxis()->GetXmin(),fATOF->PureCoincedence.GetYaxis()->GetXmax());
		}
	}
	//SpectrumHist2D
	int AType=0;//0- BordersFit, 1-PeakSigmaFit,2-BordersGraph, 3-PeakSigmaGraph
	if(AnalysisType == "BordersFit")
	{
		AType=0;
	}
	else if(AnalysisType == "PeakSigmaFit")
	{
		AType=1;
	}
	else if(AnalysisType == "BordersGraph")
	{
		AType=2;
	}
	else if(AnalysisType == "PeakSigmaGraph")
	{
		AType=3;
	}
	//cout<<"AType="<<AType<<"\n";
	
	for(int i=1;i<fATOF->FullSpectrum.GetNbinsX()+1;i++)
	{
		double E=fATOF->FullSpectrum.GetXaxis()->GetBinCenter(i);
		double LeftBorderValue=0,RightBorderValue=0;
		
		if(fATOF)
		{
			if(E<fATOF->TOFDependenceLeft)
			{
				E=fATOF->TOFDependenceLeft;
			}
			else if(E>fATOF->TOFDependenceRight)
			{
				E=fATOF->TOFDependenceRight;
			}
		}
		
		
		if(AType==0)
		{
			LeftBorderValue=LeftBordersFit.Eval(E);
			RightBorderValue=LeftBordersFit.Eval(E);
		}
		else if(AType==1)
		{
			double Pos=0;
			double Sig=0;
			
			Pos=PeakPositions.Eval(E);
			Sig=SigmaValues.Eval(E);
			
			LeftBorderValue=Pos-Sig*NSigmaLeft;
			RightBorderValue=Pos+Sig*NSigmaRight;
		}
		else if(AType==2)
		{
			LeftBorderValue=LeftBorderGraph.Eval(E);
			RightBorderValue=RightBorderGraph.Eval(E);
		}
		else if(AType==3)
		{
			double Pos=0;
			double Sig=0;
			
			Pos=PeakPositionGraph.Eval(E);
			Sig=SigmaGraph.Eval(E);
			
			LeftBorderValue=Pos-Sig*NSigmaLeft;
			RightBorderValue=Pos+Sig*NSigmaRight;
		}
		/*double Pos=0;
		double Sig=0;
		
		if(Fitted)
		{
			Pos=PeakPositions.Eval(E);
			Sig=SigmaValues.Eval(E);
		}
		else
		{
			Pos=PeakPositionGraph.Eval(E);
			SigmaGraph.Eval(E);
		}
		

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
		*/
		
		int YMin=fATOF->PureCoincedence.GetYaxis()->FindBin(LeftBorderValue);
		int YMax=fATOF->PureCoincedence.GetYaxis()->FindBin(RightBorderValue);
		
		//cout<<"YMin="<<LeftBorderValue<<" YMax="<<RightBorderValue<<"\n";
		
		double Value=0,Error=0;
		for(int j=YMin;j<=YMax;j++)
		{
			//cout<<"Value="<<Value<<" Error="<<Error<<"\n";
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

void ATOFProcess::AttachFitFunction(TF1 *PrevFit,double PosDelta,double kSigma_min,double kSigma_max)
{
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		TOFWindows[i].AttachFitFunction(PrevFit,PosDelta,kSigma_min,kSigma_max);
	}
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

bool CheckEscapeInterval(double E,vector<vector<double> > *Escape=0)
{
	if(!Escape)
	{
		return false;
	}
	for(unsigned int i=0;i<Escape->size();i++)
	{
		if(Escape->at(i)[0]>E && Escape->at(i)[1]<E)
		{
			return true;
		}
		if(Escape->at(i)[1]>E)
		{
			return false;
		}
	}
	return false;
}

void ATOFProcess::GenerateTOFWindows(vector<vector<double> > *Windows,vector<vector<double> > *Escape)
{
	TOFWindows.resize(0);
	TH1D *pX=PureCoincedence.ProjectionX();
	TH1D pY=*(PureCoincedence.ProjectionY());
	
	for(unsigned int i=0;i<Windows->size();i++)
	{
		
		pY.Reset();
		int MinX=pX->GetXaxis()->FindBin(Windows->at(i)[0]);
		int MaxX=pX->GetXaxis()->FindBin(Windows->at(i)[1]);
		
		//cout<<"borders: "<<MinX<<" "<<MaxX<<"\n";
		double ECentroid=0,ValuesSum=0;
		for(int p=MinX;p<=MaxX;p++)
		{
			if(!CheckEscapeInterval(pX->GetBinCenter(p),Escape))
			{
				ECentroid+=pX->GetBinContent(p)*pX->GetBinCenter(p);
				ValuesSum+=pX->GetBinContent(p);
				
				for(int j=1;j<=pY.GetNbinsX();j++)
				{
					double Val=pY.GetBinContent(j), Err=pY.GetBinError(j);
					pY.SetBinContent(j,Val+PureCoincedence.GetBinContent(p,j));
					pY.SetBinError(j,Err+pow(PureCoincedence.GetBinError(p,j),2));
					
				}
			}
			
		}
		for(int j=1;j<=pY.GetNbinsX();j++)
		{
			pY.SetBinError(j,sqrt(pY.GetBinError(j)));
		}
		ECentroid=ECentroid/ValuesSum;
		TOFWindow TW;
		TW.fATOF=this;
		TW.TOFSpectrum=pY;
		TW.TOFSpectrum.SetName(TString::Format("%s_TW_%d",FullSpectrum.GetName(),int(TOFWindows.size())));
		TW.MinE=Windows->at(i)[0]; TW.MaxE=Windows->at(i)[1]; TW.Centroid=ECentroid;
		TOFWindows.push_back(TW);
	}
}

void ATOFProcess::GenerateTOFWindows(double WidthValue)
{
	TOFWindows.resize(0);
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
			//double bkg0=BkgValue;
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
	//GenerateTOFWindows(500);
}


void ATOFProcess::GenerateAntiCoincedence_GaussianSubstrate(double LeftBorder,double RightBorder, double CoinLeftBorder,double CoinRightBorder,double Thr_min,double Thr_max)
{
	double Xmin=FullSpectrum.GetXaxis()->GetXmin();
	double Xmax=FullSpectrum.GetXaxis()->GetXmax();
	
	TF1 Gauss_substrate("gs","gaus(0)+pol1(3)",-1000,1000);
	FullSpectrum.GetXaxis()->SetRangeUser(Thr_min,Thr_max);
	TH1D hp=*(FullSpectrum.ProjectionY());
	hp.GetXaxis()->SetRangeUser(LeftBorder,RightBorder);
	double Integral_S=hp.Integral();//среднее значение в бине
	double c0=hp.Interpolate(LeftBorder+1);
	double h=hp.GetMaximum()-c0;
	if(h<0)
	{
		h=0;
	}
	Gauss_substrate.SetParameters(h,hp.GetBinCenter(hp.GetMaximumBin()),(RightBorder-LeftBorder)/6,c0,0);
	hp.Fit(&Gauss_substrate,"QR","",LeftBorder,RightBorder);
	FullSpectrum.GetXaxis()->SetRangeUser(Xmin,Xmax);
	/*hp.Draw("e hist");
	hp.GetXaxis()->SetRangeUser(LeftBorder,CoinRightBorder);
	Gauss_substrate.Draw("same");
	gPad->GetCanvas()->Print("test.pdf","pdf");*/
	
	
	int YMin=FullSpectrum.GetYaxis()->FindBin(LeftBorder);
	int YMax=FullSpectrum.GetYaxis()->FindBin(RightBorder);
	
	int YMinC=FullSpectrum.GetYaxis()->FindBin(CoinLeftBorder);
	int YMaxC=FullSpectrum.GetYaxis()->FindBin(CoinRightBorder);
	
	//нужен коэффициент, позволяющий считать значение в бине по известной функции, описывающей подложку	
	
	
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
		//BkgValue=BkgValue/(YMax-YMin);
		//cout<<"YMax-YMin="<<YMax-YMin<<" "<<iterator<<" "<<x_values.size()<<"\n";
		
		for(int j=YMinC+1;j<=YMaxC;j++)
		{
			double BinCenter=FullSpectrum.GetYaxis()->GetBinCenter(j);
			double Corr=Gauss_substrate.Eval(BinCenter)*BkgValue/Integral_S;
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
	//GenerateTOFWindows(500);
}

void ATOFProcess::ProcessComponents()
{
	/*for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		TOFWindows[i].FitWindow();
	}*/
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		//if(FittedTOF)
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

void ATOFProcess::GenerateComponents(double ReferencePosition)
{
	for(unsigned int i=0;i<TOFWindows.size();i++)
	{
		SelectPeaksFCN(&(TOFWindows[i]),ReferencePosition);
	}
	for(unsigned int i=0;i<TOFComponents.size();i++)
	{
		TOFComponents[i].FillComponent("");
	}
}

void ATOFProcess::ReadFromTFile_Full(TFile *f,TString hName)
{
	FullSpectrum=*((TH2F*)f->Get(hName));
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

void ATOFProcess::AddReferencePeak(double XMin, double XMax,double PeakMin,double PeakMax,double Energy,int CompNumber)//метод, добавляющий опорный гамма-пик во временную компоненту
{
	//cout<<"TOFComponents.size(): "<<TOFComponents.size()<<"\n";
	
	ReferenceGammaPeak RGP;
	RGP.XMin=XMin;
	RGP.XMax=XMax;
	RGP.PeakMin=PeakMin;
	RGP.PeakMax=PeakMax;
	RGP.Energy=Energy;
	RGP.fProcess=this;
	RGP.Use2Dhist=Use2DhistForReferencePeaks;
	RGP.GenerateSubstrateHistogram();
	if(CompNumber>-1)
	{
		if(int(TOFComponents.size())<=CompNumber)
		{
			TOFComponents.resize(CompNumber+1);
		}
		TOFComponents[CompNumber].ReferencePeaks.push_back(RGP);
	}
}

void MoveTH2F(TH2 *f1,double Mv)
{
	vector<vector<double> > Content;
	vector<vector<double> > Error;
	Content.resize(f1->GetNbinsX());
	Error.resize(f1->GetNbinsX());
	double BinWidthY=f1->GetYaxis()->GetBinWidth(1);
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
				continue;
				//BinNumber=1;
			}
			if(BinNumber==f1->GetNbinsY()+1)
			{
				//BinNumber=f1->GetNbinsY();
				continue;
			}
			
			Content[i][j]=f1->GetBinContent(i+1,BinNumber);
			Error[i][j]=f1->GetBinError(i+1,BinNumber);
		}
	}
	f1->Reset();
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

TH2F CutTH2(TH2F *h1,double x1,double x2,double y1,double y2)
{
	int x1Bin=h1->GetXaxis()->FindBin(x1);
	int x2Bin=h1->GetXaxis()->FindBin(x2);
	int y1Bin=h1->GetYaxis()->FindBin(y1);
	int y2Bin=h1->GetYaxis()->FindBin(y2);
	
	if(x1Bin==0)
	{
		x1Bin=1;
	}
	if(x2Bin==h1->GetNbinsX())
	{
		x2Bin=h1->GetNbinsX()-1;
	}
	if(y1Bin==0)
	{
		y1Bin=1;
	}
	if(y2Bin==h1->GetNbinsY())
	{
		y2Bin=h1->GetNbinsY()-1;
	}
	double BinWidthX=h1->GetXaxis()->GetBinWidth(1);
	double BinWidthY=h1->GetYaxis()->GetBinWidth(1);
	
	x1=h1->GetXaxis()->GetBinCenter(x1Bin)-BinWidthX*0.5;
	x2=h1->GetXaxis()->GetBinCenter(x2Bin)+BinWidthX*0.5;
	
	y1=h1->GetYaxis()->GetBinCenter(y1Bin)-BinWidthY*0.5;
	y2=h1->GetYaxis()->GetBinCenter(y2Bin)+BinWidthY*0.5;
	
	int NX=x2Bin-x1Bin+1;
	int NY=y2Bin-y1Bin+1;
	
	TH2F result(TString(h1->GetName())+"_cut",TString(h1->GetTitle())+"_cut",NX,x1,x2,NY,y1,y2);
	int x_iter=1;
	for(int i=x1Bin;i<=x2Bin;i++)
	{
		int y_iter=1;
		for(int j=y1Bin;j<=y2Bin;j++)
		{
			result.SetBinContent(x_iter,y_iter,h1->GetBinContent(i,j));
			result.SetBinError(x_iter,y_iter,h1->GetBinError(i,j));
			y_iter++;
		}
		x_iter++;
	}
	return result;
}

TOFComponent* ATOFProcess::GetOrCreateTOFComponent(int CompNumber)
{
	if(TOFComponents.size()<=CompNumber)
	{
		TOFComponents.resize(CompNumber+1);
		TOFComponents[CompNumber].CompNumber=CompNumber;
		
	}
	TOFComponents[CompNumber].fATOF=this;
	return &TOFComponents[CompNumber];
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

	LeftBorder=new TLine(); RightBorder=new TLine(); Centroid=new TLine();
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
		//double MinX=CurrentHist->GetXaxis()->GetXmin();
		double MaxX=CurrentHist->GetXaxis()->GetXmax();
		if(event==11)
		{
			c2->cd();
			TLine *l=new TLine();
			int Bin=CurrentHist->GetXaxis()->FindBin(c->AbsPixeltoX(x));
			//double BinContent=CurrentHist->GetBinContent(Bin);
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


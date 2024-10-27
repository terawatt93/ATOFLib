#include "ATOFLib.hh"
class TH2Test
{
	public:
	TH2F h_full, h_substrate,h_peak;
	void GenerateSubstrateHistogram(double x1,double x2,int Averaging=6);//x1-левая граница пика, x2-правая
};

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


void TH2Test::GenerateSubstrateHistogram(double x1,double x2,int Averaging)
{
	h_substrate=h_full;
	h_substrate.SetName("substrate");
	h_peak=h_full;
	h_peak.SetName("peak");
	int x1Bin=h_substrate.GetXaxis()->FindBin(x1);
	int x2Bin=h_substrate.GetXaxis()->FindBin(x2);
	
	int NBinsX=h_substrate.GetNbinsX();
	int NBinsY=h_substrate.GetNbinsY();
	for(int i=1;i<NBinsY+1;i++)
	{
		vector<double> x_val,y_val,y_err;
		for(int j=1;j<NBinsX+1;j++)
		{
			if(j>=x1Bin && j<=x2Bin)
			{
				continue;
			}
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
					Y+=h_full.GetBinContent(j,p);
					Y_err+=pow(h_full.GetBinError(j,p),2);
					Count++;
				}
				Y=Y/Count;
				Y_err=sqrt(Y_err)/sqrt(Count);
			}
			else
			{
				Y=h_full.GetBinContent(j,i);
				Y_err=h_full.GetBinError(j,i);
			}
			
			x_val.push_back(h_full.GetXaxis()->GetBinCenter(j));
			y_val.push_back(Y);
			y_err.push_back(Y_err);
			
			h_peak.SetBinContent(j,i,0);
			h_peak.SetBinError(j,i,0);
		}
		vector<double> result;
		LinearRegression(&x_val,&y_val,0,&y_err,result);
		for(int j=x1Bin;j<=x2Bin;j++)
		{
			double x=h_full.GetXaxis()->GetBinCenter(j);
			h_substrate.SetBinContent(j,i,result[0]+h_full.GetXaxis()->GetBinCenter(j)*result[1]);
			//double Error=sqrt(pow(result[2],2)+pow(result[3]*x,2));
			double Error=h_full.GetBinError(j,i);
			h_substrate.SetBinError(j,i,Error);
			h_peak.SetBinContent(j,i,h_full.GetBinContent(j,i)-(result[0]+h_full.GetXaxis()->GetBinCenter(j)*result[1]));
		}
	}
}

void EventsInPeaks()
{
	TFile f("~/Server/ZFSRAID/Ing27-HPGe/root/pure/test159_HPGe_2_Labr_4_Ing27_84prc_Generator_100hz_Ba_center_Pb_1mm_SiO2_2cm.root");
	
	TH2F *h_sp=(TH2F*)f.Get("h_4_2_pure");
	TH2Test hh;
	hh.h_full=CutTH2(h_sp,1600,1900,-20,50);
	//hh.h_full.RebinY(20);
	hh.GenerateSubstrateHistogram(1720,1850);
	TCanvas c;
	hh.h_substrate.SetLineColor(2);
	hh.h_full.ProjectionX()->Draw("e hist");
	hh.h_substrate.ProjectionX()->Draw("e hist same");
	c.Print("test.pdf","pdf");
	TFile ff("test.root","recreate");
	ff.WriteTObject(&hh.h_full);
	ff.WriteTObject(&hh.h_substrate);
	ff.WriteTObject(&hh.h_peak);
	ff.Close();
}


#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TGTextEntry.h"
#include "TGDoubleSlider.h"
#include <TGLabel.h>
#include <TGraphErrors.h>
#include <TGNumberEntry.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <sstream>
//#include <TGCheckButton.h>
#pragma once
#define VERSION 8
using namespace std;

class GUIclass;

class ATOFProcess;

#include "TalysLib.hh"


class TOFComponent:public TNamed
{
	public:
	int CompNumber=0;
	TH1D SpectrumHist;
	TH2F SpectrumHist2D;
	TGraphErrors PeakPositionGraph, SigmaGraph, AmplitudeGraph;
	TF1 PeakPositions, SigmaValues;
	TGraphErrors LeftBorderGraph, RightBorderGraph;
	double NSigmaLeft=2.5,NSigmaRight=2.5;//потом добавить поле для редактирования
	ATOFProcess *fATOF=0;//!
	bool Fitted=false;
	void FillComponent();
	void SaveToRoot(TFile *f)
	{
		f->WriteTObject(&SpectrumHist);
		f->WriteTObject(&PeakPositionGraph);
		f->WriteTObject(&SigmaGraph);
		f->WriteTObject(&AmplitudeGraph);
	}
	void FitTOFComponent(string Type="HPGe");
	ClassDef(TOFComponent,VERSION);

};


class TOFWindow:public TNamed
{
	public:
	TH1D TOFSpectrum;
	TF1 FitFunction;//!
	double LeftBorder=0, RightBorder=0;
	bool Fitted=false;
	bool ManualFit=false;
	double MinE=0, MaxE=0;
	bool FitNow=false;
	vector<double> PosValues;//PosValues[0]-левая граница, PosValues[PosValues.size()-1], остальные-положения пиков
	ATOFProcess *fATOF=0;//!
	int NPeaks=0;
	int WindowNumber=0;
	void FindPeaksWithTSpectrum();
	bool UseTSpectrum=false;
	void Draw();
	void BuildFitFunction();
	void CreateFitFunction();
	void GenerateNames();
	void SaveToRoot(TFile *f);
	void RemoveFunctionFromHist();
	int FindTOFComponent(double PosValue, double &_diff_value);
	void GetParametersFromTOFWindow(TOFWindow *w, bool UseTSpectrum=false);
	void AddPointToCompGraphAuto();
	ClassDef(TOFWindow,VERSION);
};

class ATOFProcess:public TNamed
{
	public:
	TH2F FullSpectrum, Anticoincedence, Coincedence, PureCoincedence;
	int ChannelAlpha=0, ChannelGamma=0;

	bool UseFixedLimitsForSigma=false;//использовать фиксированные пределы для сигма на временных пиках
	double MinTOFSigma=0,MaxTOFSigma=0;//пределы по сигма для временных спектров ()
	bool GenerateTH2=false;
	TH1D NormAnticoincedence;
	double MaxPosition;
	string DetType;
	vector<TOFComponent> TOFComponents;
	vector<TOFWindow> TOFWindows;
	TOFWindow *CurrentWindow=0;
	bool GeneratedAnti=false;
	bool FittedTOF=false;
	void GenerateTOFWindows(double WidthValue);
	double AntiLeft,AntiRight,CoinLeft,CoinRight;
	void GenerateAntiCoincedence(double LeftBorder,double RightBorder, double CoinLeftBorder,double CoinRightBorder,bool UseLinearRegression=false);
	void DrawInGUI();
	void Draw();
	int CurrentWindowNumber=0;
	int ClickIterator=0;
	bool ProcIterator=0;//0-ничего не делать, 1-выбор подложки
	void SaveToRoot(TFile *f);
	int NFitted=0;
	void AutomaticTOFFit(ATOFProcess *proc);
	GUIclass *fGUI;//!
	int BeamNumber, DetNumber;
	void ParseName();
	void GenerateAntiCoincedence(ATOFProcess *proc);
	void ProcessComponents();
	void GetParametersFromATOF(ATOFProcess *p);
	void AssignPointers(GUIclass *GUI=0);
	int HasManualFits();
	void RefitAutoFitted();
	void Add(ATOFProcess &p,double k=-1);
	double PeakSubstrateRatio=0;
	double TOFDependenceLeft=0,TOFDependenceRight=10000;//границы для временных окон по энергии: если энергия лежит вне этих окон, то ширина и положение фиксируются по ближайшей границе 
	ATOFProcess *ReferenceATOF=0;//!
	ClassDef(ATOFProcess,VERSION);
};
vector<ATOFProcess> ReadATOFFromRootFile(TFile *f);

class GUIclass:public TGMainFrame
{
	public:
	TRootEmbeddedCanvas *fCanvas, *fGraphCanvas;
	TGLayoutHints       *fLcan;
	TF1                 *fFitFcn;
	TGHorizontalFrame   *fHframe0, *fHframe1, *fHframe2, *fHframe3, *fHframe4, *fMainFrame;
	TGVerticalFrame *fVFrame1,*fVFrame2;
	TGLayoutHints       *SliderLayout, *fMainLayout;
	TGDoubleHSlider *Slider;
	TGTextButton *Select_CACWin,*PrevTimeWindow,*NextTimeWindow, *FitTOF, *exit, *NextSpectrum, *PrevSpectrum,*ResetCalibration, *Save, *AutoTOFFit, *ProcessComponents,*RefitButton;
	TGNumberEntry       *EnergyForm;
	TGLabel             *LabelFitResult;
	TGCheckButton       *fDrawGraphs;
	TH1 *CurrentHist=0;
	vector<ATOFProcess> ATOFP;
	ATOFProcess *CurrentATOF=0;
	int CurrentSpectrumNumber=0;
	TString fName, f_time_calibr;
	TLine *LeftBorder=0, *RightBorder=0, *Centroid=0;
	int DrawType=0;//что рисовать: 0-гистограмму и проекцию, 1-временные окна
	int ActionType=0;//0-ничего не делать, 1-поиск совпадений/антисовпадений, 2-фит временных спектров
	//GUIclass(vector<TH2F> histograms, TString _fName);
	GUIclass(vector<ATOFProcess> &_ATOFP, TString _fName, TString _f_time_calibr="");
	virtual ~GUIclass();
	void CloseWindow();
	void DoText(const char *text);
	void DoSlider();
	bool FitNow=false;
	void DoCanvas(Int_t event, Int_t x, Int_t y, TObject *selected);
	void SelectCoincedenceAnticoincedence();
	void PrevTimeWindowFunction();
	void NextTimeWindowFunction();
	void ResetCalibrationFunction();
	void PrevSpectrumFunction();
	void NextSpectrumFunction();
	void FitTOFFunction();
	void CleanupLines();
	void SaveFunction();
	void AutoFitFunction();
	void RefitButtonFunction();
	void ProcessComponentsFunction();
	bool SavingNow=false;
	ClassDef(GUIclass, 0)
};

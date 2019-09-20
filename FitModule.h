#ifndef __WIDTHFITTER__
#define __WIDTHFITTER__

#include <TNamed.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TPolyMarker.h>
#include <TList.h>

#include <vector>
using namespace std;

const int N_BINS = 400;
const int WIDTH_MIN = 0;
const int WIDTH_MAX = 2500;
const int CHARGE_MIN = 0;
const int CHARGE_MAX = 350;
const int N_ITERATIONS_BG = N_BINS/10;

class WidthFitter : public TNamed {

	int 						id;					
	int 						day;			
	int 						hour;	
	int 						hType;	
	const char 					*namePostfix;		
	int 						counts;

	TH1 						*initH;				
	TPolyMarker					*polyMarker;
	TH1 						*bgH;				
	TH1 						*woBgH;				
	vector<float> 				*rawPeakPositions;	
	vector<float> 				*peakPositions;	
	vector<TF1*> 				*fitFuncs;			
	vector<float> 				*hWidthes;
	TH1 						*stub;				
	TF1 						*totFit;			
	TH1							*approxH;

	int 						firstAppr;		//index of first 'physical' peak
	TF1 						*z1_z2_peaks[2];
	float 						relativePairArea;   //!	
	TH1 						*tempH;				//!

	TH1 						*freqHist;
	// TSpectrum 					*sp;				
	long 						spP; 			//pointer on spectrum saved in int variable				
public:
	static const int N_PEAKS_1;	
	static const int N_PEAKS_2;	
	static const float SP_RES_1;
	static const float SP_RES_2;
	static const float SP_SIGMA_1;	
	static const float SP_SIGMA_2;	
	static const float SP_THRESH_1;	
	static const float SP_THRESH_2;	
	static const int N_ITERATIONS_FLATTER_1;
	static const int N_ITERATIONS_FLATTER_2;
	
	static const float HW_COEF;	
	static const float FREQ_PERCENTAGE;	
	void setFrequencyHist(TH1 *h) {
		freqHist = h;
	}
	bool isCentralCell() {
		// if (id==52 || id==53 || id==54 || id==55 || id==64 || id==67 || id==76 || id==79 || id==88 || id==89 || id==90 || id==91) return true;
		if (counts>0.9*freqHist->GetMaximum()) return true; 
		return false;
	}
	WidthFitter() {
		id=day=hour=0;
		namePostfix=0;
		counts = 0;
		hType = 0;

		polyMarker=0;
		initH=bgH=woBgH=stub=approxH=0;
		rawPeakPositions=peakPositions=0;
		fitFuncs=0;
		hWidthes=0;
		totFit=0;

		z1_z2_peaks[0]=0;
		z1_z2_peaks[1]=0;
	
		freqHist=0;
		spP=0;
	}
	// type==0 - width from be dst
	// type==1 - charge from be dst
	// type==2 - width from co dst
	WidthFitter(int cell, int day, int hour, TSpectrum *sp, int hType=0) {
		this->id=cell;
		this->day=day;
		this->hour=hour;
		this->spP=(long)sp;
		this->hType = hType;

		counts=0;
		const char *temp;
		switch (hType) {
			case 0: { 
				temp = "";
				break; 
			}
			case 1: { 
				temp = "Charge"; 	
				break; 
			}
			case 2: { 
				temp = "Cosmic";	
				break; 
			}
		}
		namePostfix = Form("cell%s%d_day%d_hour%d",temp,id,day,hour);
		this->SetName(namePostfix);

		initH=new TH1F("mock1","",10,0,10);
		bgH=new TH1F("mock2","",10,0,10);
		woBgH=new TH1F("mock3","",10,0,10);
		stub=new TH1F("mock4","",10,0,10);
		approxH=new TH1F("mock5","",10,0,10);
	
		rawPeakPositions=new vector<float>;
		peakPositions=new vector<float>;
		fitFuncs = new vector<TF1*>;
		hWidthes = new vector<float>;
	
		totFit = new TF1("mock6","gaus");
		initH = new TH1F(namePostfix, "", N_BINS, (hType==1)?CHARGE_MIN:WIDTH_MIN, (hType==1)?CHARGE_MAX:WIDTH_MAX);
		printf("initH initialized with name: %s\n",initH->GetName());
		polyMarker = new TPolyMarker;

		z1_z2_peaks[0]=new TF1("mock7","gaus");
		z1_z2_peaks[1]=new TF1("mock8","gaus");
	
		freqHist = new TH1F("mock9","",10,0,10);
	}
	int getType() {
		return hType;
	}
	void initSp(TSpectrum *p) {
		spP = (long)p;
	}
	void fill(float w) {
		initH->Fill(w);
	}
	void setCounts(int c) {
		counts = c;
	}
	void flatInit() {
		printf("\t-------> FLAT INIT <---------\n");
		
		TH1* old = initH;
		initH=((TSpectrum*)spP)->Background(initH,(hType==2)?N_ITERATIONS_FLATTER_2:N_ITERATIONS_FLATTER_1);
		initH->SetName(old->GetName());
		delete old;
	}
	void subtractBg() {
		printf("\t-------> subtractBg <---------\n");
		bgH = ((TSpectrum*)spP)->Background(initH,N_ITERATIONS_BG);
		woBgH=(TH1*)initH->Clone(Form("woBg_%s",namePostfix));
		woBgH->Add(bgH,-1);
	} 
	void doSpectr() {
		printf("\t-------> doSpectr <---------\n");
		int n = ((TSpectrum*)spP)->Search(woBgH, (hType==2)?SP_SIGMA_2:SP_SIGMA_1, "", (hType==2)?SP_THRESH_2:SP_THRESH_1);
		float *pos = ((TSpectrum*)spP)->GetPositionX();	
		int ind[N_PEAKS_1];
		TMath::Sort(n,pos,ind,0);
		rawPeakPositions->clear();
		for (int i=0; i<n; i++) rawPeakPositions->push_back(pos[ind[i]]);
		polyMarker=(TPolyMarker*)woBgH->GetListOfFunctions()->At(0);
		
		printf("Cell %d rawPeakPositions: ", id);
		for (float &rawPos : *rawPeakPositions) printf("%7.2f",rawPos);
		printf("\n");
	}
	void prepareFitFuncs() {
		printf("\t-------> prepareFitFuncs <---------\n");
		
		int n = rawPeakPositions->size();
		fitFuncs->clear();
		hWidthes->clear();
		for (int p=0; p<n; p++) {
			fitFuncs->push_back(new TF1(Form("peak%d_%s",p+1,namePostfix), "gaus",WIDTH_MIN,WIDTH_MAX));
			float ledge = (p==0) ? WIDTH_MIN : (rawPeakPositions->at(p-1) + rawPeakPositions->at(p))/2;
			float hw = rawPeakPositions->at(p)-ledge;
			float redge = (p==n-1) ? WIDTH_MAX : (rawPeakPositions->at(p) + rawPeakPositions->at(p+1))/2;
			hw = (redge-rawPeakPositions->at(p)<hw) ? redge-rawPeakPositions->at(p) : hw;
			fitFuncs->at(p)->SetParLimits(1, rawPeakPositions->at(p)-hw*HW_COEF, rawPeakPositions->at(p)+hw*HW_COEF);
			fitFuncs->at(p)->SetParLimits(2, hw*0.1, hw*HW_COEF);
			printf("ledge: %.2f \t redge: %.2f \t hw: %.2f\n",ledge, redge, hw);
			printf("setted mean limits: [%.2f, %.2f]\n",rawPeakPositions->at(p)-hw*HW_COEF, rawPeakPositions->at(p)+hw*HW_COEF);
			// fitFuncs->at(p)->SetParLimits(2, 0.01*hw, hw*HW_COEF);
			// fitFuncs->at(p)->SetParLimits(2, 1, 100);
			hWidthes->push_back(hw);
		}
	}	

	TF1* nextFitStep(int p) {
		printf("\t-------> nextFitStep <---------\n");
		
		TF1 *clone = (TF1*)fitFuncs->at(p)->Clone();
		printf("Cell %d Fitting in range [%.2f,%.2f]\n",id,rawPeakPositions->at(p)-hWidthes->at(p)*HW_COEF, rawPeakPositions->at(p)+hWidthes->at(p)*HW_COEF);
		double meanMin, meanMax;
		fitFuncs->at(p)->GetParLimits(1,meanMin,meanMax);
		printf("            meanRange:   [%.2f,%.2f]\n",meanMin, meanMax);
		printf("entries in stub: %f\n",stub->GetEntries());
		printf("name of fitFuncs[%d]: %s\n",p, fitFuncs->at(p)->GetName());
		stub->Fit(fitFuncs->at(p),"qB","",rawPeakPositions->at(p)-hWidthes->at(p)*HW_COEF, rawPeakPositions->at(p)+hWidthes->at(p)*HW_COEF);
		printf("height of fitFuncs[%d]: %.2f\n",p , fitFuncs->at(p)->GetParameter(0));
		printf("mean of stub: %.2f\n", stub->GetMean());
		peakPositions->push_back(fitFuncs->at(p)->GetParameter(1));
		stub->Add(fitFuncs->at(p),-1);
		return clone;
	}
	TH1 *getStub() {
		return stub;
	}

	void doFit(int steps=0) {

		peakPositions->clear();
		int n = rawPeakPositions->size();
		printf("-----doFit(): n = %d\n",n);
		stub = (TH1*)initH->Clone(Form("stub_%s",namePostfix));
		int limit=(steps==0)?n:steps;
		for (int p=0; p<limit; p++) {
			nextFitStep(p);
		}
		printf("NUMBER OF PEAKS: %d\n",peakPositions->size());
		if (!peakPositions->size()) return;
		
		if (hType==2) {
			firstAppr=1;
			return;
		}
		int shift=0;
		shift += (peakPositions->at(0)<400) ? 1 : 0;
		if (peakPositions->size()-shift == 0) {
			firstAppr=1;
			return;
		}
   		if (isCentralCell()) {
   			firstAppr=shift;
   		} else {
   			firstAppr=shift+1;
   		}
   		/*
		int maxPpos = 0;
		float maxPval = 0;
		for (int p = shift; p<peakPositions->size(); p++) if (fitFuncs->at(p)->GetMaximum()>maxPval) {maxPval=fitFuncs->at(p)->GetMaximum(); maxPpos=p;}
		shift += (shift == maxPpos) ? 0 : 1;
		firstAppr = shift;
		*/
	}
	void constructTotFit() {
		printf("\t-------> constructTotFit <---------\n");
		
		string formula;
		int n = rawPeakPositions->size();
		double pars[N_PEAKS_1*3];
		for (int p=0; p<n; p++) {
			formula+= (p==0) ? Form("gaus(%d)",p*3) : Form("+gaus(%d)",p*3);
			fitFuncs->at(p)->GetParameters(&pars[p*3]);
			float mean = fitFuncs->at(p)->GetParameter(1);
			fitFuncs->at(p)->SetRange(mean-hWidthes->at(p), mean+hWidthes->at(p));
			fitFuncs->at(p)->SetLineColor(kPink);
			fitFuncs->at(p)->SetLineWidth(1.5);
		}
		totFit = new TF1(Form("totFit_%s",namePostfix), formula.data(), WIDTH_MIN, WIDTH_MAX);
		totFit->SetParameters(pars);
		printf("totFit->Integral(0,2500): %.2f\n",totFit->Integral(0,2500));
		approxH = (TH1*)bgH->Clone(Form("approx_%s",namePostfix));
		approxH->Add(totFit);
	}

	int countAppropPeaks() {
		// printf("\n peakPositions->size()=%d\n",peakPositions->size());
		return peakPositions->size() - firstAppr;
	}

	float approxDiff() {
		TH1 *tempH = (TH1*)woBgH->Clone("temp");
		tempH->Add(totFit,-1);
		for (int bin=1; bin<=tempH->GetXaxis()->GetNbins(); bin++) {
			float val = tempH->GetBinContent(bin);
			if (val<0.0) tempH->SetBinContent(bin, -val);
		}
		float res = tempH->Integral();
		delete tempH;
		return res;
	}

	// search in all fitted peaks laying to the right of 'magic' peak for the pair integral of which is max
	float findBetterCombination() {
		if ( id<208 && countAppropPeaks()<2 || countAppropPeaks()<1 ) return 0;

		float area1=0;
		float area2=0;
		float totalArea = woBgH->Integral()*(WIDTH_MAX-WIDTH_MIN)/N_BINS;
		float maxPairArea=0;
		if (firstAppr+1==peakPositions->size()) {
			z1_z2_peaks[0]= fitFuncs->at(firstAppr);
			float mean = fitFuncs->at(firstAppr)->GetParameter(1);
			maxPairArea = fitFuncs->at(firstAppr)->Integral(mean-hWidthes->at(firstAppr), mean+hWidthes->at(firstAppr));
		} else {
			for (int p1=firstAppr; p1<peakPositions->size(); p1++) {
				float mean = fitFuncs->at(p1)->GetParameter(1);
				area1 = fitFuncs->at(p1)->Integral(mean-hWidthes->at(p1), mean+hWidthes->at(p1));
				for (int p2=p1+1; p2<peakPositions->size(); p2++) {
					mean = fitFuncs->at(p2)->GetParameter(1);
					area2 = fitFuncs->at(p2)->Integral(mean-hWidthes->at(p2), mean+hWidthes->at(p2));
					if (area1+area2>maxPairArea) {
						z1_z2_peaks[0]=fitFuncs->at(p1); 
						z1_z2_peaks[1]=fitFuncs->at(p2); 
						maxPairArea=area1+area2;
					}
				}
			}
		}
		// printf("Cell %d: maxPairArea: %.2f total: %.2f\n",id,maxPairArea, totalArea);
		relativePairArea = maxPairArea/totalArea;
		return relativePairArea;
	}

	TH1 *getInit() {
		return initH;
	}
	int getCounts() {
		return counts;
	}
	TH1 *getApprox() {
		return approxH;
	}
	TH1 *getBg() {
		return bgH;
	}
	TH1 *getWoBg() {
		return woBgH;
	}
	TPolyMarker *getPM() {
		return polyMarker;
	}
	TF1 **getZ1Z2peaks() {
		return z1_z2_peaks;
	}
	TF1 *getFitFunc(int p) {
		return fitFuncs->at(p);
	}
	vector<float> *getPeakPositions() {
		return peakPositions;
	}

	bool checkNullPointer() {
		if (!namePostfix) return false;
		if (!initH) return false;
		if (!bgH) return false;
		if (!woBgH) return false;
		if (!rawPeakPositions) return false;
		if (!fitFuncs) return false;
		if (!stub) return false;
		if (!totFit) return false;
		if (!approxH) return false;
		return true;
	}

	ClassDef(WidthFitter,1);
};


#endif 

#include <TFile.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TH2F.h>

#include "hloop.h"
#include "hcategorymanager.h"
#include "hwallraw.h"
#include "hwallhit.h"

#include "FitModule.h"
const int N_PADS = 9;
const int GREEN = TColor::GetColor("#c4ffd4");
const int RED = TColor::GetColor("#ffc4ed");

class Calibration {
	static const int N_CELLS = 304;
	const char *DEFAULT_FILE_NAME = "calibration.root";

	TFile 			*ioF;
	WidthFitter  	*cells[N_CELLS];
	WidthFitter  	*cells_charge[N_CELLS];
	WidthFitter  	*cells_cosmic[N_CELLS];
	float  			adc_slopes[N_CELLS];
	float  			adc_offsets[N_CELLS];

	HCategory 		*rawCat;
	HCategory 		*hitCat;
	int 			evts;
	TH1				*relPairAreaH;
	TH1 			*freqHist;

	static TSpectrum sp1;
	static TSpectrum sp2;
	TCanvas			*canva;
public:
	Calibration(const char *fname) {
		ioF=new TFile(fname,"update");
		for (int id=0; id<N_CELLS; id++) {
			cells[id]=(WidthFitter*)ioF->FindObjectAny(Form("cell%d_day63_hour0",id));
			cells_charge[id]=(WidthFitter*)ioF->FindObjectAny(Form("cellCharge%d_day63_hour0",id));
			cells_cosmic[id]=(WidthFitter*)ioF->FindObjectAny(Form("cellCosmic%d_day63_hour0",id));
		}
		relPairAreaH = new TH1F("relPairArea","",100, 0, 1);
		freqHist = new TH1F("frequency", "", N_CELLS, 0,N_CELLS);
		for (int id=0; id<N_CELLS; id++) {
			cells[id]		->initSp(&sp1);
			cells_charge[id]->initSp(&sp1);
			cells_cosmic[id]->initSp(&sp2);
		}
		canva = new TCanvas("c1","",1200,800);
	}
	Calibration() {
		ioF=new TFile(DEFAULT_FILE_NAME,"recreate");
		for (int id=0; id<N_CELLS; id++) {
			cells[id]=new WidthFitter(id,63,0,&sp1);
			cells_charge[id]=new WidthFitter(id,63,0,&sp1,1);
			cells_cosmic[id]=new WidthFitter(id,63,0,&sp2,2);
		}
		freqHist = new TH1F("frequency", "", N_CELLS, 0,N_CELLS);
		canva = new TCanvas("c1","",1200,800);
	}	
	~Calibration() {
		ioF->Close();
		if (gLoop) {
			gLoop->Delete();
			gLoop=NULL;
		}
	}
	void initLoop(const char *fname) {
		new HLoop(1);
		gLoop->addMultFiles(fname);
		gLoop->setInput("-*,+HWallRaw,+HWallHit");
		evts = gLoop->getChain()->GetEntries();
		rawCat=gLoop->getCategory("HWallRaw");
		hitCat=gLoop->getCategory("HWallHit");
	}
	void loopBe() {
		HWallRaw *raw;
		for (int e=0; e<evts; e++) {
			if (e%5000==0) printf("Event %d\n",e);
			gLoop->nextEvent(e);
			for (int h=0; h<rawCat->getEntries(); h++) {
				raw=(HWallRaw*)HCategoryManager::getObject(raw,rawCat,h);
				int id = raw->getCell();
				if (id>=N_CELLS) continue;
				cells[id]->fill(raw->getWidth(1));
			}
		}

	}
	void saveCells() {
		ioF->cd();
		for (int id=0; id<N_CELLS; id++) {
			cells[id]->setCounts(cells[id]->getInit()->GetEntries());
			cells[id]->Write("",TObject::kOverwrite);
			cells_charge[id]->Write("",TObject::kOverwrite);
			cells_cosmic[id]->Write("",TObject::kOverwrite);
		}
	}
	void fitBe() {
		fillFrequencyHist();
		doFit(cells);
		for (int id=0; id<N_CELLS; id++) relPairAreaH->Fill(cells[id]->findBetterCombination());
		
	}
	void fitCo() {
		doFit(cells_cosmic);
		for (int id=0; id<N_CELLS; id++) relPairAreaH->Fill(cells[id]->findBetterCombination());
	}

	void doFit(WidthFitter *fitters[]) {
		for (int id=0; id<N_CELLS; id++) {
			printf("\n\n\n\n===========FITTING OF CELL %d ===============\n",id);			
			fitters[id]->flatInit();
			fitters[id]->subtractBg();
			fitters[id]->doSpectr();
			fitters[id]->prepareFitFuncs();
			fitters[id]->doFit();
			fitters[id]->constructTotFit();
			
		}
	}
	TList *getInitList(WidthFitter *fitters[]) {
		TList *res=new TList;
		int minId=0;
		if (fitters[0]->getType()==2) minId=208; 
		for (int id=minId; id<N_CELLS; id++) {
			TH1* h =fitters[id]->getInit();
			h->GetListOfFunctions()->Add(fitters[id]->getPM());		
			h->SetLineColor(kBlue);
			h->SetTitle(Form("CELL %d",id));
			res->Add( h );
		}
		return res;
	}
	TList *getApproxList(WidthFitter *fitters[]) {
		TList *res=new TList;
		int minId=0;
		if (fitters[0]->getType()==2) minId=208; 
		for (int id=minId; id<N_CELLS; id++) {
			TH1 *h=fitters[id]->getApprox();
			h->GetListOfFunctions()->Add(fitters[id]->getZ1Z2peaks()[0]);
			h->GetListOfFunctions()->Add(fitters[id]->getZ1Z2peaks()[1]);
			// h->SetLineColor(kRed);
			h->SetLineColor(4000);
			h->SetTitle(Form("CELL %d",id));
			res->Add( h );
		}
		return res;
	}
	TList *getBgList(WidthFitter *fitters[]) {
		TList *res = new TList;
		int minId=0;
		if (fitters[0]->getType()==2) minId=208; 
		for (int id=minId; id<N_CELLS; id++) {
			TH1 *h=fitters[id]->getBg();
			h->SetLineColor(kMagenta);
			h->SetTitle(Form("CELL %d",id));
			res->Add(h); 
		}
		return res;
	}
	TPolyMarker **getPMs(WidthFitter *fitters[]) {
		TPolyMarker **res = new TPolyMarker*[N_CELLS];
		int minId=0;
		if (fitters[0]->getType()==2) minId=208; 
		for (int id=minId; id<N_CELLS; id++) {
			res[id]= fitters[id]->getPM();
		}
		return res;		
	}

	TH1 *getApproxDiffHist(WidthFitter *fitters[]) {
		TH1 *res=new TH1F("approxDif","",200,0,50000);
		int minId=0;
		if (fitters[0]->getType()==2) minId=208; 
		for (int id=minId; id<N_CELLS; id++) res->Fill(fitters[id]->approxDiff());
		return res;
	}
	float *getAppropPeaksNum(WidthFitter *fitters[]) {
		float *res = new float[N_CELLS];
		int minId=0;
		if (fitters[0]->getType()==2) minId=208; 
		for (int id=minId; id<N_CELLS; id++) res[id]=fitters[id]->countAppropPeaks();
		return res;
	}

	void fillFrequencyHist() {
		for (int id=0; id<N_CELLS; id++) {
			int counts = cells[id]->getCounts();
			// printf("Counts in cell%d: %d\n\n", id, counts);
			freqHist->SetBinContent(id+1,counts);
			cells[id]->setFrequencyHist(freqHist);
		}
	}

	void fitQA_be(const char *pdfname) {
		fitQA(pdfname,cells);
		canva->Print(Form("%s]",pdfname));
	}
	void fitQA_co(const char *pdfname) {
		fitQA(pdfname,cells_cosmic);
		canva->Clear();
		getZ1ratioHist()->Draw();
		canva->Print(Form("%s)",pdfname));
	}

	void fitQA(const char *pdfname, WidthFitter *fitters[]) {
		TList *il = getInitList(fitters);
		TList *al = getApproxList(fitters);
		TList *bl = getBgList(fitters);
		TPolyMarker **pms = getPMs(fitters);
		TH1 *approxDiff = getApproxDiffHist(fitters);
		float *apprNum = getAppropPeaksNum(fitters);

		THStack st1, st2, st3;
		int pad=0;
		int hist= (fitters[0]->getType()==2) ? 208 : 0;
		canva->Print(Form("%s[",pdfname));


		for (int i=0; i<il->GetEntries(); i++) {
			st1.Add((TH1*)il->At(i));
			st2.Add((TH1*)al->At(i),"same");
			st3.Add((TH1*)bl->At(i),"same");
			if (++pad==N_PADS) {
				st1.Draw("PADS");
				st2.Draw("PADS,same");
				st3.Draw("PADS,same");
				canva->Update();
				
				
				for (int p=0; p<N_PADS; p++) {
					if (hist+p<208 && apprNum[hist+p]<2 || apprNum[hist+p]<1) canva->GetPad(p+1)->SetFillColor(RED);
					else  canva->GetPad(p+1)->SetFillColor(GREEN);
				}
				hist+=N_PADS;
				
				canva->Print(pdfname);
				pad=0;
				st1.GetHists()->Clear();
				st2.GetHists()->Clear();
				st3.GetHists()->Clear();
			}
		}

		if (st1.GetNhists()) {
			st1.Draw("PADS");
			st2.Draw("PADS,same");
			st3.Draw("PADS,same");
			canva->Update();
			for (int p=0; p<st1.GetNhists(); p++) {
				if (hist+p<208 && apprNum[hist+p]<2 || apprNum[hist+p]<1) canva->GetPad(p+1)->SetFillColor(RED);
				else  canva->GetPad(p+1)->SetFillColor(GREEN);
			}
			canva->Print(pdfname);
		}

		//distribution of integrated difference fit/init
		canva->Clear();
		approxDiff->Draw();
		canva->Print(pdfname);

		//distribution of relative area under Z1/Z2 peaks
		canva->Clear();
		relPairAreaH->Draw();
		canva->Print(pdfname);

		//disribution of each cell freq
		canva->Clear();
		freqHist->Draw();
		canva->Print(pdfname);

	}

	void checkCells_be() {
		checkCells(cells);
	}
	void checkCells_co() {
		checkCells(cells_cosmic);
	}

	void checkCells(WidthFitter *fitters[]) {
		for (int id=0; id<N_CELLS; id++) {
			printf("%s\n",string(122,'-').c_str());
			printf("|%9s %19s %17s   %17s   %17s   %17s   %17s   |\n",     "Cell","num of apprPeaks","diff","relPairArea","p1_sigma","p2_sigma","sigmaUpperLimit");
			double temp,sigmaUpper;
			fitters[id]->getZ1Z2peaks()[1]->GetParLimits(2,temp,sigmaUpper);
			printf("|%9d %19d %19.2f %19.2f %19.2f %19.2f %19.2f |\n",id,fitters[id]->countAppropPeaks(),fitters[id]->approxDiff(),fitters[id]->findBetterCombination(),
				fitters[id]->getZ1Z2peaks()[0]->GetParameter(2),fitters[id]->getZ1Z2peaks()[1]->GetParameter(2), sigmaUpper);
		
			//PEAKS
			vector<float> *positions = fitters[id]->getPeakPositions();
			printf("PEAKS: ");
			for (float &pos : *positions) printf("%10.2f",pos);
			printf("\n");
			////////
			// printf("counts: %d\t freqHist->max: %.2f\n",cells[id]->getCounts(),freqHist->GetMaximum());
			printf("INTENCITY:\t%.2f\n",1.0*fitters[id]->getCounts()/freqHist->GetMaximum());
			printf("%s\n",string(122,'-').c_str());
		}
	}

	void doCal() {
		for (int id=0; id<N_CELLS; id++) {
			float slope=0, offset=0;
			if (id<208) {
				if (cells[id]->countAppropPeaks()<2) continue;
				slope=100 / ( cells[id]->getZ1Z2peaks()[1]->GetParameter(1)-cells[id]->getZ1Z2peaks()[0]->GetParameter(1) );
				offset = 100 - slope*cells[id]->getZ1Z2peaks()[0]->GetParameter(1);
			} else {
				if (cells[id]->countAppropPeaks()<1 || cells_cosmic[id]->countAppropPeaks()<1) continue;
				slope=240 / ( cells_cosmic[id]->getZ1Z2peaks()[0]->GetParameter(1)-cells[id]->getZ1Z2peaks()[0]->GetParameter(1) );
				offset= 100 - slope*cells[id]->getZ1Z2peaks()[0]->GetParameter(1);
			}
			adc_slopes[id]=slope;
			adc_offsets[id]=offset;
			printf("---->doCal(): Cell%d. adc_slope: %.2f \t adc_offsets: %.2f \t apprBeam: %d \t apprCosm: %d\n",id,slope,offset, cells[id]->countAppropPeaks(),cells_cosmic[id]->countAppropPeaks());
		}
	}

	void loopBeRecal() {
		HWallHit *hit;
		HWallRaw *raw;
		for (int id=0; id<N_CELLS; id++) cells_charge[id]->getInit()->Reset();
		for (int e=0; e<evts; e++) {
			if (e%5000==0) printf("Event %d\n",e);
			gLoop->nextEvent(e);
			// for (int h=0; h<hitCat->getEntries(); h++) {
			for (int h=0; h<rawCat->getEntries(); h++) {
				// hit=(HWallHit*)HCategoryManager::getObject(hit,hitCat,h);
				raw=(HWallRaw*)HCategoryManager::getObject(raw,rawCat,h);
				int id = raw->getCell();
				if (id>=N_CELLS) continue;
				float val = raw->getWidth(1)*adc_slopes[id]+adc_offsets[id];
				// printf("---loop2(): filling with %.2f\n",val);
				cells_charge[id]->fill(val);
			}
		}		
	}

	TH2 *getQAplot() {
		TH2 *qaPlot = new TH2F("qaPlot", "", N_BINS, CHARGE_MIN,CHARGE_MAX, N_CELLS, 0, N_CELLS);
		for (int id=1; id<=N_CELLS; id++) {
			if (id<208 && cells[id]->countAppropPeaks()<2) continue;
			else if (cells[id]->countAppropPeaks()<1) continue;
			for (int bin=1; bin<=N_BINS;bin++) {
				 qaPlot->SetBinContent(bin,id,cells_charge[id-1]->getInit()->GetBinContent(bin));
			} 
		}
		float maxVal = qaPlot->GetMaximum();
		for (int id=1; id<=N_CELLS; id++) {
			if (id<208 && cells[id]->countAppropPeaks()>=2) continue;
			else if (cells[id]->countAppropPeaks()>=1) continue;
			for (int bin=1; bin<=N_BINS;bin++) {
				 qaPlot->SetBinContent(bin,id,maxVal);
			}
		}
		return qaPlot;
	}

	// hist with distribution of z1_cosm/z1_beam values for good fitted outter cells
	TH1 *getZ1ratioHist() {
		TH1 *res = new TH1F("z1Ratio", "", 25,0,5);
		for (int id=208; id<N_CELLS; id++) {
			if (cells[id]->countAppropPeaks()<1 || cells_cosmic[id]->countAppropPeaks()<1) continue;
			float z1Cosm = cells_cosmic[id]->getZ1Z2peaks()[1]->GetParameter(1);
			float z1Beam = cells[id]->getZ1Z2peaks()[0]->GetParameter(1);
			printf("----->Cell%d: z1Cosm / z1Beam = %.2f\n",id,z1Cosm/z1Beam);
			res->Fill(z1Cosm/z1Beam);
		}	
		return res;
	}

	void checkCharge() {
		printf("Addreses of charge initH: \n");
		for (int id=0; id<N_CELLS; id++) {
			printf("\tCell %d: %p\n",id,cells_charge[id]->getInit());
		}
	}
	void resetCo() {
		for (int id=0; id<N_CELLS; id++) {
			cells_cosmic[id]->getInit()->Reset();
		}
	}

	void loopCo() {
		HWallRaw *raw;
		for (int e=0; e<evts; e++) {
			if (e%5000==0) printf("Event %d\n",e);
			gLoop->nextEvent(e);
			for (int h=0; h<rawCat->getEntries(); h++) {
				raw=(HWallRaw*)HCategoryManager::getObject(raw,rawCat,h);
				int id = raw->getCell();
				if (id>=N_CELLS) continue;
				cells_cosmic[id]->fill(raw->getWidth(1));
			}
		}		
	}
};

TSpectrum Calibration::sp1(WidthFitter::N_PEAKS_1, WidthFitter::SP_RES_1);
TSpectrum Calibration::sp2(WidthFitter::N_PEAKS_2, WidthFitter::SP_RES_2);
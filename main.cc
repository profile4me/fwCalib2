#include "Calibration.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
//
const char *DST_NAME = "/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020101.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020102.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020103.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020104.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020105.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020106.root";
// const char *DST_NAME = "/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020101.root";
const char *DST_COSMIC_DAY57="/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905708102201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905708074701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700094501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700220801.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905708090601.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700423701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700034201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700281201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700155101.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700444601.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905708110301.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905708052401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700505501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700323701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/057/co1905700383801.root";
const char *DST_COSMIC_DAY56="/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905600343501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905600120601.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905600262801.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905600201501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905623170001.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905623105901.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905623230401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905623040701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905600035701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/056/co1905623100801.root";
const char *DST_COSMIC_DAY71="/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112210001.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112211101.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112400201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112472401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112185001.root";
const char *DST_COSMIC_DAY52="/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905223300801.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905217170201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905223212401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905217260501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905217081801.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905223121701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905217350001.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905223032601.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/052/co1905223390001.root";
const char *DST_FW="/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112210001.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112211101.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112400201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112472401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/071/fw1907112185001.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/067/fw1906710500501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/058/fw1905823192501.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/058/fw1905820201701.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/058/fw1905820163201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/058/fw1905823581801.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/058/fw1905822543401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/059/fw1905909184401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/059/fw1905901193901.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/059/fw1905909561201.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/059/fw1905908441401.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/059/fw1905909134601.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/059/fw1905900390101.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/066/fw1906615513901.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/066/fw1906615460501.root";


void init() {
	Calibration *cal=new Calibration();
	cal->saveCells();
	delete cal;	
}
void loopBe() {
	Calibration *cal=new Calibration("calibration.root");
	cal->initLoop(DST_NAME);
	cal->loopBe();
	cal->saveCells();
	delete cal;
}
void fitBe() {
	Calibration *cal=new Calibration("calibration.root");
	cal->fitBe();
	cal->checkCells_be();
	cal->fitQA_be("pdfs/fitQA.pdf");
	cal->saveCells();
	delete cal;
}

void calib_loopBeRecal() {
	Calibration *cal=new Calibration("calibration.root");
	cal->doCal();
	cal->initLoop(DST_NAME);
	cal->loopBeRecal();
	cal->saveCells();
	delete cal;
}

void recalQA() {
	Calibration *cal=new Calibration("calibration.root");
	// cal->checkCharge();
	TH2 *qaPlot = cal->getQAplot();
	TCanvas canva;
	qaPlot->Draw("colz");
	canva.SetLogz();
	canva.Print("pdfs/qaPlot.pdf)");
	TFile *outF=new TFile("test.root","recreate");
	qaPlot->Write();
	outF->Close();
}
void loopCo() {
	Calibration *cal=new Calibration("calibration.root");
	cal->initLoop(DST_FW);
	cal->resetCo();
	cal->loopCo();
	cal->saveCells();
	delete cal;
}
void fitCo() {
	Calibration *cal=new Calibration("calibration.root");
	cal->fitCo();
	cal->checkCells_co();
	cal->fitQA_co("pdfs/fitQA_cosm.pdf");
	cal->saveCells();
	delete cal;
}
int main(int argc, char const *argv[])
{	
	gROOT->SetBatch(1);
	TSystemDirectory dir("baseDir",".");
    TObject *obj = dir.GetListOfFiles()->FindObject("pdfs");
    if(!obj) gSystem->mkdir("pdfs");


    init();
	
	loopBe();
	fitBe();
	
	loopCo();
	fitCo();

	calib_loopBeRecal();
	recalQA();
	return 0;
}
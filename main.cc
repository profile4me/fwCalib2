#include "Calibration.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
//
const char *DST_NAME = "/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020101.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020102.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020103.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020104.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020105.root,/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020106.root";
// const char *DST_NAME = "/lustre/nyx/hades/user/dborisen/dstProduction/dst_mar19/063/be1906300020101.root";

void step1() {
	Calibration *cal=new Calibration();
	cal->initLoop(DST_NAME);
	cal->loop();
	cal->saveCells();
	delete cal;
}
void step2() {
	Calibration *cal=new Calibration("calibration.root");
	cal->doFit();
	cal->checkCells();
	cal->fitQA("pdfs/fitQA.pdf");
	cal->saveCells();
	delete cal;
}

void step3() {
	Calibration *cal=new Calibration("calibration.root");
	cal->doCal();
	cal->initLoop(DST_NAME);
	cal->loop2();
	cal->saveCells();
	delete cal;
}

void step4() {
	Calibration *cal=new Calibration("calibration.root");
	// cal->checkCharge();
	TH2 *qaPlot = cal->getQAplot();
	TCanvas canva;
	qaPlot->Draw("colz");
	canva.SetLogz();
	canva.Print("pdfs/qaPlot.pdf)");
}

int main(int argc, char const *argv[])
{	
	gROOT->SetBatch(1);
	TSystemDirectory dir("baseDir",".");
    TObject *obj = dir.GetListOfFiles()->FindObject("pdfs");
    if(!obj) gSystem->mkdir("pdfs");


	step1();
	step2();
	step3();
	step4();

	return 0;
}
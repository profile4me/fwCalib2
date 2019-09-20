#include "FitModule.h"

const int WidthFitter::N_PEAKS_1 = 5;
const int WidthFitter::N_PEAKS_2 = 2;
const float WidthFitter::SP_RES_1 = 1.0;
const float WidthFitter::SP_RES_2 = 1.0;
const float WidthFitter::SP_SIGMA_1 = 3.0;
const float WidthFitter::SP_SIGMA_2 = 15.0;
const float WidthFitter::SP_THRESH_1 = 0.005;
const float WidthFitter::SP_THRESH_2 = 0.025;
const int WidthFitter::N_ITERATIONS_FLATTER_1 = N_BINS/100;
const int WidthFitter::N_ITERATIONS_FLATTER_2 = N_BINS/40;

const float WidthFitter::HW_COEF = 0.5;
const float WidthFitter::FREQ_PERCENTAGE = 0.9;

ClassImp(WidthFitter);

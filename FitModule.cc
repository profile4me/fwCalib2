#include "FitModule.h"

const int WidthFitter::N_PEAKS = 5;
const float WidthFitter::HW_COEF = 0.5;
const float WidthFitter::SP_SIGMA = 3.0;
const float WidthFitter::SP_THRESH = 0.005;
const float WidthFitter::FREQ_PERCENTAGE = 0.9;
TSpectrum WidthFitter::sp(N_PEAKS);

ClassImp(WidthFitter);

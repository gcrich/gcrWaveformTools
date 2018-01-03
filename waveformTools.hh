#ifndef GCRWAVEFORMTOOLS_HH
#define GCRWAVEFORMTOOLS_HH

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef BATCH
#define BATCH 1
#endif



#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTemplWaveform.hh"
#include "SystemOfUnits.hh"
#include <deque>
//#include "/Users/grayson/Documents/quenchingFactors/tunlAnalysisTools/csiGlobalAnalysisParams.hh"
//
//#ifndef CSIGLOBALANALYSISPARAMS_H

// constants for the FIR filter
// could maybe make them as static consts rather than with defines
// either way, the goal is to globalize them and improve performance with
// compiletime setting/memory stuff
#define CSI_PEFIR_WIDTH 20
#define CSI_PEFIR_DELAY 15

#define BD_PEFIR_WIDTH 50
#define BD_PEFIR_DELAY 35







struct localIntegralResult_t {
    Double_t integral;
    Int_t integrationLength;
};



struct baseline_t {
    Double_t baseline;
    Double_t FWHM;
    Double_t chisquare;
};



static const Int_t trapFilter_scat_width = 5;
static const Int_t trapFilter_scat_delay = 3;
static const Int_t trapFilter_backingDet_width = 5;
static const Int_t trapFilter_backingDet_delay = 1;

// this "pretrace" is really just for SPE calibration
static const Int_t pretraceStart_index = 25;
static const Int_t pretraceEnd_index = 525;

static const Int_t posttraceStart_index = 8000;
static const Int_t posttraceEnd_index = 8500;

static const Double_t posttrace_start_time = 16500;
static const Double_t posttrace_end_time = 20200;

// this pretrace actually covers the region before we'd expect signals
// in other words, this is the one we probably want to use for pretrace\
//   signal cuts
static const Double_t longPretrace_start_time = 500.;
static const Double_t longPretrace_end_time = 4200.;


static const Double_t longPretrace_rollingBaseline_start_time = 500;
static const Double_t longPretrace_rollingBaseline_end_time = 4200;
static const Double_t pretrace_rollingBaseline_start_time = 500;
static const Double_t pretrace_rollingBaseline_end_time = 1500;


// time in ns relative to BD CFD time that we should START integrating scatter
static const Double_t scatIntegration_start_relToBDCFD_time = -800;
// time in ns relative to BD CFD time that we should STOP integrating scatter
static const Double_t scatIntegration_stop_relToBDCFD_time = 9200;


// in nanoseconds, with respect to where CFD occurs
static const Double_t LSCELLS_INTEGRATION_START_WRTCFD = -6;
// next, the END of the integration region for both the tail and full..
static const Double_t LSCELLS_INTEGRATION_END_WRTCFD = 200;
// finally, the start of the tail region
static const Double_t LSCELLS_TAIL_START_WRTCFD = 16;

//#endif

//
// this function produces the filtered and jump-rejecting baseline
// it populates the array passed into it as address rollingPwave
// assumes rollingPwave is already allocated and appropriate in size
// width sets the width of the averaging window
//
template<typename _Tp>
void getRollingPwave( TTemplWaveform<_Tp>* waveform, double* rollingPwave,
                     int halfWidth = 50, double preloadValue = -7777.7,
                     double rejectThreshold = 20);


/*
 getCMAfilteredWave
 
 perform CMA filtering - this returns a pointer to the waveform resulting from subtraction of the CMA-determined baseline
 inversion can be performed as well
 
 assumes that rollingPwave and cmaArray point to pre-allocated arrays with sufficient size to include the entire length of waveform
 */
template<typename _Tp>
TTemplWaveform<Double_t>* getCMAfilteredWave( TTemplWaveform<_Tp>* waveform,
                                             double* rollingPwave,
                                             double* cmaArray,
                                             int halfWidth = 50,
                                             double preloadValue = -7777.7,
                                             double rejectThreshold = 20,
                                             bool invert = true);


//
// this performs a triangle filter on the supplied waveform
// populates the array also supplied
//--- assumes the supplied array is properly allocated already
// order specifies how many samples are included
// --- order 3 includes 2 samples on either side of the index
// order must be >= 1
template<typename _Tp>
void getTriangleFilteredWave( TTemplWaveform<_Tp>* waveform,
                             double* triangleFilteredWave, int order = 3);




// this function is independent of what KIND of twaveform you've got (e.g. TDoubleWaveform, TShortWaveform)
template<typename _Tp>
void getIndicesFromTimeInWave( TTemplWaveform<_Tp>* wave, Double_t startTime, Double_t endTime, Int_t* indices );


template<typename _Tp>
Double_t getBaseline( TTemplWaveform<_Tp>* tWaveform, Double_t baselineGuess = 15600 );




// this function is effectively a way of evaluating the baseline subtracted waveform
// BUT it protects against out-of-range indexing
// so if the index you request is either too big (index >= wave->GetLength()) or too small (index < 0)
// then this function returns 0
// otherwise, it returns the baseline subtracted INVERTED value (ie use this for negative going signals)
// could be updated to handle negative & positive going signals, but just negative for now
template<typename _Tp>
Double_t getWaveformValueSafely( TTemplWaveform<_Tp>* wave, baseline_t baseline, Int_t index, bool invert = true );





Double_t getBaselineDouble( TTemplWaveform<Double_t>* tWaveform, Short_t baselineGuess = 15600 );


template<typename _Tp>
baseline_t getBaselineData( TTemplWaveform<_Tp>* tWaveform, Double_t baselineGuess = 14000, Int_t startIndex = 10, Int_t nBaselineCalcLength = 1000 );





template<typename _Tp>
Double_t* getBaselineFWHM( TTemplWaveform<_Tp>* tWaveform, Double_t baselineGuess = 15600 );




template<typename _Tp>
baseline_t baselineFit( TTemplWaveform<_Tp>* wave, Int_t baselineCalcLength = 14000, bool useHalfGaussian = false );







TTemplWaveform<Double_t>* getBaselineSuppressedWaveform( TTemplWaveform<Short_t>* tWaveform, baseline_t baselineData );




template<typename _Tp>
Double_t getPeakValue( TTemplWaveform<_Tp>* tWaveform, baseline_t baselineData, Double_t startTime = -.20, Double_t endTime = -2, bool invert = true );









Double_t getWaveformValueForTrapFilter( TTemplWaveform<Short_t>* wave, Int_t index, bool invert = true, Double_t baseline = -777.77 );




// this is effectively... some kind of .. zero crossing discrim?
// was originally called getCFDtime but seemed to work on trap filtered waveforms...
// or... something...
/*
 Double_t getCFDtimeTrapFilter( TTemplWaveform<Double_t>* tWaveform, Double_t threshold, Double_t fraction = 0.2, bool invert = true, Double_t searchStartTime = 0, Double_t searchEndTime = 0 ) {
 
 // to use getBaseline, we need a guess of the baseline
 // assume the average of the first 10 samples gives us a good guess for the baseline
 Double_t baselineGuess = 0;
 for( Int_t i = 0; i < 10; i++ ) {
 baselineGuess += tWaveform->At(i);
 }
 baselineGuess /= 10;
 
 Double_t baseline = getBaselineDouble( tWaveform, baselineGuess );
 
 Double_t maxValue = 0;
 Int_t index = 0;
 
 for( Int_t i = 1000; i < tWaveform->GetLength(); i++ ) {
 if( invert ) {
 if( -1 * (tWaveform->At(i) - baseline) > maxValue ) {
 maxValue = -1 * (tWaveform->At(i) - baseline);
 index = i;
 }
 }
 else {
 if( (tWaveform->At(i) - baseline) > maxValue ) {
 maxValue = (tWaveform->At(i) - baseline);
 index = i;
 }
 }
 
 }
 
 //    printf("found max of %.2f at index %i\n", maxValue, index );
 
 //    return maxValue;
 
 
 
 //    printf("peak value found %.2f, searched from index %i to index %i\n", maxValue, 0, tWaveform->GetLength());
 if( maxValue < threshold ) {
 return 0;
 }
 
 
 Double_t sample, prevSample;
 
 for( Int_t k = index - 20; k < index; k++ ) {
 sample = getWaveformValueForTrapFilter( tWaveform, k, invert, baseline );
 prevSample = getWaveformValueForTrapFilter( tWaveform, k-1, invert, baseline );
 if( sample > fraction * maxValue ) {
 //            return tWaveform->GetTimeAtIndex(k - 1) + (0.2 * maxValue - prevSample)/(sample - prevSample);
 //            return (k-1) + (0.2 * maxValue - prevSample)/(sample - prevSample);
 return (tWaveform->GetTimeAtIndex(k-1) + tWaveform->GetTimeAtIndex(k))/2;
 
 }
 }
 return 0;
 }
 */




/*
 getCFDtime
 
 gets CFD-like timing of a pulse in a given waveform
 finds time at which signal passes through a given fraction of its max amplitude
 
 search region should be specified  - if not, will default to entire waveform
 will only find 1st pulse
 fxn isn't really intended to do a pulse SEARCH but instead to do a timing refinement
 
 args
 wave - waveform to investigate
 baseline - baseline_t object representing waveform's baseline
 threshold (optional, default 0) - threshold for peak detection, fxn does no calc if max value in region is < threshold
 fraction (optional, default 0.2) - fraction at which to calc CFD time
 invert (optional, def. true) - specify positive going (false) vs negative going (true) signals
 searchStartTime (optional, default 0) - time in waveform (ns) at which to start search for peak to CFD
 searchEndTime (optional, default is end of waveform ) - time in waveform (ns) at which to END search
 */
Double_t getCFDtime( TTemplWaveform<Short_t>* wave, baseline_t baseline, Double_t threshold = 0, Double_t fraction = 0.2, bool invert = true, Double_t searchStartTime = 0, Double_t searchEndTime = -1 );


/*
 integrateWaveform
 
 returns integral of waveform over a specified range, handling baseline subtraction and inversion (if needed) as well as offering ability to 'suppress' inclusion of values within a certain range around baseline
 
 don't recommend using suppression
 
 args
    wave - pointer to wave to integrate
    baseline - wave baseline
    startTime - start time for integration (ns)
    endTime - end time for integration (ns)
    invert (optional, def. true) - negative-going signals (true) or positive-going signals (false)
    suppress (optional, def. false) - suppress around baseline
    suppressThreshold_baselineRelative (optional, def 0) - threshold, relative to baseline, for suppression
 */
template<typename _Tp>
Double_t integrateWaveform( TTemplWaveform<_Tp>* wave, baseline_t baseline, Double_t startTime, Double_t endTime, bool invert = true, bool suppress = false, Double_t suppressThreshold_baselineRelative = 0. );







Double_t integrateWaveformSuppress( TTemplWaveform<Short_t>* wave, baseline_t baseline, Double_t startTime = 0, Double_t endTime = -1, bool invert = true );







template<typename _Tp>
localIntegralResult_t integrateWaveformLocal( TTemplWaveform<_Tp>* wave, baseline_t baseline, Double_t triggerTime, Double_t preTime = 10., Double_t postTime = 20., bool invert = true );







// calculate PSD for BD event
// during calibration runs with these detectors at TUNL in Jan 2015
// used TAIL integral boundaries of CFD-trigger index +41 samples through end of full integral
// used full integral boundaries of CFD-trigger -4 through CFD-trigger + 150
// calculate PSD for BD event
Double_t getBackingDetectorPSD( TTemplWaveform<Short_t>* bdWaveform, baseline_t baselineInfo, Double_t triggerTime, bool invert = false, bool suppress = false, Double_t integralStartWRTcfd = LSCELLS_INTEGRATION_START_WRTCFD, Double_t tailStartWRTcfd = LSCELLS_TAIL_START_WRTCFD, Double_t integralEndWRTcfd = LSCELLS_INTEGRATION_END_WRTCFD );





Double_t getWaveformValueForTrapFilterRaw( TTemplWaveform<Short_t>* wave, Int_t index );






TTemplWaveform<Double_t>* getTrapFilteredWaveform( TTemplWaveform<Short_t>* wave,
                                                  Int_t sumWidth = 100,
                                                  Int_t delayWidth = 20,
                                                  std::deque<Double_t>* leadingWindow = NULL,
                                                  std::deque<Double_t>* trailingWindow = NULL,
                                                  Double_t* firResponse = NULL );

TTemplWaveform<Double_t>* getTrapFilteredEnergyWaveform( TTemplWaveform<Short_t>* wave,
                                                        Int_t sumWidth = 100,
                                                        Int_t delayWidth = 20,
                                                        std::deque<Double_t>* leadingWindow = NULL,
                                                        std::deque<Double_t>* trailingWindow = NULL,
                                                        Double_t* firResponse = NULL );





// returns a vector of pulse times in ns when FIR wave exceeds specified threshold
std::vector<Double_t> getPulseTimesFromFIR( TTemplWaveform<Double_t>* firWave, Double_t startTime = 0, Double_t endTime = 0, Double_t threshold = 45.0, bool positiveGoing = false );



// as of now (march 16 2016) this skips 40 ns after finding a pulse
// so it won't effectively discern between closely correlated pulses, but this kind of pileup would be tricky anyway
template<typename _Tp>
std::vector<Double_t> getPulseTimesFromEdgeCounting( TTemplWaveform<_Tp>* wave, baseline_t baseline, bool invert = true, Double_t threshold = 30, Double_t startTime = 0, Double_t endTime = -1, Double_t skipTime = 40. );




// returns time (in ns) of next BPM zero crossing, starting from startSearchTime, in ns
Double_t getNextBPMzeroCrossing( TTemplWaveform<Short_t>* bpmWaveform, Double_t startSearchTime = 5500, Double_t baseline = 8150 );


/*
 sharpEventCheck_integral
 
 inspects waveform for presence of a 'sharp' event, which is likely interactions in PMT glass
 
 compares integrals of two regions
    region 1: [startTime, startTime + regionLength)
    region 2: [startTime + regionLength, startTime + 2*regionLength )
 if ratio r2/r1 < ratioThreshold, return true
 else, return false
 */
template<typename _Tp>
bool sharpEventCheck_integral( TTemplWaveform<_Tp>* wave, baseline_t baseline, bool invert=true, Double_t startTime= 5000, Double_t regionLength=200, Double_t ratioThreshold =0.1 );


/*
 sharpEventCheck_aq
 
 inspects waveform for presence of a 'sharp' event, which is likely interactions in PMT glass
 
 looks at ratio of peak height to integral
 if a / q > threshold, return true
 else false
 */
template<typename _Tp>
bool sharpEventCheck_aq( TTemplWaveform<_Tp>* wave, baseline_t baseline, bool invert=true, Double_t startTime= 5000, Double_t regionLength=50, Double_t ratioThreshold =0.1 );


#endif
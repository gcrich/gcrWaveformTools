#ifndef GCRWAVEFORMTOOLS_CC 
#define GCRWAVEFORMTOOLS_CC

#define DEBUG 0
//#define DEBUGFIR 0
#define DEBUGMEMLEAK 1
#define BATCH 1

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTemplWaveform.hh"
#include "SystemOfUnits.hh"
#include <deque>

#ifndef BATCH
    #include "TCanvas.h"
#endif

#include "waveformTools.hh"







//==============================================================================
//==============================================================================


// this function is independent of what KIND of twaveform you've got (e.g. TDoubleWaveform, TShortWaveform)
template<typename _Tp>
void getIndicesFromTimeInWave( TTemplWaveform<_Tp>* wave, Double_t startTime, Double_t endTime, Int_t* indices ) {
    
    //Int_t* indices = new Int_t[2];
    
    if( startTime <= 0 ) {
        indices[0] = 1;
    }
    else if( startTime >= (wave->GetLength()-1) * wave->GetSamplingPeriod() ) {
        //indices[0] = wave->GetLength()-1;
        indices[0] = 1;
    }
    else
    indices[0] = wave->GetIndexAtTime( startTime );
    
    // assume negative end times means go to end of waveform
    if( endTime <= 0 ) {
        indices[1] = wave->GetLength() - 1;
    }
    else if( endTime >= (wave->GetLength() - 1)*wave->GetSamplingPeriod() ) {
        indices[1] = wave->GetLength() - 1;
    }
    else
    indices[1] = wave->GetIndexAtTime( endTime );
    
    
    if( DEBUG >= 2 )
    printf("getIndicesFromTimeInWave - returning indices of %i and %i\n", indices[0], indices[1] );
    //return indices;
}


template<typename _Tp>
Double_t getBaseline( TTemplWaveform<_Tp>* tWaveform, Double_t baselineGuess ) {
    
    // number of samples of the beginning of the waveform to use to calculate the baseline
    Int_t nBaselineCalcLength = 200;
    
    // number of bins around max to consider in calculation
    Int_t nCalcWidth = 10;
    Double_t weight = 0;
    Double_t baseline = 0;
    
    // number of bins to have in the histogram we're working with
    Int_t nCalcBins = 500;
    
    TH1D* baselineHist = new TH1D( "baselineHist", "baselineHist;ADC units;Counts",
                                  nCalcBins, baselineGuess - nCalcBins/2,
                                  baselineGuess + nCalcBins/2 - 1 );
    baselineHist->SetDirectory(0);
    
    // populate the histogram
    for( Int_t i = 0; i < nBaselineCalcLength; i++ ) {
        baselineHist->Fill( tWaveform->At(i) );
    }
    
    for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
        i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
        baseline += baselineHist->GetBinCenter(i) * baselineHist->GetBinContent(i);
        weight += baselineHist->GetBinContent(i);
    }
    baseline = baseline / weight;
    
    baselineHist->Delete();
    
    return baseline;
}




// this function is effectively a way of evaluating the baseline subtracted waveform
// BUT it protects against out-of-range indexing
// so if the index you request is either too big (index >= wave->GetLength()) or too small (index < 0)
// then this function returns 0
// otherwise, it returns the baseline subtracted INVERTED value (ie use this for negative going signals)
// could be updated to handle negative & positive going signals, but just negative for now
template<typename _Tp>
Double_t getWaveformValueSafely( TTemplWaveform<_Tp>* wave, baseline_t baseline, Int_t index, bool invert ) {
    
    
    if (index < 0) {
        return 0;
    }
    if( index >= wave->GetLength() ) {
        return 0;
    }
    
    if( !invert ) {
        return ((Double_t)wave->At(index) - baseline.baseline);
    }
    //    fail through - return the inverted, baselinesubbed value
    return -1*((Double_t)wave->At(index) - baseline.baseline);
}





Double_t getBaselineDouble( TTemplWaveform<Double_t>* tWaveform, Short_t baselineGuess ) {
    
    // number of samples of the beginning of the waveform to use to calculate the baseline
    Int_t nBaselineCalcLength = 200;
    
    // number of bins around max to consider in calculation
    Int_t nCalcWidth = 10;
    Double_t weight = 0;
    Double_t baseline = 0;
    
    // number of bins to have in the histogram we're working with
    Int_t nCalcBins = 500;
    
    TH1D* baselineHist = new TH1D( "baselineHist", "baselineHist;ADC units;Counts",
                                  nCalcBins, baselineGuess - nCalcBins/2,
                                  baselineGuess + nCalcBins/2 - 1 );
    baselineHist->SetDirectory(0);
    
    // populate the histogram
    for( Int_t i = 0; i < nBaselineCalcLength; i++ ) {
        baselineHist->Fill( tWaveform->At(i) );
    }
    
    for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
        i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
        baseline += baselineHist->GetBinCenter(i) * baselineHist->GetBinContent(i);
        weight += baselineHist->GetBinContent(i);
    }
    baseline = baseline / weight;
    
    baselineHist->Delete();
    
    return baseline;
}


template<typename _Tp>
baseline_t getBaselineData( TTemplWaveform<_Tp>* tWaveform, Double_t baselineGuess, Int_t startIndex, Int_t nBaselineCalcLength ) {
    
    
    baseline_t retBaseline;
    
    // number of samples  of the waveform to use to calculate the baseline
    //    Int_t nBaselineCalcLength = 10000;
    //    index at which to start the calculation
    //Int_t startIndex = 10;
    
    // number of bins around max to consider in calculation
    Int_t nCalcWidth = 200;
    Double_t weight = 0;
    Double_t baseline = 0;
    
    Double_t baselineFWHM = 0;
    Int_t lowerSide_HM_index = 0;
    Int_t upperSide_HM_index = 0;
    
    // number of bins to have in the histogram we're working with
    Int_t nCalcBins = 500;
    
    TH1D* baselineHist = new TH1D( "baselineHist", "baselineHist;ADC units;Counts",
                                  nCalcBins, baselineGuess - nCalcBins/2,
                                  baselineGuess + nCalcBins/2 - 1 );
    baselineHist->SetDirectory(0);
    
    // populate the histogram
    for( Int_t i = startIndex; i < startIndex + nBaselineCalcLength; i++ ) {
        baselineHist->Fill( tWaveform->At(i) );
    }
    
    for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
        i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
        baseline += baselineHist->GetBinCenter(i) * baselineHist->GetBinContent(i);
        weight += baselineHist->GetBinContent(i);
    }
    baseline = baseline / weight;
    
    retBaseline.baseline = baseline;
    
    //    see if rebinning gets better FWHM estimate
    //    apr 2016 - noticed that the FWHM is not very accurate when compared to
    //    an actual gaussian fit
    //    though the mean is pretty dead on
    baselineHist->Rebin();
    nCalcWidth /= 2;
    
    //    now get the FWHM
    Double_t maxBinContent = baselineHist->GetBinContent( baselineHist->GetMaximumBin() );
    for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
        i <= baselineHist->GetMaximumBin(); i++ ) {
        if( lowerSide_HM_index == 0 && baselineHist->GetBinContent(i) >= 0.5 * maxBinContent ) {
            //            set the lower side index to the PREVIOUS bin
            //            just to be conservative
            lowerSide_HM_index = i-1;
            break;
        }
    }
    
    for( Int_t i = baselineHist->GetMaximumBin();
        i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
        if( upperSide_HM_index == 0 && baselineHist->GetBinContent(i) <= 0.5 * maxBinContent ) {
            upperSide_HM_index = i;
            break;
        }
    }
    
    baselineFWHM = baselineHist->GetBinCenter(upperSide_HM_index) - baselineHist->GetBinCenter(lowerSide_HM_index);
    
    baselineHist->Delete();
    //
    //    baselineData[0] = baseline;
    //    baselineData[1] = baselineFWHM;
    
    retBaseline.FWHM = baselineFWHM;
    
    return retBaseline;
}





template<typename _Tp>
Double_t* getBaselineFWHM( TTemplWaveform<_Tp>* tWaveform, Double_t baselineGuess ) {
    
    // number of samples of the beginning of the waveform to use to calculate the baseline
    Int_t nBaselineCalcLength = 200;
    
    // number of bins around max to consider in calculation
    Int_t nCalcWidth = 30;
    Double_t weight = 0;
    Double_t baseline = 0;
    
    Double_t baselineFWHM = 0;
    Int_t lowerSide_HM_index = 0;
    Int_t upperSide_HM_index = 0;
    
    // number of bins to have in the histogram we're working with
    Int_t nCalcBins = 500;
    
    TH1D* baselineHist = new TH1D( "baselineHist", "baselineHist;ADC units;Counts",
                                  nCalcBins, baselineGuess - nCalcBins/2,
                                  baselineGuess + nCalcBins/2 - 1 );
    baselineHist->SetDirectory(0);
    
    // populate the histogram
    for( Int_t i = 0; i < nBaselineCalcLength; i++ ) {
        baselineHist->Fill( tWaveform->At(i) );
    }
    
    for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
        i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
        baseline += baselineHist->GetBinCenter(i) * baselineHist->GetBinContent(i);
        weight += baselineHist->GetBinContent(i);
    }
    baseline = baseline / weight;
    
    //    now get the FWHM
    Double_t maxBinContent = baselineHist->GetBinContent( baselineHist->GetMaximumBin() );
    for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
        i <= baselineHist->GetMaximumBin(); i++ ) {
        if( lowerSide_HM_index == 0 && baselineHist->GetBinContent(i) >= 0.5 * maxBinContent ) {
            //            set the lower side index to the PREVIOUS bin
            //            just to be conservative
            lowerSide_HM_index = i-1;
            break;
        }
    }
    
    for( Int_t i = baselineHist->GetMaximumBin();
        i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
        if( upperSide_HM_index == 0 && baselineHist->GetBinContent(i) <= 0.5 * maxBinContent ) {
            upperSide_HM_index = i;
            break;
        }
    }
    
    baselineFWHM = upperSide_HM_index - lowerSide_HM_index;
    
    baselineHist->Delete();
    
    
    Double_t* returnValue = new Double_t[2];
    returnValue[0] = baseline;
    returnValue[1] = baselineFWHM;
    
    return returnValue;
}




template<typename _Tp>
baseline_t baselineFit( TTemplWaveform<_Tp>* wave, Int_t baselineCalcLength, bool useHalfGaussian ) {
    
    char buff[256];
    
    
    Double_t baselineGuess = wave->At(0);
    
    
    // number of bins around max to consider in calculation
    Int_t nCalcWidth = 10;
    Double_t weight = 0;
    Double_t baseline = 0;
    
    // number of bins to have in the histogram we're working with
    Int_t nCalcBins = 500;
    
    TH1D* baselineHist = new TH1D( "baselineHist", "baselineHist;ADC units;Counts",
                                  nCalcBins, baselineGuess - nCalcBins/2,
                                  baselineGuess + nCalcBins/2 - 1 );
    baselineHist->SetDirectory(0);
    
    // populate the histogram
    for( Int_t i = 0; i < baselineCalcLength; i++ ) {
        baselineHist->Fill( wave->At(i) );
    }
    
    baselineHist->Rebin(4);
    
    
    baselineGuess = baselineHist->GetBinCenter( baselineHist->GetMaximumBin() );
    
    TF1* baselineFit;
    if( useHalfGaussian ) {
        baselineFit = new TF1("baselineFit", "gaus(0)", baselineGuess - 5, baselineGuess + 75 );
    }
    else {
        baselineFit = new TF1("baselineFit", "gaus(0)", baselineGuess - 75, baselineGuess + 75 );
    }
    
    
    baselineHist->Fit( baselineFit, "LMREQN" );

    
    
    baseline_t retBaseline;
    retBaseline.baseline = baselineFit->GetParameter(1);
    retBaseline.FWHM = baselineFit->GetParameter(2);
    retBaseline.chisquare =baselineFit->GetChisquare();
    
    baselineHist->Delete();
    baselineFit->Delete();
    
    return retBaseline;
}







TTemplWaveform<Double_t>* getBaselineSuppressedWaveform( TTemplWaveform<Short_t>* tWaveform, baseline_t baselineData ) {
    /*
     // number of samples of the beginning of the waveform to use to calculate the baseline
     Int_t nBaselineCalcLength = 200;
     
     // number of bins around max to consider in calculation
     Int_t nCalcWidth = 30;
     Double_t weight = 0;
     Double_t baseline = 0;
     
     Double_t baselineFWHM = 0;
     Int_t lowerSide_HM_index = 0;
     Int_t upperSide_HM_index = 0;
     
     // number of bins to have in the histogram we're working with
     Int_t nCalcBins = 500;
     
     Double_t baselineGuess = tWaveform->At(0);
     
     TH1D* baselineHist = new TH1D( "baselineHist", "baselineHist;ADC units;Counts",
     nCalcBins, baselineGuess - nCalcBins/2,
     baselineGuess + nCalcBins/2 - 1 );
     baselineHist->SetDirectory(0);
     
     // populate the histogram
     for( Int_t i = 0; i < nBaselineCalcLength; i++ ) {
     baselineHist->Fill( tWaveform->At(i) );
     }
     
     for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
     i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
     baseline += baselineHist->GetBinCenter(i) * baselineHist->GetBinContent(i);
     weight += baselineHist->GetBinContent(i);
     }
     baseline = baseline / weight;
     
     //    now get the FWHM
     Double_t maxBinContent = baselineHist->GetBinContent( baselineHist->GetMaximumBin() );
     for( Int_t i = baselineHist->GetMaximumBin() - 1 - nCalcWidth;
     i <= baselineHist->GetMaximumBin(); i++ ) {
     if( lowerSide_HM_index == 0 && baselineHist->GetBinContent(i) >= 0.5 * maxBinContent ) {
     //            set the lower side index to the PREVIOUS bin
     //            just to be conservative
     lowerSide_HM_index = i-1;
     break;
     }
     }
     
     for( Int_t i = baselineHist->GetMaximumBin();
     i <= baselineHist->GetMaximumBin() + nCalcWidth - 1; i++ ) {
     if( upperSide_HM_index == 0 && baselineHist->GetBinContent(i) <= 0.5 * maxBinContent ) {
     upperSide_HM_index = i;
     break;
     }
     }
     
     baselineFWHM = upperSide_HM_index - lowerSide_HM_index;
     
     baselineHist->Delete(); */
    
    //    Double_t upperbound = baseline + baselineFWHM / 2.355 * 7;
    //    Double_t lowerbound = baseline - baselineFWHM / 2.355 * 7;
    
    
    Double_t upperbound = baselineData.baseline + baselineData.FWHM / 2.355 * 5;
    Double_t lowerbound = baselineData.baseline - baselineData.FWHM / 2.355 * 5;
    
    Double_t* newWaveArray = new Double_t[tWaveform->GetLength()];
    for( Int_t i = 0; i < tWaveform->GetLength(); i++ ) {
        //        if( tWaveform->At(i) > lowerbound && tWaveform->At(i) < upperbound ) {
        if( tWaveform->At(i) > lowerbound ) {
            newWaveArray[i] = 0;
        }
        else {
            newWaveArray[i] = -1 * (tWaveform->At(i) - baselineData.baseline );
        }
    }
    
    TTemplWaveform<Double_t>* baselineSuppressedWave = new TTemplWaveform<Double_t>( newWaveArray, tWaveform->GetLength() );
    baselineSuppressedWave->SetSamplingFreq( 500 * CLHEP::megahertz );
    
    delete [] newWaveArray;
    
    return baselineSuppressedWave;
    
}




template<typename _Tp>
Double_t getPeakValue( TTemplWaveform<_Tp>* tWaveform, baseline_t baselineData, Double_t startTime, Double_t endTime, bool invert ) {
    
    Int_t indices[2];
    getIndicesFromTimeInWave( tWaveform, startTime, endTime, indices );
    
    Int_t startIndex = indices[0], endIndex = indices[1];
    
    
//    Double_t baselineGuess = 0;
//    for( Int_t i = 0; i < 10; i++ ) {
//        baselineGuess += tWaveform->At(i);
//    }
//    baselineGuess /= 10;
//    
//    Double_t baseline = getBaseline( tWaveform, baselineGuess );
    
    Double_t maxValue = 0;
    Int_t index = 0;
    
    for( Int_t i = startIndex; i < endIndex; i++ ) {
        if( invert ) {
            if( -1 * (tWaveform->At(i) - baselineData.baseline) > maxValue ) {
                maxValue = -1 * (tWaveform->At(i) - baselineData.baseline);
                index = i;
            }
        }
        else {
            if( (tWaveform->At(i) - baselineData.baseline) > maxValue ) {
                maxValue = (tWaveform->At(i) - baselineData.baseline);
                index = i;
            }
        }
        
    }
    
    //    printf("found max of %.2f at index %i\n", maxValue, index );
    
    return maxValue;
}



/*
 // returns max value found over specified region, no baseline subtraction
 Double_t getPeakValueDouble( TTemplWaveform<Double_t>* tWaveform, Int_t startIndex = 0, Int_t endIndex = 0) {
 
 // to use getBaseline, we need a guess of the baseline
 // assume the average of the first 10 samples gives us a good guess for the baseline
 if( startIndex < 0 ) {
 startIndex = 0;
 }
 if( startIndex >= tWaveform->GetLength() ) {
 startIndex = tWaveform->GetLength() - 1; // trivial - means we do nothing
 }
 if( endIndex == 0 ) {
 endIndex = tWaveform->GetLength();
 }
 if( endIndex >= tWaveform->GetLength() ) {
 endIndex = tWaveform->GetLength()-1;
 }
 
 
 
 Double_t maxValue = 0;
 Int_t index = 0;
 
 for( Int_t i = startIndex; i < endIndex; i++ ) {
	if( tWaveform->At(i) > maxValue ) {
 maxValue = tWaveform->At(i);
 index = i;
	}
 
 }
 
 //    printf("found max of %.2f at index %i\n", maxValue, index );
 
 return maxValue;
 }
 */






Double_t getWaveformValueForTrapFilter( TTemplWaveform<Short_t>* wave, Int_t index, bool invert, Double_t baseline ) {
    if (index < 0) {
        return 0;
    }
    if( index >= wave->GetLength() ) {
        return 0;
    }
    
    if( baseline == -777.77 ) {
        baseline = getBaseline( wave );
    }
    
    if( invert )
    return -1*(wave->At(index)-baseline);
    return ( wave->At(index) - baseline );
}




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
Double_t getCFDtime( TTemplWaveform<Short_t>* wave, baseline_t baseline, Double_t threshold, Double_t fraction, bool invert, Double_t searchStartTime, Double_t searchEndTime ) {
    
    
    Double_t maxValue = 0;
    Int_t maxIndex = -1;
    
    Int_t searchStartIndex, searchEndIndex;
    Int_t indices[2];
    getIndicesFromTimeInWave( wave, searchStartTime, searchEndTime, indices );
    searchStartIndex = indices[0];
    searchEndIndex = indices[1];
    
    double prevValue = getWaveformValueSafely( wave, baseline, searchStartIndex-1, invert);
    Double_t value;
    for( int i = searchStartIndex; i <= searchEndIndex; i++ ) {
        value = getWaveformValueSafely( wave, baseline, i, invert );
        if( value > maxValue ) {
            maxValue = value;
            maxIndex = i;
        }
        if( value < prevValue  && maxValue >= threshold){
            break;
        }
    }
    
    
    if( DEBUG > 2 ) {
        printf( "getCFDtime -  found max value of %.2f at index %i or time %.0f\n", maxValue, maxIndex, wave->GetTimeAtIndex(maxIndex) );
    }
    
    /*
     now we've found the max value and the index at which it occurs
     step back.. 30 samples and search for CFD trigger
     */
    Double_t cfdTime= -777.77;
    if( maxValue < threshold ) {
        return cfdTime;
    }
    
    
    Double_t sampleValuePrevious = getWaveformValueSafely( wave, baseline, maxIndex - 30, invert );
    for( int i = maxIndex - 29; i <= maxIndex; i++ ) {
        value = getWaveformValueSafely( wave, baseline, i, invert );
        
        if( value > fraction * maxValue ) {
            cfdTime = ( i - 1 ) + ( fraction*maxValue - sampleValuePrevious)/(value - sampleValuePrevious);
            break;
        }
        sampleValuePrevious = value;
    }
    
    if(DEBUG>0 && cfdTime < 0) {
        printf("getCFDtime - WARNING - time < 0 found\n");
    }
    Double_t retVal = cfdTime * wave->GetSamplingPeriod();
    if( std::isinf(retVal) ) {
        if(DEBUG >= 1 ) {
            printf( "getCFDtime - WARNING - infinite CFD time found\n");
        }
        retVal = wave->GetTimeAtIndex(maxIndex);
    }
    return retVal;
}



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
Double_t integrateWaveform( TTemplWaveform<_Tp>* wave, baseline_t baseline, Double_t startTime, Double_t endTime, bool invert, bool suppress, Double_t suppressThreshold_baselineRelative ) {
    
    Int_t indices[2];
    getIndicesFromTimeInWave( wave, startTime, endTime, indices );
    Int_t startIndex = indices[0], endIndex = indices[1];
    
    
    
    Double_t integral = 0;
    Double_t sample;
    
    for( Int_t index = startIndex; index <= endIndex; index++ ) {
        
        sample = (Double_t)wave->At(index);
        
        // if we are suppressing things "under" the baseline...
        if( suppress ) {
            
            // if INVERT (negative going signals), then skip sample if it is HIGHER than our baseline-relative threshold
            if( invert && sample > (baseline.baseline - suppressThreshold_baselineRelative ) ) {
                continue;
            }
            
            // if NOT invert (positive going signals), skip sample if it is BELOW our baseline
            if( !invert && sample < (baseline.baseline + suppressThreshold_baselineRelative) ) {
                continue;
            }
        }
        
        // now increment the integral
        // we won't reach this if suppression is enabled and the sample was judged to be suppression-worthy
        if( invert ) {
            integral += -1 * (sample - baseline.baseline);
        }
        else {
            integral += sample - baseline.baseline;
        }
    }
    
    return integral;
}







Double_t integrateWaveformSuppress( TTemplWaveform<Short_t>* wave, baseline_t baseline, Double_t startTime, Double_t endTime, bool invert ) {
    
    //    if( startIndex < 0 ) {
    //        startIndex = 0;
    //    }
    //    if( endIndex > wave->GetLength() - 1 ) {
    //        endIndex = wave->GetLength() - 1;
    //    }
    
    Int_t indices[2];
    getIndicesFromTimeInWave( wave, startTime, endTime, indices );
    Int_t startIndex = indices[0], endIndex = indices[1];
    
    
    
    //    define cut region within +/- 3 sigma of baseline
    //    if a sample is within that region, call it 0
    Double_t upperbound = baseline.baseline + baseline.FWHM / 2.355 * 5;
    Double_t lowerbound = baseline.baseline - baseline.FWHM / 2.355 * 5;
    //    baseline = baseline.baseline;
    //    printf( "baseline %.2f, upper / lower bounds are  %.2f and %.2f\n", baselineInfo[0], upperbound, lowerbound );
    
    Double_t integral = 0;
    Double_t sample = 0;
    for( Int_t index = startIndex; index <= endIndex; index++ ) {
        sample = wave->At(index);
        //        printf("sample %.2f\n", sample );
        if( sample < lowerbound ) {
            integral += -1 * (wave->At(index) - baseline.baseline);
            //            printf("integral now %.2f after sample %i\n", integral, index);
        }
    }
    //    delete [] baselineInfo;
    //    printf("integral is %.2f\n", integral);
    return integral;
}







template<typename _Tp>
localIntegralResult_t integrateWaveformLocal( TTemplWaveform<_Tp>* wave, baseline_t baseline, Double_t triggerTime, Double_t preTime, Double_t postTime, bool invert ) {
    
    localIntegralResult_t result;
    result.integral =  0.;
    result.integrationLength = 0;
    Double_t integral = 0.;
    
    //    if( baseline < 0 ) {
    ////        _Tp guess = (_Tp)wave->At(5);
    //        baseline = getBaseline( wave );
    //    }
    //    make a double for the baseline info if it wasn't supplied
    //    keep track of whether or not we need to clean up after ourselves
    
    
    Int_t countdown = (Int_t) (postTime / wave->GetSamplingPeriod());
    
    Int_t index = wave->GetIndexAtTime( triggerTime - preTime );
    while( index < wave->GetLength() ) {
        
        if( invert ) {
            integral += -1 * ( wave->At(index) - baseline.baseline );
        }
        else {
            integral += wave->At(index) - baseline.baseline;
        }
        
        if( TMath::Abs(wave->At(index) - baseline.baseline ) < 5 * baseline.FWHM ) {
            //            if the signal goes under the 5FWHM threshold, count down
            countdown--;
        }
        else {
            //            reset countdown if we're over threshold
            countdown = (Int_t) (postTime / wave->GetSamplingPeriod());
        }
        
        if( countdown <= 0 ) {
            break;
        }
        result.integrationLength++;
        index++;
    }
    
    
    //    return integral;
    result.integral = integral;
    return result;
}







// calculate PSD for BD event
// during calibration runs with these detectors at TUNL in Jan 2015
// used TAIL integral boundaries of CFD-trigger index +41 samples through end of full integral
// used full integral boundaries of CFD-trigger -4 through CFD-trigger + 150
// calculate PSD for BD event
Double_t getBackingDetectorPSD( TTemplWaveform<Short_t>* bdWaveform, baseline_t baselineInfo, Double_t triggerTime, bool invert, bool suppress, Double_t integralStartWRTcfd, Double_t tailStartWRTcfd, Double_t integralEndWRTcfd ) {
    
    //    baseline = getBaseline( bdWaveform, baseline );
    //    baseline_t baselineInfo = getBaselineData( bdWaveform, bdWaveform->At(10), 10, 5000 );
    
    //    convert trigger time in NS to samples
    //    Int_t fullIntegralStart = bdWaveform->GetIndexAtTime( triggerTime + LSCELLS_INTEGRATION_START_WRTCFD );
    Double_t fullIntegralStart_time = triggerTime + integralStartWRTcfd;
    //    Int_t fullIntegralEnd = bdWaveform->GetIndexAtTime( triggerTime + LSCELLS_INTEGRATION_END_WRTCFD );
    Double_t fullIntegralEnd_time = triggerTime + integralEndWRTcfd;
    //    Int_t shortIntegralStart = bdWaveform->GetIndexAtTime( triggerTime + LSCELLS_TAIL_START_WRTCFD );
    Double_t shortIntegralStart_time = triggerTime + tailStartWRTcfd;
    
    //    Double_t fullIntegralStart = triggerTime - 20;
    //    Double_t fullIntegralEnd = triggerTime + 100;
    //    Double_t shortIntegralStart = triggerTime + 10;
    Double_t psd = -100;
    if( !suppress ) {
        psd = integrateWaveform( bdWaveform, baselineInfo, shortIntegralStart_time, fullIntegralEnd_time, false ) / integrateWaveform( bdWaveform, baselineInfo, fullIntegralStart_time, fullIntegralEnd_time, false );
    }
    if( suppress ) {
        psd = integrateWaveform( bdWaveform, baselineInfo, shortIntegralStart_time, fullIntegralEnd_time, false, true ) / integrateWaveform( bdWaveform, baselineInfo, fullIntegralStart_time, fullIntegralEnd_time, false, true );
    }
    return psd;
}





Double_t getWaveformValueForTrapFilterRaw( TTemplWaveform<Short_t>* wave, Int_t index ) {
    if (index < 0) {
        return 0;
    }
    if( index >= wave->GetLength() ) {
        return 0;
    }
    
    return (Double_t)wave->At(index);
}





/*
TTemplWaveform<Double_t>* getTrapFilteredWaveform( TTemplWaveform<Short_t>* wave,
                                                  Int_t sumWidth,
                                                  Int_t delayWidth,
                                                  std::deque<Double_t>* leadingWindow,
                                                  std::deque<Double_t>* trailingWindow,
                                                  Double_t* firResponse ) {
    
    bool allocatedFIRresponse = false;
    if( firResponse == NULL ) {
        firResponse = new Double_t[wave->GetLength()];
        allocatedFIRresponse = true;
    }
    //    Double_t* firResponse = new Double_t[wave->GetLength()];
    
    Double_t beginSum = 0;
    Double_t beginDelay = 0;
    
    
    
    Double_t leadingSum = 0;
    Double_t trailingSum = 0;
    Double_t temp;
    
    //    populate the initial values for the running sums
    for( Int_t i = 0; i < sumWidth; i++ ) {
        temp = getWaveformValueForTrapFilterRaw( wave, i );
        leadingSum += temp;
        leadingWindow->push_back( temp );
        leadingWindow->pop_front();
    }
    for( Int_t i = -1*delayWidth; i <= sumWidth - delayWidth; i++ ) {
        temp = getWaveformValueForTrapFilterRaw( wave, i );
        trailingSum += temp;
        trailingWindow->push_back(temp);
        trailingWindow->pop_front();
    }
    
    //    print initial sums
    //    only if debugging...
    //    printf("\ninitial leading window\n");
    //    for( int dbi = 0; dbi < leadingWindow->size(); dbi++ ) {
    //        printf("%.2f\t", leadingWindow->at(dbi));
    //    }
    //    printf("\nintial trailing window\n");
    //    for( int dbi = 0; dbi < trailingWindow->size(); dbi++ ) {
    //        printf("%.2f\t", trailingWindow->at(dbi));
    //    }
    
    for( Int_t index = 0; index < wave->GetLength(); index++ ) {
        
        //        get rid of old values
        leadingSum -= leadingWindow->front();
        leadingWindow->pop_front();
        
        trailingSum -= trailingWindow->front();
        trailingWindow->pop_front();
        
        
        //        get new values
        temp = getWaveformValueForTrapFilterRaw( wave, index + sumWidth -1 );
        leadingSum += temp;
        leadingWindow->push_back(temp);
        
        temp = getWaveformValueForTrapFilterRaw( wave, index + sumWidth - delayWidth -1 );
        trailingSum += temp;
        trailingWindow->push_back(temp);
        
        //        compute FIR response at index
        firResponse[index] = leadingSum - trailingSum;
        
    }
    
    
    
    TTemplWaveform<Double_t>* firWave = new TTemplWaveform<Double_t>( firResponse, wave->GetLength() );
    firWave->SetSamplingFreq( 500 * CLHEP::megahertz );
    
    if( allocatedFIRresponse ) {
        //        clean up after ourselves if we made the array we used...
        delete [] firResponse;
    }
    
    return firWave;
}
*/





/*
 getTrapFilteredWaveform
 
 produces and returns a trapezoidal-filtered waveform from supplied waveform
 uses and algorithm very much like that implemented on Struck digitizers
 knobs to turn are sumWidth and delayWidth
 specify length of two rolling integration windows and their starting-index-offset, respectively
 uses deque's to do the calculation
 MUCH faster if they're globally allocated (ie define them once in your analysis loop and pass the pointers into this function, don't let this fxn allocate new ones every time it runs)
 also uses a scratch space - a double array with length equal to waveform
 same story about speed - faster if scratch space is allocated outside of this fxn
 
 args
 wave - waveform to trap filter
 sumWidth - width of the two rolling integrals to use for trap filter calc (in samples)
 delayWidth - offset between leading edges of rolling integrals (in samples)
 leadingWindow (optional) - deque used for calculation
 trailingWindow (optional) - deque used for calculation
 firResponse (optional) - scratch space
 */
TTemplWaveform<Double_t>* getTrapFilteredWaveform( TTemplWaveform<Short_t>* wave,
                                                  Int_t sumWidth,
                                                  Int_t delayWidth,
                                                  std::deque<Double_t>* leadingWindow,
                                                  std::deque<Double_t>* trailingWindow,
                                                  Double_t* firResponse ) {
    
    bool allocatedFIRresponse = false;
    if( firResponse == NULL ) {
#ifdef DEBUGMEMLEAK
        printf("allocated FIR response array\n");
#endif
        firResponse = new Double_t[wave->GetLength()];
        allocatedFIRresponse = true;
    }
    
    Double_t beginSum = 0;
    Double_t beginDelay = 0;
    
    
    
    Double_t leadingSum = 0;
    Double_t trailingSum = 0;
    Double_t temp;
    
    //    populate the initial values for the running sums
    for( Int_t i = 0; i < sumWidth; i++ ) {
        temp = getWaveformValueForTrapFilterRaw( wave, i );
        leadingSum += temp;
        leadingWindow->push_back( temp );
        leadingWindow->pop_front();
    }
    for( Int_t i = -1*delayWidth; i <= sumWidth - delayWidth; i++ ) {
        temp = getWaveformValueForTrapFilterRaw( wave, i );
        trailingSum += temp;
        trailingWindow->push_back(temp);
        trailingWindow->pop_front();
    }
    
    /*
     print initial sums
     only if debugging...
     */
#ifdef DEBUGFIR
     printf("\ninitial leading window\n");
     for( int dbi = 0; dbi < leadingWindow->size(); dbi++ ) {
     printf("%.2f\t", leadingWindow->at(dbi));
     }
     printf("\nintial trailing window\n");
     for( int dbi = 0; dbi < trailingWindow->size(); dbi++ ) {
     printf("%.2f\t", trailingWindow->at(dbi));
     }
#endif
    
    
    for( Int_t index = 0; index < wave->GetLength(); index++ ) {
        
        //        get rid of old values
        leadingSum -= leadingWindow->front();
        leadingWindow->pop_front();
        
        trailingSum -= trailingWindow->front();
        trailingWindow->pop_front();
        
        
        //        get new values
        temp = getWaveformValueForTrapFilterRaw( wave, index + sumWidth -1 );
        leadingSum += temp;
        leadingWindow->push_back(temp);
        
        temp = getWaveformValueForTrapFilterRaw( wave, index + sumWidth - delayWidth -1 );
        trailingSum += temp;
        trailingWindow->push_back(temp);
        
        //        compute FIR response at index
        firResponse[index] = leadingSum - trailingSum;
        
    }
    
    
    
    TTemplWaveform<Double_t>* firWave = new TTemplWaveform<Double_t>( firResponse, wave->GetLength() );
    firWave->SetSamplingFreq( 500 * CLHEP::megahertz );
    
    if( allocatedFIRresponse ) {
#ifdef DEBUGMEMLEAK
        printf("deleted fir response \n");
#endif
        //        clean up after ourselves if we made the array we used...
        delete [] firResponse;
    }
    
    return firWave;
}





/*
 getTrapFilteredEnergyWaveform
 
 not sure exactly what this does - can't remember off the top of my head 
 */
TTemplWaveform<Double_t>* getTrapFilteredEnergyWaveform( TTemplWaveform<Short_t>* wave,
                                                        Int_t sumWidth,
                                                        Int_t delayWidth,
                                                        std::deque<Double_t>* leadingWindow ,
                                                        std::deque<Double_t>* trailingWindow,
                                                        Double_t* firResponse ) {
    
    bool allocatedFIRresponse = false;
    if( firResponse == NULL ) {
        firResponse = new Double_t[wave->GetLength()];
        allocatedFIRresponse = true;
    }
    //    Double_t* firResponse = new Double_t[wave->GetLength()];
    
    Double_t beginSum = 0;
    Double_t beginDelay = 0;
    
    
    
    Double_t leadingSum = 0;
    Double_t trailingSum = 0;
    Double_t temp;
    
    //    populate the initial values for the running sums
    for( Int_t i = sumWidth+delayWidth; i < 2*sumWidth+delayWidth; i++ ) {
        temp = getWaveformValueForTrapFilterRaw( wave, i );
        leadingSum += temp;
        leadingWindow->push_back( temp );
        leadingWindow->pop_front();
    }
    for( Int_t i = 0; i <= sumWidth; i++ ) {
        temp = getWaveformValueForTrapFilterRaw( wave, i );
        trailingSum += temp;
        trailingWindow->push_back(temp);
        trailingWindow->pop_front();
    }
    
    
    for( Int_t index = 0; index < wave->GetLength(); index++ ) {
        
        //        get rid of old values
        leadingSum -= leadingWindow->front();
        leadingWindow->pop_front();
        
        trailingSum -= trailingWindow->front();
        trailingWindow->pop_front();
        
        
        //        get new values
        temp = getWaveformValueForTrapFilterRaw( wave, index + 2*sumWidth +delayWidth );
        leadingSum += temp;
        leadingWindow->push_back(temp);
        
        temp = getWaveformValueForTrapFilterRaw( wave, index + sumWidth );
        trailingSum += temp;
        trailingWindow->push_back(temp);
        
        //        compute FIR response at index
        firResponse[index] = leadingSum - trailingSum;
        
    }
    
    
    
    TTemplWaveform<Double_t>* firWave = new TTemplWaveform<Double_t>( firResponse, wave->GetLength() );
    firWave->SetSamplingFreq( 500 * CLHEP::megahertz );
    
    if( allocatedFIRresponse ) {
        //        clean up after ourselves if we made the array we used...
        delete [] firResponse;
    }
    
    return firWave;
}





// returns a vector of pulse times in ns when FIR wave exceeds specified threshold
std::vector<Double_t> getPulseTimesFromFIR( TTemplWaveform<Double_t>* firWave, Double_t startTime, Double_t endTime, Double_t threshold, bool positiveGoing ) {
    
    std::vector<Double_t> pulseTimes;
    
    bool okToCount = true;
    
    
    Int_t startSample = firWave->GetIndexAtTime( startTime );
    
    Int_t endSample = firWave->GetIndexAtTime( endTime );
    if( endTime == 0 ) {
        endSample = firWave->GetLength();
    }
    
    if( !positiveGoing ) {
        
        for( Int_t i = startSample; i < endSample; i++ ) {
            if( okToCount && firWave->At(i) >= threshold ) {
                while( i < firWave->GetLength() - 1 ) {
                    if( firWave->At(i) < 0 ) {
                        //                    find zero crossing
                        //                    for now just approximate
                        pulseTimes.push_back( (firWave->GetTimeAtIndex(i-1) + firWave->GetTimeAtIndex(i))/2 );
                        break;
                    }
                    i++;
                }
                //            pulseTimes.push_back( firWave->GetTimeAtIndex(i) );
                okToCount = false;
                //            could speed up run time by advancing a bit
                //            i += 10;
            }
            if( !okToCount && firWave->At(i) < threshold ) {
                okToCount = true;
            }
        }
    }
    else {
        for( Int_t i = startSample; i < endSample; i++ ) {
            if( okToCount && firWave->At(i) <= threshold ) {
                while( i < firWave->GetLength()-1 ) {
                    if( firWave->At(i) > 0 ) {
                        pulseTimes.push_back( ( firWave->GetTimeAtIndex(i-1) + firWave->GetTimeAtIndex(i))/2 );
                        break;
                    }
                    i++;
                }
                okToCount = false;
            }
            if( !okToCount && firWave->At(i) > threshold ) {
                okToCount = true;
            }
        }
    }
    
    return pulseTimes;
}



// as of now (march 16 2016) this skips 40 ns after finding a pulse
// so it won't effectively discern between closely correlated pulses, but this kind of pileup would be tricky anyway
template<typename _Tp>
std::vector<Double_t> getPulseTimesFromEdgeCounting( TTemplWaveform<_Tp>* wave, baseline_t baseline, bool invert, Double_t threshold, Double_t startTime, Double_t endTime, Double_t skipTime) {
    
    std::vector<Double_t> pulseTimes;
    
    //    if( baseline < 0 ) {
    //        baseline = getBaseline( wave, wave->At(0) );
    //    }
    
    
    /*
     Int_t startIndex;
     if( startTime < 0 ) {
     startIndex = 1;
     }
     else {
     startIndex = wave->GetIndexAtTime( startTime );
     }
     
     Int_t endIndex;
     if( endTime < 0 ) {
     endIndex = wave->GetLength() - 1;
     }
     if( endTime > wave->GetLength() * wave->GetSamplingPeriod() ) {
     endIndex = wave->GetLength() - 1;
     }
     else {
     endIndex = wave->GetIndexAtTime( endTime );
     }
     */
    Int_t startIndex, endIndex;
    Int_t indices[2];
    getIndicesFromTimeInWave( wave, startTime, endTime, indices );
    startIndex = indices[0];
    endIndex = indices[1];
    
    bool okToCount = true;
    Double_t sampleValue  = 0;
    
    for( Int_t index = startIndex; index < endIndex; index++ ) {
        
        if( invert )
        sampleValue = -1*(wave->At(index)-baseline.baseline);
        else
        sampleValue = wave->At(index)-baseline.baseline;
        
        if( !okToCount && (sampleValue < (threshold * 0.2)) )
        okToCount = true;
        if( okToCount && (sampleValue > threshold) ) {
            pulseTimes.push_back( wave->GetTimeAtIndex(index) );
            okToCount = false;
            // jump ahead 40 ns
            index += (Int_t) skipTime / wave->GetSamplingPeriod();
        }
    }
    
    
    return pulseTimes;
}




// returns time (in ns) of next BPM zero crossing, starting from startSearchTime, in ns
Double_t getNextBPMzeroCrossing( TTemplWaveform<Short_t>* bpmWaveform, Double_t startSearchTime, Double_t baseline ) {
    
    baseline = getBaseline( bpmWaveform, 8150 );
    Double_t zeroTime = 0;
    
    Int_t startIndex = bpmWaveform->GetIndexAtTime( startSearchTime );
    if( startIndex < 0 ) {
        startIndex = 0;
    }
    if( startIndex >= bpmWaveform->GetLength() ) {
        startIndex = bpmWaveform->GetLength() - 1;
    }
    bool searchForZeroCross = false;
    for( Int_t i = startIndex; i < bpmWaveform->GetLength(); i++ ) {
        if( searchForZeroCross ) {
            if( bpmWaveform->At(i) -baseline < 0 ) {
                //                do a simple interpolation for now
                zeroTime = bpmWaveform->GetTimeAtIndex(i-1) + 1;
                return zeroTime;
            }
        }
        else {
            if( bpmWaveform->At(i) - baseline > 500 ) {
                searchForZeroCross = true;
            }
        }
    }
    
    return 0;
}


//
// this function produces the filtered and jump-rejecting baseline
// it populates the array passed into it as address rollingPwave
// assumes rollingPwave is already allocated and appropriate in size
// width sets the width of the averaging window
//
template<typename _Tp>
void getRollingPwave( TTemplWaveform<_Tp>* waveform, double* rollingPwave,
                     int halfWidth, double preloadValue,
                     double rejectThreshold ) {
    
    std::deque<double> movingBaselineFilter;
    baseline_t baselineDummy;
    baselineDummy.baseline = 0;
    baselineDummy.FWHM = 0;
    baselineDummy.chisquare = 0;
    

    
    double movingBaselineValue = 0;
    

    if( preloadValue > -7777 ) {
        /* check if preloadValue argument is a 'valid' sample value */
        /* default argument is -7777 */
        /* mildly arbitrary, but decidedly invalid for samples between 0 and 16383 */
        /* if a good preload value is supplied: */
        /*      then we treat the preload value as the "first half" of the of the boxcar at the start of the wave */
        /*      this helps stabilize the filter against the presence of signal in the first bit of the waveform */
        for( int i = 0; i < halfWidth; i++ ) {
            movingBaselineFilter.push_back( preloadValue );
            movingBaselineValue += preloadValue;
        }
    }
    
    if( DEBUG >= 2 ) {
        printf("after preload, deque has %lu entries\n", movingBaselineFilter.size() );
        printf("average value is %.2f\n", movingBaselineValue / movingBaselineFilter.size() );
    }
    
    
    
    for( int i = 0; i < halfWidth; i++ ) {
        if( preloadValue > -7777 ) {
            if( fabs(getWaveformValueSafely( waveform, baselineDummy, i, false ) - movingBaselineValue/movingBaselineFilter.size() ) >= rejectThreshold ) {
                if( DEBUG >= 2 ) {
                    printf("skipping early value %.2f at sample %i\n", getWaveformValueSafely( waveform, baselineDummy, i, false), i );
                }
                continue;
            }
        }
        movingBaselineFilter.push_back( getWaveformValueSafely(waveform, baselineDummy, i, false ) );
        movingBaselineValue += getWaveformValueSafely( waveform, baselineDummy, i, false );
    }
    
    
    for( int i = 0; i < waveform->GetLength(); i++ ) {
        
        if( fabs(getWaveformValueSafely( waveform, baselineDummy, i+halfWidth, false ) - movingBaselineValue/movingBaselineFilter.size() ) < rejectThreshold  ) {
            if( i + halfWidth < waveform->GetLength() ) {
                /* if we're here, we're still within valid lengths of the waveform */
                if( movingBaselineFilter.size() >= halfWidth*2 + 1 ) {
                    /* here, the filter is fully occupied */
                    /* so, pop values off the back as we go */
                    movingBaselineValue -= movingBaselineFilter.front();
                    movingBaselineFilter.pop_front();
                }
                
                movingBaselineValue += getWaveformValueSafely( waveform, baselineDummy, i+halfWidth, false );
                movingBaselineFilter.push_back(getWaveformValueSafely( waveform, baselineDummy, i+halfWidth, false ) );
            }
            else {
                movingBaselineValue -= movingBaselineFilter.front();
                movingBaselineFilter.pop_front();
            }
        }
        rollingPwave[i] = movingBaselineValue / movingBaselineFilter.size();
    }
    
    
    if( DEBUG >= 2 ) {
        printf("moving baseline deque size: %lu\n", movingBaselineFilter.size() );
        printf("getRollingPwave - dumping first 10 values of waveform...\n");
        for( int i = 0; i < 10; i++ ) {
            printf("%.2f\t", rollingPwave[i] );
        }
        printf("\n");
    }
}


/*
 getCMAfilteredWave
 
 perform CMA filtering - this returns a pointer to the waveform resulting from subtraction of the CMA-determined baseline
 inversion can be performed as well
 
 assumes that rollingPwave points to a pre-allocated array with sufficient size to include the entire length of waveform
 */
template<typename _Tp>
TTemplWaveform<Double_t>* getCMAfilteredWave( TTemplWaveform<_Tp>* waveform,
                                             double* rollingPwave,
                                             double* cmaArray,
                                             int halfWidth,
                                             double preloadValue,
                                             double rejectThreshold, bool invert ) {
    getRollingPwave( waveform, rollingPwave, halfWidth, preloadValue, rejectThreshold );
    int coeff = 1;
    if( invert ) {
        coeff = -1;
    }
    for( int idx = 0; idx < waveform->GetLength(); idx++ ) {
        cmaArray[idx] = coeff * (waveform->At(idx) - rollingPwave[idx]);
    }
    
    TTemplWaveform<Double_t>* cmaWave = new TTemplWaveform<Double_t>( cmaArray, waveform->GetLength() );
    cmaWave->SetSamplingFreq(waveform->GetSamplingFreq());
    
    return cmaWave;
}


//
// this performs a triangle filter on the supplied waveform
// populates the array also supplied
//--- assumes the supplied array is properly allocated already
// order specifies how many samples are included
// --- order 3 includes 2 samples on either side of the index
// order must be >= 1
template<typename _Tp>
void getTriangleFilteredWave( TTemplWaveform<_Tp>* waveform, double* triangleFilteredWave, int order ) {
    
//    if order 1 or lower, just populate the array with the raw waveform
    if( order <= 1 ) {
        for( int i = 0; i < waveform->GetLength(); i++ ) {
            triangleFilteredWave[i] = waveform->At(i);
        }
        return;
    }
    
    baseline_t baselineDummy;
    baselineDummy.baseline = 0;
    baselineDummy.chisquare = 0;
    baselineDummy.FWHM = 0;
    
    //    std::vector<double> coeffs;
    //    for( int i = 1; i <= order; i++ ) {
    //        coeffs.push_back
    //    }
    
    for( int i = 0; i < waveform->GetLength(); i++ ) {
        
        double sampleValue = order * getWaveformValueSafely(waveform, baselineDummy, i, false );
        
        for( int idx = order - 1; idx > 0; idx--  ) {
            //            get the sample on the positive side
            sampleValue += idx * getWaveformValueSafely( waveform, baselineDummy, i + idx, false );
            //            get the sample on the negative side
            sampleValue += idx * getWaveformValueSafely( waveform, baselineDummy, i - idx, false );
        }
        //        now we need to scale based on the order
        //        this is a division by the sum of all the coefficients we've included
        //        so scale = order + 2 * (order - 1) + 2 * (order - 2) .. 0
        //        there's probably a simple mathematical expression for this
        //        but just use a for loop
        int coefficientSum = order;
        for( int idx = order - 1; idx > 0; idx-- ) {
            coefficientSum += idx;
            //            only add the coefficient twice if it is inside the waveform
            if( i - idx >= 0 && i + idx < waveform->GetLength() ) {
                coefficientSum += idx;
            }
        }
        
        //        do the scaling
        sampleValue /= coefficientSum;
        
        //        and set the waveform value equal to our filtered value
        triangleFilteredWave[i] = sampleValue;
        
        //        triangleFilteredWave[i] = (getWaveformValueSafely( waveform, baselineDummy, i-2, false ) + 2 * getWaveformValueSafely( waveform, baselineDummy, i-1, false ) + 3 * getWaveformValueSafely( waveform, baselineDummy, i, false ) + 2 * getWaveformValueSafely( waveform, baselineDummy, i+2, false ) + getWaveformValueSafely( waveform, baselineDummy, i+1, false ) )/9;
    }
    
    if( DEBUG >= 2 ) {
        printf("getTriangleFilteredWave - dumping first 10 values, starting at 0, of waveform...\n");
        for( int i = 0; i < 10; i++ ) {
            printf("%.2f\t", triangleFilteredWave[i] );
        }
        printf("\n");
    }
}


template<typename _Tp>
bool sharpEventCheck_integral( TTemplWaveform<_Tp>* wave, baseline_t baseline, bool invert, Double_t startTime, Double_t regionLength, Double_t ratioThreshold ) {
    
    Double_t q1 = integrateWaveform( wave, baseline, startTime, startTime + regionLength, invert );
    Double_t q2 = integrateWaveform( wave, baseline, startTime + regionLength, startTime + 2 * regionLength, invert );
    
//    printf("region 1 %.2f\tregion 2 %.2f\n", q1, q2 );
    
    if( q2 / q1 < ratioThreshold ) {
        return true;
    }
    return false;
}



template<typename _Tp>
bool sharpEventCheck_aq( TTemplWaveform<_Tp>* wave, baseline_t baseline, bool invert, Double_t startTime, Double_t regionLength, Double_t ratioThreshold ) {
    
    Double_t q = integrateWaveform( wave, baseline, startTime, startTime + regionLength, invert );
    Double_t a = getPeakValue( wave, baseline, startTime, startTime+regionLength, invert );
    
//    printf("region Q %.2f\tregion A %.2f\n", q, a );
    
    if( a / q > ratioThreshold ) {
        return true;
    }
    return false;
}


#endif

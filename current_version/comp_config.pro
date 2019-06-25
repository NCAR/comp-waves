FUNCTION comp_config

; I/O
root_dir='analysis/'
numberWaveLengths = 3
filenameExtra=''


; Cross-correlation parameters
ccThreshold =  0.04
maxNumberIterations = 200

;Interploation of missing frame
;Number of coefficients to use for Maximum entropy method
maxentCoefficientNumber=10

; Sun parameters
radiusSunMm = 695.7 ;solar radius in Mm

; Wave angle calculations
coherenceBoxHalfLength  = 10
coherenceSmoothing = 15
coherenceLimit = 0.5
filterWidth  = 0.0015
filterCentralFrequency = 0.0035
minNumberCoherentPixels=10

; Basic mask parameters 
upperMaskRadius = 280 ; pixel radius just greater than 1.3 Rsun
lowerMaskRadius = 220 ; pixel radius just greater than 1.02 Rsun

;Space time run parameters
maxTrackLength=25


return,DICTIONARY('numberWaveLengths', numberWaveLengths, 'radiusSunMm', radiusSunMm, 'lowerMaskRadius', lowerMaskRadius,$ 
                  'upperMaskRadius', upperMaskRadius, 'ccThreshold', ccThreshold,'maxNumberIterations',maxNumberIterations, $
                   'coherenceBoxHalfLength',coherenceBoxHalfLength,'coherenceSmoothing',coherenceSmoothing,'coherenceLimit',coherenceLimit,$
                   'minNumberCoherentPixels',minNumberCoherentPixels, $
                   'filterWidth', filterWidth,'filterCentralFrequency',filterCentralFrequency,'maxentCoefficientNumber',maxentCoefficientNumber, $
                   'root_dir',root_dir,'filenameExtra',filenameExtra,'maxTrackLength',maxTrackLength)


END

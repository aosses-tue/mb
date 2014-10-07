% AMT - Single stages of auditory models
%
%  The AMT team, 2011 - 2014.
%
%  Peripheral stages
%     HEADPHONEFILTER    - FIR filter to model headphones
%     MIDDLEEARFILTER    - FIR filter to model the middle ear
%     AUDITORYFILTERBANK - Linear auditory filterbank.
%     IHCENVELOPE        - Inner hair cell envelope extration.
%     ADAPTLOOP          - Adaptation loops.
%     KARJALAINEN1996    - Adaptation network.
%     DRNL               - Dual resonance non-linear filterbank.
%     MODFILTERBANK      - Modulation filter bank.
%     MODFILTERBANKEPSM  - Modulation filter bank from Joergensen model.
%
%  Binaural processing stages
%     LANGENDIJK2002COMP     - Comparison process from Langendijk 20002.
%     LINDEMANN1986BINCORR - Running cross-correlation between two signals.
%     EICELL             - Excitation-inhibition cell model by Breebaart.
%     CULLING2005BMLD    - BMLD calculation from Culling et al. (2005).
%     DIETZ2011UNWRAPITD - IPD to ITD transformation for the Dietz model
%
%  The Takanen 2013 model   - Takanen 2013 binaural auditory model
%     TAKANEN2013CONTRACOMPARISON - Enhance contrast between hemispheres
%     TAKANEN2013CUECONSISTENCY - Check consistency before cue combination
%     TAKANEN2013DIRECTIONMAPPING - Map the directional cues to directions
%     TAKANEN2013FORMBINAURALACTIVITYMAP - Steer cues on a topographic map
%     TAKANEN2013LSO        - Model of the lateral superior olive
%     TAKANEN2013MSO        - Model of the medial superior olive
%     TAKANEN2013ONSETENHANCEMENT - Emphasize onsets on direction analysis
%     TAKANEN2013PERIPHERY  - Process input through the model of periphery
%     TAKANEN2013WBMSO      - Wideband medial superior olive model
%
%  Utility
%     PMV2PPP            - PMV to PPP conversion
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/Contents.php

% Copyright (C) 2009-2014 Peter L. Søndergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.



function f = processWindow(obj, dataIn)
% PROCESSWINDOW This is the core workhorse of this analyser
%
% NOTE: Every window is the same size and comes appropriately
%       zero-padded both in front and at rear, whenever neccessary,
%       and so that is why the following is just a straightforward
%       call.

% Create a function handle for roughness
fs = get(obj, 'fs');
fH = roughness(fs);

% function handle to return
f = @run;

  % Local nested function:
  function dataOut = run(dataIn)
  % function dataOut = run(dataIn)
  %
  % 1. Description: 
  %         Local nested function.
  %
  %         The following multiplication was added to set levels correctly for roughness.
  %         A 100% amplitude modulated tone, fmod = 70 Hz, 60 dB SPL comes out at 0.91 asper,
  %         However, it should be 1 asper. The discrepancy is a gain issue, and identical 
  %         results can be obtained by using an offset of approx 4.867dB. We (PsySound team)
  %         don't want to change SLM, so offset is applied here.
        dataOut = fH(dataIn * 1.71); 
  end

end % processWindow
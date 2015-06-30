function out = Autoc_opt2_Alt_proc(map, input)
% function out = Autoc_opt2_Alt_proc(map, input)
%
% Autoc_proc_opt2: Main autocorrelation-based F0 extractor. (2nd
% optimization)
%
% Inputs:
% map:      Autocorrelation map
%           use map.CompareSimulink = 1 for direct comparison to Simulink
%           (using a normal delay instead of a non-delay buffer)
%
% inputs:   Audio file to process (Processing only)
% 
% Outputs:
%   out     = {input, F0, map.is_silence, tF0}, passes the input and gives 
%           extracted F0 (registered at time tF0)
%
% % Stand-alone example:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: K.U.Leuven, Lab. Exp. ORL, Matthias Milczynski
%        $Date: 04/09/2006 $
%      Authors: Matthias Milczynski, Koen Eneman
% Edited by Alejandro Osses, ExpORL, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contador = 0;

switch nargin  
    case 0
        return;
    case 1
    %%%%%%%%%%%%%%%%
    % Initialization
    %%%%%%%%%%%%%%%%
    % Setup parameters
        map = Ensure_field(map,'audio_sample_rate',15675.375);
    
        map.sr      = map.audio_sample_rate;
        
        map = Ensure_field(map, 'CompareSimulink', 0);
        map = Ensure_field(map, 'min_f0', 75);  % Hz
        map = Ensure_field(map, 'max_f0', 300); % Hz
        map = Ensure_field(map, 'min_lag', round(map.sr/map.max_f0));
        map = Ensure_field(map, 'max_lag', round(map.sr/map.min_f0));
        map = Ensure_field(map, 'vUv_level', 0);
        map = Ensure_field(map, 'vUv_thr', 0.55);
        
        % Set to 1 when interpolation is desired
        map = Ensure_field(map, 'interpol_level', 0);
        
        % Most feasible robustness level 
        map = Ensure_field(map, 'robust_level', 1);
        
        % Set to 1 when to look for dominant subharmonic
        map = Ensure_field(map, 'subharm_level', 1);
        
        % Set to 1 when quadratic interpolation in subharmonic detection is
        % desired
        map = Ensure_field(map, 'subharm_interpol', 0);
        map = Ensure_field(map, 'subharm_thr', 0.45); 
        
        % In the first step the incoming signal is lowpass filterd to remove 
        % "disturbing" Formants -> oldfashioned ?        
        map = Ensure_field(map, 'fir_ord', 32); 
        map = Ensure_field(map, 'fir_cf', 900);
        
        % Silence threshold and maximal frequency (discrete) gradient
        map = Ensure_field(map, 'silence_thr', 0.05);  
        map = Ensure_field(map, 'max_freq_grad', 30);  
                
        % Used for clipping
        map = Ensure_field(map, 'peak_fac', 0.64);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate filter coefficients -----------------------------------
        fir_ord = map.fir_ord;
        fir_cf  = map.fir_cf;
   
        map.fir_coeffs  = fir1(fir_ord, fir_cf*2/map.sr);
        map.ac_gd = (length(map.fir_coeffs) - 1)/2;
        map = Ensure_field(map, 'autoc_block_size', 432);
        map = Ensure_field(map, 'autoc_block_shift', 144);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        map = Ensure_field(map, 'delta_lag_rob', round(0.001*map.sr));
        map = Ensure_field(map, 'delta_lag_subh', round(0.000125*map.sr));
        % do_lpf = 1 equals true
        map = Ensure_field(map, 'do_lpf', 1);
        map = Ensure_field(map, 'do_plot_F0', 0);
        map = Ensure_field(map, 'plot_pause', 0);
        map = Ensure_field(map, 'gain_dB_ACF_input', 0); % Added by AOV
        out = map;
        
        %disp('Init done')
    case 2
    %%%%%%%%%%%%
    % Processing
    %%%%%%%%%%%%    
        if isfield(map,'DEBUG')
            map.DEBUG=0;
        end
        input = From_dB(map.gain_dB_ACF_input)*input; % Added by AOV, for calibration dBFS
        
        % LPF at 900 Hz cutoff
        t = (0:length(input)-1)/map.sr;
        if(map.do_lpf)
            audio_in = filter( map.fir_coeffs, 1, ... FIR Coeffs
                              [input(:);zeros(map.ac_gd, 1)]); % Signal + 0 padding            
            %Remove transient
            audio_in = audio_in((map.ac_gd + 1):end);
        else
            audio_in = input;
        end            
        
        audio_in = audio_in - mean(audio_in); % Makes the new mean(audio_in) = 0...
        block_shift = map.autoc_block_shift;
        block_size  = map.autoc_block_size;
        max_amp     = max(abs(audio_in));
        audio_len   = length(audio_in);
        res_len     = round(audio_len/block_shift);
        
        tF0         = (1:audio_len) / map.audio_sample_rate;
        tF0         = buffer(tF0, block_size, block_size - block_shift, 'nodelay');
        tF0         = tF0(1,:);
        
        F0prev      = 0;
        
        % Use buffer function to simplify further computation (especially AM)
        if map.CompareSimulink == 0
            sig_blocks = buffer(audio_in, ...
                                block_size, ...
                                block_size - block_shift, ...
                                'nodelay'   ); % No delay supresses all the zero padding before the first column
        else
            sig_blocks = buffer(audio_in, ...
                                block_size, ...
                                block_size - block_shift);
        end
        
        num_blocks = size(sig_blocks, 2);
        ci_clip    = zeros(size(sig_blocks,1), size(sig_blocks,2));
        auto_corr  = zeros(num_blocks, map.max_lag);
        ac_idx_pairs = zeros(2, num_blocks);
        F0           = nan(1, num_blocks);
        if map.plot_pause
           pdat = [];
           pdat.ac_blockShift = map.autoc_block_shift;
           pdat.ac_blockLength = map.autoc_block_size;
           pdat.vUvThr = map.vUv_thr;
           pdat.subharmThr = map.subharm_thr;
        end
        
        map.is_silence = zeros(num_blocks,1);
        
        for i=1:num_blocks
            section = sig_blocks(:, i);
            if ( max(abs(section)) > map.silence_thr )
                map.is_silence(i) = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Detect peaks at the beginning (first 10ms = 144/16000) and 
                % at the end (last 10ms) of the current section
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                peak_1 = max(abs(section(1:block_shift))); 
                peak_2 = max(abs(section((block_size - block_shift + 1):block_size)));
        
                clip_lev = map.peak_fac*min(peak_1,peak_2);  % See formula Page 29, Matthias' thesis
                % In Simulink model clip_lev will be renamed to beta...
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Center and inifinte peak clipping
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ci_clip(:, i) = sign(section).*(abs(section) > clip_lev);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Core Autocorrelation algorithm -- begin
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                auto_norm = max(eps, ci_clip(:, i)'*ci_clip(:, i)/block_size);
                for l = map.min_lag:map.max_lag
                    auto_corr(i, l) = ci_clip(1:block_size-l, i)'*ci_clip(l+1:block_size, i)/(block_size-l);
                end
        
                auto_corr(i, :) = auto_corr(i, :)/auto_norm;
                [auto_val, auto_idx] = max(auto_corr(i, :));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Core Autocorrelation algorithm -- end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Detect dominant subharmonics
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                map.subharm_level;
                if map.subharm_level
                    [auto_val, auto_idx] = subharm_detect(map, auto_corr(i, :), auto_val, auto_idx);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % Interpolate when desired
                %%%%%%%%%%%%%%%%%%%%%%%%%%
        
                if map.interpol_level
                    [auto_val, auto_idx] = quad_interpol(auto_corr(i, :), auto_val, auto_idx, 1, map.max_lag, res_len);
                end
                                
                ac_idx_pairs(1,i) = auto_val;
                ac_idx_pairs(2,i) = auto_idx;
  
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Determine v/uv detection approach
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                switch map.vUv_level
            
                    case 0
                        auto_thr = auto_val;    
                    case 1
                        auto_thr = clip_lev;
                    case 2
                        auto_thr = clip_lev*auto_val;
                    otherwise
                        error('vUv_level can be either 0, 1 or 2!');
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate frequency from autocorrelation index
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (auto_thr > map.vUv_thr)
                    F0(i) = map.sr / auto_idx;  % Omega / Tau_best (pag 30, Matthias' Thesis)
                                                % Fs    /
                    tmp_vUv_thr = map.vUv_thr;
                end
                
                if auto_thr > map.vUv_thr
                    contador = contador + 1;
                else
                    contador = contador + 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply optimization/postprocessing when desired
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                if ((i > 2) && map.robust_level >= 1)
                    [F0_suc, F0_prev] = post_proc_F0( F0(i    ), ...
                                                      F0(i - 1), ... 
                                                      F0(i - 2),...
                                                      map.min_f0,...
                                                      map.max_f0,...
                                                      map.max_freq_grad);
                                          
                    F0(i) = F0_suc;
                    if map.robust_level == 2
                        F0(i - 1) = F0_prev;
                    end
                end
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Handle "Fade-out" problem when desired
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                if map.robust_level >= 1
                    if ((i > 2) && (F0(i) == 0) && (F0(i - 1) > 0))
                        idx_ref = round(1/F0(i - 1)*map.sr);
                        if (abs(auto_idx - idx_ref) < map.delta_lag_rob) && ...
                           (auto_thr > tmp_vUv_thr/2)
                           % !! Here we can run out boundaries minF0 and max F0 ... !!
                           
                           F0(i) = 1/(auto_idx - 1)*map.sr;
                           
                           if (auto_thr < tmp_vUv_thr)
                               tmp_vUv_thr = tmp_vUv_thr/2;
                           end
                        else
                            tmp_vUv_thr = map.vUv_thr;
                        end
                    end  
                end
                % ... thus check again for correct range (not a nice implementation !!)
                if((F0(i) < map.min_f0) || (F0(i) > map.max_f0))
                    F0(i) = 0;
                    ac_idx_pairs(1, i) = 0;                    
                end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % Relates to upper most if-condition
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if map.plot_pause
                   F0EstPlotData(t, input, audio_in, i, ci_clip(:, i), ...
                       auto_corr(i, :), auto_thr, auto_idx, pdat, 16);
                    pause;
                    clf;
                end
            else
                F0(i)=NaN;    
            end
        end
        
        if( map.do_plot_F0 )
           
            figure;
            plot(F0);
            ylim([0 340]);
            grid on;
            
        end
        input = From_dB(-1*map.gain_dB_ACF_input)*input; % Back to original scale
        out = {input, F0, map.is_silence, tF0};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End main autoc_proc function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val, idx] = quad_interpol(auto_corr, max_val, max_idx, min_lag, max_lag, max_len)

if (max_val > min_lag) && (max_val < max_lag) && (max_idx < max_len)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perfrom quadratic interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = max_val;
    
    b =(auto_corr(max_idx + 1)- auto_corr(max_idx - 1))/2;
    a = auto_corr(max_idx + 1) - max_val - b;
    idx = -b/2/a + max_idx;
    val = -b*b/4/a + c;
    
else
    idx = max_idx;
    val = max_val;
end


%%%%%%%%%%%%%%%%%%%
% End quad_interpol
%%%%%%%%%%%%%%%%%%%


function [F0_cur, F0_prev] = post_proc_F0(F0_cur, F0_prev, F0_prevprev, min_F0, max_F0, freq_grad)

if (F0_cur < min_F0)||(F0_cur > max_F0)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Too low level -> F0 = -1
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if F0_cur~=-1
        F0_cur=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete gradient of F0
% Allow F0 to de/increase at most by freq_grad per block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (abs(F0_prev - F0_cur) > freq_grad) && (F0_cur > 0) && (F0_prev > 0)
    if F0_cur>F0_prev
        F0_cur = F0_prev + freq_grad;   
    else
        F0_cur = F0_prev - freq_grad; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set currently extracted F0 to 0 when outside given range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (F0_cur < min_F0) || (F0_cur > max_F0)
    if F0_cur~=-1
        F0_cur=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!SUCCESSIVE STEPS ARE CURRENTLY NOT APPLIED!!
% Interpolate outlier F0prev 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (F0_cur > 0) && (F0_prevprev > 0) && (F0_prev <= 0)
    F0_prev = (F0_prevprev + F0_cur)/2;
end

%%%%%%%%%%%%%%%%
% Same as before
%%%%%%%%%%%%%%%%

if (F0_cur <= 0) && (F0_prevprev <= 0) && (F0_prev > 0)
    F0_prev=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct/Interpolate too high gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (abs(F0_cur - F0_prevprev) < 5) && (abs(F0_cur - F0_prev) > 20)
    F0_prev=(F0_prevprev + F0_cur)/2;
end

%%%%%%%%%%%%%%%%%%
% End post_proc_F0
%%%%%%%%%%%%%%%%%%

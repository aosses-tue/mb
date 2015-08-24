function [out_1avg out_2avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,in_signal_pre,fs,model,sigma,Ntimes,idx_fc)
% function [out_1avg out_2avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,in_signal_pre,fs,model,sigma,Ntimes,idx_fc)
%
% 1. Description:
%       Obtains an averaged internal representation of the signal in_signal_pre
%       superimposed on a random selected sample of the buffered noise in_masker_pre.
%       
%       in_masker_pre - buffered 'noise'
%       in_signal_pre - current test signal
% 
%           ir1 - related to 'Noise alone'
%           ir2 - related to 'Suprathreshold signal' (signal well above threshold)
%           setup.fs - Sampling frequency: relevant for the normalisation process.
% 
% 2. Stand-alone example:
%           filename = AMTControl_Examples(0);
%           [insig1 fs]      = Wavread(filename{1}); 
%           [insig2supra fs] = Wavread(filename{2}); insig2supra = insig2supra(1:0.1*fs);
%           Get_internalrep_stochastic(insig1, insig2supra,fs,'dau1996',0.85,10);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/08/2015
% Last update on: 13/08/2015 
% Last use on   : 13/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    Ntimes = 10;
end

if nargin < 5
    sigma = 0.8;
end

N = size(in_signal_pre,1);

out_1 = [];
out_2 = [];
for i = 1:Ntimes
    
    in_masker = Randomise_insig(in_masker_pre); % random sample of the noise
    in_masker = in_masker(1:N,:);
    
    switch model
        case 'dau1996'
            [out_pre , fc] = dau1996preproc(in_masker,fs); % out_1pre affected by external noise
            fs_intrep = fs;
            mu = 0;
            if nargin < 7
                idx_fc = 1:length(fc);
            end
            [Ni,Mi] = size(out_pre(:,idx_fc));
            tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
            out_1 = [out_1 tmp(:)]; % out_1 affected by internal noise
            
            [out_pre , fc] = dau1996preproc(in_masker + in_signal_pre,fs); % out_2pre affected by external noise
            tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
            out_2 = [out_2 tmp(:)]; % out_1 affected by internal noise
            
        case 'dau1997'
            
            [out_1pre , fc, mfc] = dau1997preproc_1Ch(in_masker                ,fs,fc);
            [out_2pre , fc, mfc] = dau1997preproc_1Ch(in_masker + in_signal_pre,fs,fc);
            fs_intrep = fs;
            
            if bDeterministic == 1 
                % Deterministic noise
                [out_1 noise] = Add_gaussian_noise_deterministic(out_1pre,mu,sigma); 
                out_2 = out_2pre+noise; % deterministic noise
                
            else
                % 'Running' noise
                [n m] = size(out_1pre);
                out_1 = [out_1 Add_gaussian_noise(out_1pre(:),mu,sigma)]; 
                out_2 = [out_2 Add_gaussian_noise(out_2pre(:),mu,sigma)]; % Add internal noise
            end
            
            % tmp.fs = fs;
            % if bDeterministic
            % 
            %     out_1Mean = out_1;
            %     out_2Mean = out_2;
            % 
            % else
            % 
            %     out_1Mean = mean( out_1,2 );
            %     out_2Mean = mean( out_2,2 );
            % 
            %     out_1Mean = reshape(out_1Mean,n,m);
            %     out_2Mean = reshape(out_2Mean,n,m);
            % 
            %     bPlot = 0;
            %     if bPlot
            %         for i = 1:length(mfc)
            %             figure;
            %             plot(out_1Mean(:,i)), hold on;
            %             plot(out_2Mean(:,i),'r'); grid on
            %             plot(out_2Mean(:,i)-out_1Mean(:,i),'g')
            %             title(sprintf('mfc=%.1f Hz',mfc(i)));
            %         end
            %     end
            % end
            
        otherwise
            error('Model not added yet')
    end
    
end
   
out_1avg = mean(out_1,2);
out_2avg = mean(out_2,2);

out_1avg = reshape(out_1avg,Ni,Mi);
out_2avg = reshape(out_2avg,Ni,Mi);

outs.fc = fc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

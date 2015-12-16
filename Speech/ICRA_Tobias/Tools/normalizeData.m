function [out,STATS] = normalizeData(in,method,dim)
%normalizeData   Perform normalization of input data.
%
%USAGE
%   [OUT,STATS] = normalizeData(IN)
%   [OUT,STATS] = normalizeData(IN,METHOD,DIM)
%   
%INPUT ARGUMENTS
%       IN : data matrix arranged as [nSamples x nChannels]
%   METHOD : string specifying normalization method
%               'mean' - normalize data to have zero mean
%                'var' - normalize data to have unit variance
%            'meanvar' - normalize data to have zero mean and unit variance
%                'max' - normalize data to its maximum
%              'range' - normalize data such that it ranges between [0 1]
%                'heq' - histogram equalization 
%                        (default, METHOD = 'meanvar')
%      DIM : dimension along which the input data should be normalized
% 
%OUTPUT ARGUMENTS
%      OUT : normalized data [nSamples x nChannels]
%    STATS : structure array with the following fields
%            .method - normalization method
%            .dim    - dimension 
%            .data   - normalization statistics 
% 
%   See also scaleData.
% 
%ACKNOWLEDGEMENT
%   The histogram equalization has been adopted from Marc René Schädler's
%   reference feature implementation available at gitup: 
%   https://github.com/m-r-s/reference-feature-extraction
%
%NOTE
%   When training a classifier with features, this function is typically
%   used in combination with the function scaleData:  First, the feature
%   matrix is normalized during training. Then, the normalization
%   statistics measured during training are applied during testing to scale
%   the features accordingly. Sharing the same normalization factors
%   ensures that a particular feature value has the same relevance in the
%   training and the testing stage.   
% 
%   % Create 2D feature matrix
%   featureSpace = randn(10,250);
%
%   % Split data into training and testing set
%   fTrain = featureSpace(:,1:100);
%   fTest  = featureSpace(:,101:end);
% 
%   % Normalize feature matrix during training 
%   [fTrainNorm,STATS] = normalizeData(fTrain,'meanvar',2);
%
%   % Normalize feature matrix during testing using the normalization
%   % statistics measured during training
%   fTestNorm = scaleData(fTest,STATS);


%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/02/25
%   v.0.2   2015/05/27 added histogram equalization
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 2 || isempty(method); method = 'meanvar';       end
if nargin < 3 || isempty(dim);    dim    = findNonSDim(in); end

% Check dimensionality of input data
if numel(size(in)) > 2
    error('Input must be two-dimensional')
end

% Flag for standard deviation computation (scaled by N - 1)
stdFlag = 0;


%% PERFORM NORMALIZATION
% 
% 
% Create normalization structure
STATS = struct('method',method,'dim',dim,'data',[]);

% Select normalization method
switch lower(method)
    case 'mean'
        
        % Compute mean
        STATS.data = mean(in,dim);
        
        % Zero mean normalization
        out = bsxfun(@minus,in,STATS.data);
        
    case 'var'
       
        % Compute standard deviation
        STATS.data = std(in,stdFlag,dim);
        
        % Prevent division by zero
        STATS.data(STATS.data == 0) = 1;
        
        % Check normalization factors
        if any(~isfinite(STATS.data))
            error('Normalization constant is not finite')
        end
        
        % Unit variance normalization
        out = bsxfun(@rdivide,in,STATS.data);
        
    case 'meanvar'
        
        % Compute mean and standard deviation
        STATS.data = mean(in,dim);
        STATS.data = cat(dim,STATS.data,std(in,stdFlag,dim));
        
        % Prevent division by zero
        if dim == 1
            STATS.data(2,STATS.data(2,:) == 0) = 1;
        else
            STATS.data(STATS.data(:,2) == 0,2) = 1;
        end
        
        % Zero mean and unit variance normalization
        if dim == 1
            out = bsxfun(@minus,in,STATS.data(1,:));
            out = bsxfun(@rdivide,out,STATS.data(2,:));
        else
            out = bsxfun(@minus,in,STATS.data(:,1));
            out = bsxfun(@rdivide,out,STATS.data(:,2));
        end
        
    case 'max'
        
        % Compute maximum value
        STATS.data = max(in,[],dim);
        
        % Prevent division by zero
        STATS.data(STATS.data == 0) = 1;
            
        % Check normalization factors
        if any(~isfinite(STATS.data))
            error('Normalization constant is not finite')
        end
        
        % Perform maximum normalization
        out = bsxfun(@rdivide,in,STATS.data);
        
    case 'range'
        
        % Compute minimum and maximum value
        STATS.data = min(in,[],dim);
        STATS.data = cat(dim,STATS.data,max(in,[],dim));
        
        % Perform range normalization
        if dim == 1
            out = bsxfun(@minus,in,STATS.data(1,:));
            out = bsxfun(@rdivide,out,STATS.data(2,:)-STATS.data(1,:));
        else
            out = bsxfun(@minus,in,STATS.data(:,1));
            out = bsxfun(@rdivide,out,STATS.data(:,2)-STATS.data(:,1));
        end
             
    case 'heq'
        
        % Number of uniform intervals
        nQuantiles = 100;
                         
        % Calculate the expected minimum and maximum quantiles
        % when drawing 'context' samples from the unknown distribution
        % C.f. Equation 4 in [1]
        context = size(in,dim);
        
        exptected_min = 1-context./(context+1);
        exptected_max = context./(context+1);
        
        % Define source and target quantiles
        STATS.data.qS = linspace(0, 1, nQuantiles);
        STATS.data.qT = linspace(exptected_min, exptected_max, nQuantiles);
        
        % Get source quantiles from data
        STATS.data.q = quantile(in, STATS.data.qS, dim);
        
        % Allocate memory
        out = zeros(size(in));
        
        if dim == 1
            % Loop over all rows
            for ii = 1:size(in,2)
                % Handle the case if all input values are almost equal
                if (STATS.data.q(end,ii) - STATS.data.q(1,ii)) < 100 * eps
                    out(:,ii) = 0.5;
                else
                    mask = [true; diff(STATS.data.q(:,ii)) > 0];
                    out(:,ii) = interp1q(STATS.data.q(mask,ii), ...
                        STATS.data.qT(mask)', in(:,ii));
                end
            end
        else
            % Loop over all rows
            for ii = 1:size(in,1)
                % Handle the case if all input values are almost equal
                if (STATS.data.q(ii,end) - STATS.data.q(ii,1)) < 100 * eps
                    out(ii,:) = 0.5;
                else
                    mask = [true diff(STATS.data.q(ii,:)) > 0];
                    out(ii,:) = interp1q(STATS.data.q(ii,mask)', ...
                        STATS.data.qT(mask)', in(ii,:)');
                end
            end
        end
        
        % Map the quantiles to the Gaussian distribution
        % using the inverse error function
        out = erfinv(out * 2 - 1);
        
    otherwise
        error('Normalization ''%s'' is not supported',method);
end 


function dim = findNonSDim(input)
% Helper function to find the first non-singleton dimension
dim = find(size(input)~=1,1);
if isempty(dim), dim = 1; end


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************
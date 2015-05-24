% 
% Offline simulator for two-dimensional Dynamic Neural Fields,
% 2-layer architecture with memory trace
%
% The model consists of three fields:
% field_u - excitatory layer
% memTrace_u - memory trace of field u
% field_v - inhibitory layer
%
% Parameters for each layer have the corresponding letter as index, with mu
% for the parameters of the memory trace
% 
% Connections in the model:
% u->u (exc.), u->v (exc.), v->u (inh.),
% u->mu (exc.), mu->u (exc.)
% 
% Connection parameters and kernels between two layers are designated with
% the indices for both layers, with the layer that receives the input 
% written first.
%
% The connection from a field to the associated memory trace does not use 
% an interaction kernel, the output of the field is added to the memory
% trace directly. Parameters with index mu refer to the connection from the
% memory trace back to the field.
%
% Examples:
% kernel_uv: interaction kernel for connection v->u
% kernel_mu: interaction kernel from u's memory trace to u


%%%%%%%%%%%%%%%%%%%%
% model parameters %
%%%%%%%%%%%%%%%%%%%%

fieldSize = 101; % must be odd

tau_u = 20; tau_v = 5; % time constants of dynamic fields
tau_mu_build = 1000; tau_mu_decay = 5000; % time constant of memory trace
beta_u = 5; beta_v = 5; % steepness parameters of sigmoid function
h_u = -5; h_v = -5; % resting levels

% c: interaction strength
% sigma: interaction width
% g: strength of global inhibition
c_uu = 25; sigma_uu = 5;
c_uv = 25; sigma_uv = 10; g_uv = 0;
c_vu = 20; sigma_vu = 5;
c_mu = 0; sigma_mu = 5;

q_u = 0; q_v = 0; % noise levels
sigma_q = 0; % width of the noise kernel

kernelWidthMultiplier = 3; % computation faster for lower values, but less precise
% kernelWidthMultiplier should not be set to values below 3

%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time course %
%%%%%%%%%%%%%%%%%%%%%%%%%%

nTrials = 1;
tMax = 250;

% tStoreFields = [100, 200]; % at these time steps the field activities will be stored
tStoreFields = 1:tMax; % use this to store field activities at all time steps
% tStoreFields = 1:5:tMax; % use this to store field activities at every 5th time step

% for example
tStimulusStart = 20;
tStimulusEnd = 100;
stim1 = 8*gauss2d(1:fieldSize, 1:fieldSize, 25, 75, 5, 5); % a localized hill of input
stim2 = 6*gauss2d(1:fieldSize, 1:fieldSize, 75, 25, 5, 5); % a localized hill of input
stim3 = repmat(2*gauss(1:fieldSize, 25, 5), fieldSize, 1); % an input ridge along the 1st dimension
stim4 = repmat(2*gauss(1:fieldSize, 25, 5)', 1, fieldSize); % an input ridge along the 2nd dimension


%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

halfField = floor(fieldSize/2);

% create matrices for field activities
field_u = zeros(fieldSize, fieldSize);
field_v = zeros(fieldSize, fieldSize);
memTrace_u = zeros(fieldSize, fieldSize);

% create 3d matrices to store field activities
history_u = zeros(fieldSize, fieldSize, nTrials * length(tStoreFields));
history_v = zeros(fieldSize, fieldSize, nTrials * length(tStoreFields));
history_mu = zeros(fieldSize, fieldSize, nTrials * length(tStoreFields));

% index of the current position in the history matrices
iHistory = 1;

% set up the interaction kernels
kernelRange_uu = min(round(kernelWidthMultiplier * sigma_uu), floor((fieldSize-1)/2));
kernel_uu = sqrt(c_uu) * gaussNorm(-kernelRange_uu:kernelRange_uu, 0, sigma_uu);
extIndex_uu = [fieldSize - kernelRange_uu + 1 : fieldSize, 1 : fieldSize, 1 : kernelRange_uu];

kernelRange_uv = min(round(kernelWidthMultiplier * sigma_uv), floor((fieldSize-1)/2));
kernel_uv = sqrt(c_uv) * gaussNorm(-kernelRange_uv:kernelRange_uv, 0, sigma_uv);
extIndex_uv = [fieldSize - kernelRange_uv + 1 : fieldSize, 1 : fieldSize, 1 : kernelRange_uv];

kernelRange_vu = min(round(kernelWidthMultiplier * sigma_vu), floor((fieldSize-1)/2));
kernel_vu = sqrt(c_vu) * gaussNorm(-kernelRange_vu:kernelRange_vu, 0, sigma_vu);
extIndex_vu = [fieldSize - kernelRange_vu + 1 : fieldSize, 1 : fieldSize, 1 : kernelRange_vu];

kernelRange_mu = min(round(kernelWidthMultiplier * sigma_mu), floor((fieldSize-1)/2));
kernel_mu = sqrt(c_mu) * gaussNorm(-kernelRange_mu:kernelRange_mu, 0, sigma_mu);
extIndex_mu = [fieldSize - kernelRange_mu + 1 : fieldSize, 1 : fieldSize, 1 : kernelRange_mu];

kernelRange_q = ceil(kernelWidthMultiplier * sigma_q);
kernel_q = gaussNorm(-kernelRange_q : kernelRange_q, 0, sigma_q);
extIndex_q = [fieldSize-kernelRange_q+1:fieldSize, 1:fieldSize, 1:kernelRange_q];


%%%%%%%%%%%%%%
% simulation %
%%%%%%%%%%%%%%

tic

% loop over trials
for i = 1 : nTrials
  
  % set stimuli for this trial
  
  % prepare matrix that holds the stimulus for each time step
  stimulus = zeros(fieldSize, fieldSize, tMax);
  
  % if needed, create a new stimulus pattern for every trial
  stimPosX = round(fieldSize * rand);
  stimPosY = round(fieldSize * rand);
  stimRand = 5*gauss2d(1:fieldSize, 1:fieldSize, stimPosX, stimPosY, 5, 5);
  
  % write the stimulus pattern into the stimulus matrix for all time steps
  % where it should be active
  for j = tStimulusStart : tStimulusEnd
    stimulus(:, :, j) = stim1;
  end
  
  
  % reset field activities
  field_u(1:fieldSize, 1:fieldSize) = h_u;
  field_v(1:fieldSize, 1:fieldSize) = h_v;
  
  % loop over time steps
  for t = 1 : tMax
    % calculating field outputs
    output_u = sigmoid(field_u, beta_u, 0);
    output_v = sigmoid(field_v, beta_v, 0);
    
    % get endogenous input to fields by convolving outputs with interaction
    % kernels
    % padding is done through the expression output_x(extIndex_xx)
    conv_uu = conv2(kernel_uu', 1, output_u(extIndex_uu, :), 'valid');
    conv_uu = conv2(1, kernel_uu, conv_uu(:, extIndex_uu), 'valid');
    
    conv_vu = conv2(kernel_vu', 1, output_u(extIndex_vu, :), 'valid');
    conv_vu = conv2(1, kernel_vu, conv_vu(:, extIndex_vu), 'valid');
    
    conv_uv = conv2(kernel_uv', 1, output_v(extIndex_uv, :), 'valid');
    conv_uv = conv2(1, kernel_uv, conv_uv(:, extIndex_uv), 'valid');
    
    conv_mu = conv2(kernel_mu', 1, memTrace_u(extIndex_mu, :), 'valid');
    conv_mu = conv2(1, kernel_mu, conv_mu(:, extIndex_mu), 'valid');
    
    totalOutput_v = sum(sum(output_v)); % for global inhibition
    
    % create field noise for this timestep
    % prepapre field noise (with spatial correlation if required)
    noise_u = q_u * randn(fieldSize, fieldSize);
    noise_v = q_v * randn(fieldSize, fieldSize);
    if sigma_q > 0
      noise_u = conv2(kernel_q', 1, noise_u(extIndex_q, :), 'valid');
      noise_u = conv2(1, kernel_q, noise_u(:, extIndex_q), 'valid');
      noise_v = conv2(kernel_q', 1, noise_v(extIndex_q, :), 'valid');
      noise_v = conv2(1, kernel_q, noise_v(:, extIndex_q), 'valid');
    end
    
    % update field activities
    field_u = field_u + 1/tau_u * (-field_u + h_u + stimulus(:, :, t) ...
      + conv_uu - conv_uv - g_uv * totalOutput_v + conv_mu) + noise_u;
    field_v = field_v + 1/tau_v * (-field_v + h_v + conv_vu) + noise_v;
    
    % update memory trace (only if there is activity in the field)
    activeRegions_u = field_u > 0;
    if any(any(activeRegions_u))
      memTrace_u = memTrace_u + 1/tau_mu_build * (-memTrace_u + output_u) .* activeRegions_u ...
        + 1/tau_mu_decay * (-memTrace_u) .* (1-activeRegions_u);
    end
    
    % store field activities at the selected time steps
    if any(tStoreFields == t)
      history_u(:, :, iHistory) = field_u;
      history_v(:, :, iHistory) = field_v;
      history_mu(:, :, iHistory) = memTrace_u;
      iHistory = iHistory + 1;
    end
  end
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization of results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nStoredFields = nTrials * length(tStoreFields);

% view mesh plots (3d) of all stored field activities iteratively
if 0
  figure;
  for i = 1 : nStoredFields
    subplot(3, 1, 1);
    mesh(1:fieldSize, 1:fieldSize, history_u(:, :, i));
    set(gca, 'ZLim', [-10 10], 'CLim', [-10 10]);
    zlabel('activity u');
    subplot(3, 1, 2);
    mesh(1:fieldSize, 1:fieldSize, history_mu(:, :, i));
    set(gca, 'ZLim', [0 10], 'CLim', [0 10]);
    zlabel('memory trace u');
    subplot(3, 1, 3);
    mesh(1:fieldSize, 1:fieldSize, history_v(:, :, i));
    set(gca, 'ZLim', [-10 10], 'CLim', [-10 10]);
    zlabel('activity v');
    pause;
  end
end


% view image plots (2d) of all stored field activities iteratively
if 1
  figure('Position', [50 50 900 300]);
  for i = 1 : nStoredFields
    subplot(1, 3, 1);
    imagesc(1:fieldSize, 1:fieldSize, history_u(:, :, i), [-10 10]);
    title('activity u');
    subplot(1, 3, 2);
    imagesc(1:fieldSize, 1:fieldSize, history_mu(:, :, i), [0 1]);
    title('memory trace u');
    subplot(1, 3, 3);
    imagesc(1:fieldSize, 1:fieldSize, history_v(:, :, i), [-10 10]);
    title('activity v');
    colormap(jet(256));
    drawnow;
    pause(0.01);
  end
end


% view cuts through field activities at one position iteratively
if 0
  xPos = 75;
  figure;
  for i = 1 : nStoredFields
    subplot(2, 1, 1);
    plot(1:fieldSize, zeros(1, fieldSize), ':k', 1:fieldSize, history_mu(:, xPos, i), '-c', ...
      1:fieldSize, history_u(:, xPos, i), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-10 10]);
    ylabel('activity u');
    subplot(2, 1, 2);
    plot(1:fieldSize, zeros(1, fieldSize), ':k', 1:fieldSize, history_v(:, xPos, i), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-10 10]);
    ylabel('activity v');
    pause(0.01);
  end
end


% view cuts through field activities as mesh plot over time
if 0
  yPos = 25;
  figure;
  subplot(3, 1, 1);
  mesh(1:fieldSize, 1:nStoredFields, squeeze(history_u(yPos, :, :))');
  zlabel('activity u');
  subplot(3, 1, 2);
  mesh(1:fieldSize, 1:nStoredFields, squeeze(history_mu(yPos, :, :))');
  zlabel('memory trace u');
  subplot(3, 1, 3);
  mesh(1:fieldSize, 1:nStoredFields, squeeze(history_v(yPos, :, :))');
  zlabel('activity v');
end


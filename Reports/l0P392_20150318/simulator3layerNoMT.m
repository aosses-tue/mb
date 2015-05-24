% 
% Offline simulator for Dynamic Neural Fields, 3-layer architecture
%
% The model consists of three fields:
% field_u - excitatory layer (usually perceptual field)
% field_v - shared inhibitory layer
% field_w - excitatory layer (usually working memory)
%
% Parameters for each layer have the corresponding letter as index.
% 
% Connections in the model:
% u->u (exc.), u->v (exc.), v->u (inh.),
% w->u (exc.), u->w (exc.)
% w->w (exc.), w->v (exc.), v->w (inh.),
% 
% Connection parameters and kernels between two layers are designated with
% the indices for both layers, with the layer that receives the input 
% written first.
%
% Example:
% kernel_uv: interaction kernel for connection v->u


%%%%%%%%%%%%%%%%%%%%
% model parameters %
%%%%%%%%%%%%%%%%%%%%

fieldSize = 181; % must be odd

tau_u = 20; tau_v = 5; tau_w = 20; % time constants of regular fields
beta_u = 5; beta_v = 5; beta_w = 5; % steepness parameters of sigmoid function
h_u = -5; h_v = -5; h_w = -5; % resting levels

% c: strength of interaction kernel
% sigma: width of interaction kernel
% g: strength of global inhibition
c_uu = 25; sigma_uu = 5;
c_uv = 20; sigma_uv = 10; g_uv = 0;
c_uw = 0; sigma_uw = 5;
c_vu = 10; sigma_vu = 5;
c_vw = 10; sigma_vw = 5;
c_wu = 20; sigma_wu = 5;
c_wv = 20; sigma_wv = 10; g_wv = 0;
c_ww = 25; sigma_ww = 5;
c_ws = 0; % strength of stimulus input to w layer (as fraction of input to u)

q_u = 0; q_v = 0; q_w = 0; % noise levels
sigma_q = 0; % width of the noise kernel


%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time course %
%%%%%%%%%%%%%%%%%%%%%%%%%%

nTrials = 1;
tMax = 250;

% tStoreFields = [100, 200]; % at these time steps the field activities will be stored
tStoreFields = 1:tMax; % use this to store field activities at all time steps
% tStoreFields = 1:5:tMax; % use this to store field activities at every 5th time step

% set times of stimulus presentation
% for multiple stimuli in one trial, repeat this for each one
% (e. g. tStimulusStart1 = ..., tStimulusStart2 = ...)
tStimulusStart = 50;
tStimulusEnd = 150;

% set up stimuli to be used during the simulation,
% for example:
stim1 = 6*gauss(1:fieldSize, 60, 5); % a localized input
stim2 = 6*gauss(1:fieldSize, 120, 10); % a broad localized input
stim3 = 2 * ones(1, fieldSize); % a homogeneous boost
stim4 = 5*gauss(1:fieldSize, 90, 5) + 1; % a localized input with homogeneous boost


%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

halfField = floor(fieldSize/2);

% create row vectors for field activities
field_u = zeros(1, fieldSize);
field_v = zeros(1, fieldSize);
field_w = zeros(1, fieldSize);

% create matrices to store field activities at different times
history_s = zeros(nTrials * length(tStoreFields), fieldSize); % stimulus
history_u = zeros(nTrials * length(tStoreFields), fieldSize);
history_v = zeros(nTrials * length(tStoreFields), fieldSize);
history_w = zeros(nTrials * length(tStoreFields), fieldSize);

% index of the current position in the history matrices
iHistory = 1;

% set up the interaction kernels
kernel_uu = c_uu * gaussNorm(-halfField:halfField, 0, sigma_uu);
kernel_uv = c_uv * gaussNorm(-halfField:halfField, 0, sigma_uv) + g_uv;
kernel_uw = c_uw * gaussNorm(-halfField:halfField, 0, sigma_uw);
kernel_vu = c_vu * gaussNorm(-halfField:halfField, 0, sigma_vu);
kernel_vw = c_vw * gaussNorm(-halfField:halfField, 0, sigma_vw);
kernel_wu = c_wu * gaussNorm(-halfField:halfField, 0, sigma_wu);
kernel_wv = c_wv * gaussNorm(-halfField:halfField, 0, sigma_wv) + g_wv;
kernel_ww = c_ww * gaussNorm(-halfField:halfField, 0, sigma_ww);

% set up the kernel for correlated noise (if required)
if sigma_q > 0
  kernel_q = gaussNorm(-halfField:halfField, 0, sigma_q);
end


%%%%%%%%%%%%%%
% simulation %
%%%%%%%%%%%%%%

% loop over trials
for i = 1 : nTrials
  
  % prepare matrix that holds the stimulus for each time step
  stimulus = zeros(tMax, fieldSize);
  
  % if needed, create a new stimulus pattern for every trial
  stimPos = round(fieldSize * rand); % random position
  stimRand = 6*gauss(1:fieldSize, stimPos, 5);
  
  % write the stimulus pattern into the stimulus matrix for all time steps
  % where it should be active
  for j = tStimulusStart : tStimulusEnd
    stimulus(j, :) = stimRand;
  end
  % repeat this if you want to use multiple stimuli in one trial
  
  
  % reset field activities to resting levels
  field_u(1:fieldSize) = h_u;
  field_v(1:fieldSize) = h_v;
  field_w(1:fieldSize) = h_w;
  
  % loop over time steps
  for t = 1 : tMax
    % calculation of field outputs
    output_u = sigmoid(field_u, beta_u, 0);
    output_v = sigmoid(field_v, beta_v, 0);
    output_w = sigmoid(field_w, beta_w, 0);
    
    % circular padding of outputs for convolution
    output_u_padded = [output_u(halfField+2:fieldSize), output_u(1:fieldSize), output_u(1:halfField)];
    output_v_padded = [output_v(halfField+2:fieldSize), output_v(1:fieldSize), output_v(1:halfField)];
    output_w_padded = [output_w(halfField+2:fieldSize), output_w(1:fieldSize), output_w(1:halfField)];

    % get endogenous input to fields by convolving outputs with interaction kernels
    conv_uu = conv2(1, kernel_uu, output_u_padded, 'valid');
    conv_vu = conv2(1, kernel_vu, output_u_padded, 'valid');
    conv_uv = conv2(1, kernel_uv, output_v_padded, 'valid');
    conv_wu = conv2(1, kernel_wu, output_u_padded, 'valid');
    conv_uw = conv2(1, kernel_uw, output_w_padded, 'valid');
    conv_ww = conv2(1, kernel_ww, output_w_padded, 'valid');
    conv_vw = conv2(1, kernel_vw, output_w_padded, 'valid');
    conv_wv = conv2(1, kernel_wv, output_v_padded, 'valid');
    
    % create field noise for this timestep
    noise_u = q_u * randn(1, fieldSize);
    noise_v = q_v * randn(1, fieldSize);
    noise_w = q_w * randn(1, fieldSize);
    if sigma_q > 0 % create spatially correlated noise by convolution
      noise_u_padded = [noise_u(halfField+2:fieldSize), noise_u, noise_u(:, 1:halfField)];
      noise_v_padded = [noise_v(halfField+2:fieldSize), noise_v, noise_v(:, 1:halfField)];
      noise_w_padded = [noise_w(halfField+2:fieldSize), noise_w, noise_w(:, 1:halfField)];
      
      noise_u = conv2(1, kernel_q, noise_u_padded, 'valid');
      noise_v = conv2(1, kernel_q, noise_v_padded, 'valid');
      noise_w = conv2(1, kernel_q, noise_w_padded, 'valid');
    end
    
    % update field activities
    field_u = field_u + 1/tau_u * (-field_u + h_u + stimulus(t, :) ...
      + conv_uu - conv_uv + conv_uw) + noise_u;
    field_v = field_v + 1/tau_v * (-field_v + h_v + conv_vu + conv_vw) + noise_v;
    field_w = field_w + 1/tau_w * (-field_w + h_w + c_ws * stimulus(t, :) ...
      + conv_wu - conv_wv  + conv_ww) + noise_w;
    
    % store field activities at the selected time steps
    if any(tStoreFields == t)
      history_s(iHistory, :) = stimulus(t, :);
      history_u(iHistory, :) = field_u;
      history_v(iHistory, :) = field_v;
      history_w(iHistory, :) = field_w;
      iHistory = iHistory + 1;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization of results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note: you can also access and plot all stored results afterwards from the
% Matlab prompt

% view plots of all stored field activities iteratively
nStoredFields = nTrials * length(tStoreFields);
if 1
  figure;
  for i = 1 : nStoredFields
    subplot(3, 1, 1);
    plot(0:fieldSize-1, zeros(1, fieldSize), ':k', 1:fieldSize, history_s(i, :), '--g', ...
      1:fieldSize, history_u(i, :), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-15 15]);
    zlabel('activity u');
    subplot(3, 1, 2);
    plot(0:fieldSize-1, zeros(1, fieldSize), ':k', 1:fieldSize, history_v(i, :), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-15 15]);
    zlabel('activity v');
    subplot(3, 1, 3);
    plot(0:fieldSize-1, zeros(1, fieldSize), ':k', 1:fieldSize, history_w(i, :), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-15 15]);
    zlabel('activity w');
    drawnow;
    pause(0.01);
  end
end

% view evolution of field activities in each trial as mesh plot
nFieldsPerTrial = length(tStoreFields);
if 0
  figure;
  for i = 1 : nTrials
    subplot(3, 1, 1);
    mesh(1:fieldSize, tStoreFields, history_u((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('activity u');
    subplot(3, 1, 2);
    mesh(1:fieldSize, tStoreFields, history_v((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('activity v');
    subplot(3, 1, 3);
    mesh(1:fieldSize, tStoreFields, history_w((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('activity w');
    pause;
  end
end


% view mesh plot of all stored field activities together
if 0
  figure;
  subplot(3, 1, 1);
  mesh(1:fieldSize, 1:nStoredFields, history_u(:, :));
  zlabel('activity u');
  subplot(3, 1, 2);
  mesh(1:fieldSize, 1:nStoredFields, history_v(:, :));
  zlabel('activity v');
  subplot(3, 1, 3);
  mesh(1:fieldSize, 1:nStoredFields, history_w(:, :));
  zlabel('activity w');
end




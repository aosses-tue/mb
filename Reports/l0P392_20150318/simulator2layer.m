% 
% Offline simulator for Dynamic Neural Fields, 2-layer architecture with
% memory trace
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

fieldSize = 181; % must be odd

tau_u = 20; tau_v = 5; % time constants of regular fields
tau_mu_build = 1000; tau_mu_decay = 5000; % time constants of memory traces
beta_u = 5; beta_v = 5; % steepness parameters of sigmoid function
h_u = -5; h_v = -5; % resting levels

% c: strength of interaction kernel
% sigma: width of interaction kernel
% g: strength of global inhibition
c_uu = 20; sigma_uu = 5;
c_uv = 15; sigma_uv = 10; g_uv = 0;
c_vu = 15; sigma_vu = 5;
c_mu = 0; sigma_mu = 5;

q_u = 0; q_v = 0; % noise levels
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
memTrace_u = zeros(1, fieldSize);
field_v = zeros(1, fieldSize);

% create matrices to store stimulus field activities at different times
history_s = zeros(nTrials * length(tStoreFields), fieldSize); % stimulus
history_u = zeros(nTrials * length(tStoreFields), fieldSize);
history_mu = zeros(nTrials * length(tStoreFields), fieldSize);
history_v = zeros(nTrials * length(tStoreFields), fieldSize);

% index of the current position in the history matrices
iHistory = 1;

% set up the interaction kernels
kernel_uu = c_uu * gaussNorm(-halfField:halfField, 0, sigma_uu);
kernel_uv = c_uv * gaussNorm(-halfField:halfField, 0, sigma_uv) + g_uv;
kernel_vu = c_vu * gaussNorm(-halfField:halfField, 0, sigma_vu);
kernel_mu = c_mu * gaussNorm(-halfField:halfField, 0, sigma_mu);

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
  
  % if needed, create a new stimulus pattern for every trial,
  % for example
  stimPos = round(fieldSize * rand); % random position
  stimRand = 6*gauss(1:fieldSize, stimPos, 5);
  
  % write the stimulus pattern into the stimulus matrix for all time steps
  % where it should be active
  for j = tStimulusStart : tStimulusEnd
    stimulus(j, :) = stim1;
  end
  % repeat this if you want to use multiple stimuli in one trial
  
  
  % reset field activities to resting levels
  field_u(1:fieldSize) = h_u;
  field_v(1:fieldSize) = h_v;
  
  % loop over time steps
  for t = 1 : tMax
    % calculation of field outputs
    output_u = sigmoid(field_u, beta_u, 0);
    output_v = sigmoid(field_v, beta_v, 0);
    
    % circular padding of outputs for convolution
    output_u_padded = [output_u(halfField+2:fieldSize), output_u(1:fieldSize), output_u(1:halfField)];
    output_v_padded = [output_v(halfField+2:fieldSize), output_v(1:fieldSize), output_v(1:halfField)];
    memTrace_u_padded = [memTrace_u(halfField+2:fieldSize), memTrace_u, memTrace_u(:, 1:halfField)];

    % get endogenous input to fields by convolving outputs with interaction kernels
    conv_uu = conv2(1, kernel_uu, output_u_padded, 'valid');
    conv_vu = conv2(1, kernel_vu, output_u_padded, 'valid');
    conv_uv = conv2(1, kernel_uv, output_v_padded, 'valid');
    conv_mu = conv2(1, kernel_mu, memTrace_u_padded, 'valid');
    
    % create field noise for this timestep
    noise_u = q_u * randn(1, fieldSize);
    noise_v = q_v * randn(1, fieldSize);
    if sigma_q > 0 % create spatially correlated noise by convolution
      noise_u_padded = [noise_u(halfField+2:fieldSize), noise_u, noise_u(:, 1:halfField)];
      noise_v_padded = [noise_v(halfField+2:fieldSize), noise_v, noise_v(:, 1:halfField)];
      
      noise_u = conv2(1, kernel_q, noise_u_padded, 'valid');
      noise_v = conv2(1, kernel_q, noise_v_padded, 'valid');
    end
    
    % update field activities
    field_u = field_u + 1/tau_u * (-field_u + h_u + stimulus(t, :) ...
      + conv_uu - conv_uv + conv_mu) + noise_u;
    field_v = field_v + 1/tau_v * (-field_v + h_v + conv_vu) + noise_v;
    
    % update memory trace (only if there is activity in the field)
    activeRegions_u = field_u > 0;
    if any(activeRegions_u)
      memTrace_u = memTrace_u + 1/tau_mu_build * (-memTrace_u + output_u) .* activeRegions_u ...
        + 1/tau_mu_decay * (-memTrace_u) .* (1-activeRegions_u);
    end
    
    % store field activities at the selected time steps
    if any(tStoreFields == t)
      history_s(iHistory, :) = stimulus(t, :);
      history_u(iHistory, :) = field_u;
      history_mu(iHistory, :) = memTrace_u;
      history_v(iHistory, :) = field_v;
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
    subplot(2, 1, 1);
    plot(0:fieldSize-1, zeros(1, fieldSize), ':k', 1:fieldSize, history_s(i, :), '--g', ...
      1:fieldSize, history_mu(i, :) + h_u, '-c', 1:fieldSize, history_u(i, :), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-10 10]);
    ylabel('activity u');
    subplot(2, 1, 2);
    plot(0:fieldSize-1, zeros(1, fieldSize), ':k', 1:fieldSize, history_v(i, :), '-b');
    set(gca, 'XLim', [0 fieldSize-1], 'YLim', [-10 10]);
    ylabel('activity v');
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
    mesh(1:fieldSize, tStoreFields, history_mu((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('memory trace u');
    subplot(3, 1, 3);
    mesh(1:fieldSize, tStoreFields, history_v((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('activity v');
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
  mesh(1:fieldSize, 1:nStoredFields, history_mu(:, :));
  zlabel('memory trace u');
  subplot(3, 1, 3);
  mesh(1:fieldSize, 1:nStoredFields, history_v(:, :));
  zlabel('activity v');
end




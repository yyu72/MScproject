function [xx_2d_pos,time_points]=WF_single_simscp_time(N,mu,sp,x0,T,ns,dt, plotfig)% Parameters
%WF_sims_single function does simulation to create frequency array of genome,time
% Initialize allele frequency array
x = zeros(T+1, 1);
if sp < 0
    x(1,:) = mu/abs(sp);
else
    x(1,:) = x0;
end
stopt=[];
for i = 1 : T
    % Calculate selection terms
    sel_pos = sp * x(i,1) .* (1 - x(i,1)); % positive selection
    
    % Calculate mutation terms
    mut_fwd = mu * (1 - x(i,1)); % forward mutation
    mut_bwd = mu * x(i,1); % backward mutation
    
    % Calculate frequency in the next generation due to selection and mutation
    x(i+1,1) = x(i,1) + sel_pos + mut_fwd - mut_bwd;
    
    % Apply genetic drift
    x(i+1,1) = binornd(N, x(i+1,1)) / N;
    
    % Check if x has reached 1
    if x(i+1,1) >= 1
        
        break; % Exit the loop
    end
    stopt= i;
end




tt=(0:T);
time_points=tt(1:dt:stopt);

%size(tt)
%size(x)

if plotfig == 1
figure;
plot(tt,x','k-')
hold on
end
%sample through every dt generations,than do binomial sampling for 10
%samples through samplex
samplex_pos=x(1:dt:stopt,:);

if plotfig == 1
plot(time_points,samplex_pos,'rx')
end

num_sites_pos = size(samplex_pos, 2);
num_timepoints = numel(time_points);

x2d_pos = zeros(num_sites_pos, num_timepoints);

xx_pos = zeros(num_sites_pos, num_timepoints);

% Iterate over the sampled frequency array
for i = 1:num_sites_pos
    for k = 1:num_timepoints
        % Sample the binomial distribution for allele 1
        x2d_pos(i, k) = samplex_pos(k, i);
        %uncomment this to cancel sampling
        xx_pos(i, k) = binornd(ns, x2d_pos(i, k)) / ns;
    end
end

xx_2d_pos = reshape(xx_pos, num_sites_pos, num_timepoints);
%use samplex to use unsampled array, xx_pos for sampled array
if plotfig == 1
plot(time_points,xx_2d_pos,'go')
hold off
end

%[Nopt,sopt,muopt,NCI,sCI,muCI,d0,alphaopt,thetaopt]=GenomeScanMutation(time_points,xx_2d_pos,I,ns,Dt,N0,mu0,samplcorrections,initialparams,sinitopt,null,fine, prior, plotfig,suppressoutput,CI,boundary,optimmethod);
end
%%% c 2002 Andrew I. Schein
%%% lpca_fitwbias.m -- code for fitting the logistic PCA model

%%% I wrote this code based off of the manuscript 
% A Generalized Linear Model for Principal Component
% Analysis of Binary Data
% to be presented at AI & Statistics 2003

%%% The function name 'wbias' is to indicate that this code includes
%%% a bias term.  The code differs from the manuscript 
%%% presentation in that the bias vector is encoded as an extra 
%%% (L+1) row in the V matrix. 

function [U,V] = lpca_fitwbias(X,L,numiterations,threshold) 

  %% INPUTS:
  %% X: the binary data matrix, rows are observations
  %% L: the size of the desired latent space, e.g. 5
  %% numiterations: the maximum number of iterations
  %% when set to 0, we use a threshold change in LL instead
  %% to determine when to stop

  %% OUTPUTS:
  %% U: the coordinates in the latent space
  %% V: the basis of the latent space

  if nargin < 3, numiterations = 0; end;
  if nargin < 4, threshold = 0.5; end; %% if we use a threshold for
                                   %stopping,
				   %this is it.
  
  [N,D] = size(X);   %% Get the dimensions of the data
  U = rand(N,L+1);    %% ...so we can create U and V matrices
  V = rand(L+1,D);
  U(:,L+1) = ones(N,1);  %% This is for the bias vector.
   
  oldU = U;
  oldV = V;

  iter = 0;          %% Used to count number of iterations in loop.
  changel = 1000;    %% keep track of the change in LL between iterations;
  oldLL = 0; %% keep track of the last LL (init to -Infinity)
  stop = 0;       %% Conditions for stopping the model-fitting

%%%
  
  i0 = find(X == 0);  % To prevent constant * NaN errors we separate
  i1 = find(X == 1);  %% 1's from 0's for likelihood calculations.
  

%%% Get a starting log likeilood.
  z = U*V;   
  oldLL = -sum(sum(log(1+exp(-z(i0))))) - sum(sum(log(1+exp(z(i1)))));
  fprintf('Iterations \t log_likelihood \n');
  fprintf('%d \t\t %f \n',iter,oldLL);

%%% The variable names below are quite similar to the Publication.  
  
while (stop <= 0)
  iter = iter + 1;
 
%%% Begin U Update
  z = U*V;
  T = tanh(z / 2)  ./ z;   
  b = (2*X - 1) * V(1:L,:)' - (T .* repmat(V(L+1,:),N,1) * V(1:L,:)') ; 
  for (n = 1:N) 
    A = (repmat(T(n,:),L,1) .* V(1:L,:)) * V(1:L,:)';  
    U(n,1:L) =  (A \ b(n,:)')' ;
  end;	

%%% Begin V Update
  z = U*V;
  T = tanh(z / 2)  ./ z;   
  b =  (2*X' - 1) * U ;
  for (d = 1:D)
    A = (repmat(T(:,d),1,L+1) .* U )' * U;         
    V(:,d) = (A \ b(d,:)');
  end

%%% Report the Log Likelihood
  z = U*V;  
  newLL = -sum(sum(log(1+exp(z(i0))))) - sum(sum(log(1+exp(-z(i1)))));

  change = newLL - oldLL;
  oldLL = newLL;
  fprintf('%d \t\t %f \n',iter,newLL);

%%% Run some checks to see if we go through the loop again.

  if (numiterations > 0) 
    if (iter >= numiterations) 
      stop = 1;
    end;
  else if (change < threshold) 
      stop =1;
       end; 

  end;	

%%% Catch numerical errors
if (isnan(newLL) )
   stop = 1;
   U = oldU;
   V = oldV;

else if (stop == 0)
    oldU = U;
    oldV = V;
  end;

end;
end;

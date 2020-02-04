
function [c, c1, c2, w1,w2] = SCCA_SCCA(X1,X2, Spar1,Spar2,niter,w1_sign,w2_sign)

%% this function calculates sparse CCA weights between two sets of variables
%
% Inputs:
% X1: first view of the data
% X2: second view of the data
% Spar1: first sparsity constraint
% Spar2: second sparsity constraint
% niter: number of iterrations script must run through
% w1_sign: sign of canonical weights associated with X1, 1 = positive, -1 =
% negative, 0 = unconstrained
% w2_sign: sign of canonical weights assocated with X2, 1 = positive, -1 =
% negative, 0 = unconstrained
%
% Outputs:
% c1: first canonical correlate
% c2: second canonical correlate
% w1: weights associated with first canonical correlate
% w2: weights associated with second canonical correlate
% 

%% initate the first canoncial weight vector

[w1s,s,w2s] = svd((X1'*X2),'econ');
w1 = w1s(:,1)/norm(w1s(:,1),2);

%% set the sparsity values

dimX1 = size(X1,2);
dimX2 = size(X2,2);

r1 = round((sqrt(dimX1)*Spar1),2);
r2 = round((sqrt(dimX2)*Spar2),2);


    if (r1 <=1)
    r1 = 1;    
    else    
    end

%%

        for n = 1:niter;
        n
        if (w2_sign == 0);
        a2 = ((X2')*X1)*w1;         
        elseif (w2_sign == 1);        
        a2 = subplus(((X2')*X1)*w1);
        else    
        a2 = -subplus(-(((X2')*X1)*w1));  
        end

        sygna2 = (((a2 > 0)*2) - 1);    
        S2_p = abs(a2);
        S2 = sygna2.*(S2_p.*(S2_p > 0));
        a2 = S2;
        w2 = a2/norm(a2,2);
        norm_wr = round(norm(w2,1),2);  

        temp = 1;
        if  (norm_wr >= r2);   
        deltamax = max(abs(a2))/2;   
        delta_temp = 0;    
          while (norm_wr ~= r2);
        temp = temp+1;
        delta_sign = (((norm(w2,1) > r2)*2) -1);    
        delta = delta_temp + ((delta_sign)*(deltamax));  
        S2_p = ((abs(a2))-delta);
        S2 = sygna2.*(S2_p.*(S2_p > 0));
        w2 = S2/norm(S2,2);
        norm_wr = round(norm(w2,1),2);  
        delta_temp = delta;
        deltamax = deltamax/2;
          end
        else  
        end
 
        if (w1_sign == 0);
            a1 = ((X1')*X2)*w2;
        elseif (w1_sign == 1);
            a1 = subplus(((X1')*X2)*w2);
        else
            a1 = -subplus(-((X1')*X2)*w2);
        end

        sygna1 = (((a1 > 0)*2) - 1);
        S1_p = abs(a1);
        S1 = sygna1.*(S1_p.*(S1_p > 0));
        a1 = S1;
        w1 = a1/norm(a1,2);
        norm_wr = round(norm(w1,1),2);
        
        temp = 1;
        if  (norm_wr >= r1);
            deltamax = max(abs(a1))/2;
            delta_temp = 0;
            while (norm_wr ~= r1);
                temp = temp+1;
                delta_sign = (((norm(w1,1) > r1)*2) -1);
                delta = delta_temp + ((delta_sign)*(deltamax));
                S1_p = ((abs(a1))-delta);
                S1 = sygna1.*(S1_p.*(S1_p > 0));
                w1 = S1/norm(S1,2);
                norm_wr = round(norm(w1,1),2);
                delta_temp = delta;
                deltamax = deltamax/2;
            end
        else
        end
        
        end

c1 = X1*w1;
c2 = X2*w2;
c = corr(c1,c2);


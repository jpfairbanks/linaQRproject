% Domenic Carr & James Fairbanks
% MATH6643 Project

% Top-level Script for Testing Speed of Updates
clear all

%%%% Parameters to Change
% mvals = [100 500 1000];
mvals = [500];
counter = 1;

        
% create data storage matrix
data = zeros(150,12);

for i = 1:1:length(mvals)
    m = mvals(i);
    steps = floor(log2(m))-1 + 1;
    
    % create appropriate n vector
    nvals = ones(1,steps+1);
    for j=2:1:steps
        nvals(j) = 2^(j-1);
    end
    nvals = 2*nvals;
    nvals(steps+1) = m;
    
    
    for ii = 1:1:length(nvals)        
        n = nvals(ii);
        stopval = n/2;
        steps = floor(log2(n)) -1 + 1;
    
        % create appropriate k vector
        if m ~= n
            kvals = ones(1,steps);
            for j = 2:1:steps
                kvals(j) = 2^(j-1);
            end
        else
            kvals = ones(1,steps+1);
            for j = 2:1:steps
                kvals(j) = 2^(j-1);
            end
            kvals(steps+1) = n/4;            
        end
        
        tvFull = zeros(30,1);
        tvEager = zeros(30,1);
        tvLazy = zeros(30,1);
        
        for iii = 1:1:length(kvals)            
            k = kvals(iii);
            columnoffset = 1;
            
            % run <m-n-k> combination 30 times            
            for jj = 1:1:30
                [tFull,tEager,tLazy] = ExecuteSpeedTest(m,n,k);
                tvFull(jj) = tFull;
                tvEager(jj) = tEager;
                tvLazy(jj) = tLazy;                
            end
            
            % store m-n-k values into data matrix
            data(counter,columnoffset) = m;
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = n;
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = k;
            columnoffset = columnoffset+1;
            
            % calculate median run times and store in matrix
            data(counter,columnoffset) = median(tvFull);
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = median(tvEager);
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = median(tvLazy);
            columnoffset = columnoffset+1;
            
            % perform signed rank test                  
            [p,h] = signrank(tvFull,tvEager);
            data(counter,columnoffset) = p;
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = h;
            columnoffset = columnoffset+1;
            
            [p,h] = signrank(tvFull,tvLazy);
            data(counter,columnoffset) = p;
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = h;
            columnoffset = columnoffset+1;
            
            [p,h] = signrank(tvLazy,tvEager);
            data(counter,columnoffset) = p;
            columnoffset = columnoffset+1;
            data(counter,columnoffset) = h;
            
            display([m n k]);
            counter = counter + 1;
            
        end
        
    end
    
end


data = data(1:counter,:);
csvwrite('workhourseoutput.csv',data);




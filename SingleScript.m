% Domenic Carr & James Fairbanks
% MATH6643 Project

% Top-level Script for Testing Speed of Updates for particular m,n,k
% combination
clear all

%%%% Parameters to Change

mnk_data = [1000 16 8; 1000 32 16; 1000 64 32; 1000 128 64; 1000 256 128; 1000 512 256; 1000 1000 256; 1000 1000 500];
% m = 500;
% n = 16;
% k = 8;

% Start Work
data = zeros(8,12);

for i = 1:1:length(mnk_data(:,1))
    
    m = mnk_data(i,1);
    n = mnk_data(i,2);
    k = mnk_data(i,3);
    
    tvFull = zeros(30,1);
    tvEager = zeros(30,1);
    tvLazy = zeros(30,1);
    
    for jj = 1:1:30
        [tFull,tEager,tLazy] = ExecuteSpeedTest(m,n,k);
        tvFull(jj) = tFull;
        tvEager(jj) = tEager;
        tvLazy(jj) = tLazy;
    end

    columnoffset = 1;

    % store m-n-k values into data matrix
    data(i,columnoffset) = m;
    columnoffset = columnoffset+1;
    data(i,columnoffset) = n;
    columnoffset = columnoffset+1;
    data(i,columnoffset) = k;
    columnoffset = columnoffset+1;

    % calculate median run times and store in matrix
    data(i,columnoffset) = median(tvFull);
    columnoffset = columnoffset+1;
    data(i,columnoffset) = median(tvEager);
    columnoffset = columnoffset+1;
    data(i,columnoffset) = median(tvLazy);
    columnoffset = columnoffset+1;

    % perform signed rank test                  
    [p,h] = signrank(tvFull,tvEager);
    data(i,columnoffset) = p;
    columnoffset = columnoffset+1;
    data(i,columnoffset) = h;
    columnoffset = columnoffset+1;

    [p,h] = signrank(tvFull,tvLazy);
    data(i,columnoffset) = p;
    columnoffset = columnoffset+1;
    data(i,columnoffset) = h;
    columnoffset = columnoffset+1;

    [p,h] = signrank(tvLazy,tvEager);
    data(i,columnoffset) = p;
    columnoffset = columnoffset+1;
    data(i,columnoffset) = h;
    
    display([m n k]);

end

csvwrite('output_remainingcases.csv',data);




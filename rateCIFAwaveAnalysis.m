%%%
clear;clc;

fig_hei=0.14;
fig_wei=0.8;
lef_cor_x=0.1;
lef_cor_y=0.1;

title_size=14;
label_size=14;
font_size=14;
legend_size = 12;

%==================================global para========
fs =10;
N = 2048;
NTime = 10;

matFile = 'rateCIFA.mat';
txtFile = 'rate-CIFA.txt';

if exist(matFile,'file')
    % delete(matFile);
    load (matFile);
else
    disp(['no ', matFile,'... checking rateCIFA.txt...']);
    if exist(txtFile,'file')
        disp([txtFile, ' exist...']);
        load (txtFile);
        rateCIFA = rate_CIFA';
        fs = 1/rateCIFA(1);
        rateCIFA = rateCIFA(2,:);
        save rateCIFA rateCIFA;
    else
        disp([txtFile, ' no exist...return;']);
        return;
    end
end



waveName = 'haar';
sig = rateCIFA(:,1:N); %[2 8 4 9 1 2 1 1 8 5 5 1 8 6 3 5 2 8 4 9 1 2 1 1 8 5 5 1 8 6 3 5];
[C,L] = wavedec(sig,5,waveName);
% app Coef
cA5 = appcoef(C,L,'db4',5);
% detail Coef
[cD1,cD2,cD3,cD4,cD5] = detcoef(C,L,[1,2,3,4,5]);
cD1 = abs(cD1);
cD2 = abs(cD2);
cD3 = abs(cD3);
cD4 = abs(cD4);
cD5 = abs(cD5);

A5 = wrcoef('a',C,L,waveName,5);
D1 = wrcoef('d',C,L,waveName,1);
D2 = wrcoef('d',C,L,waveName,2);
D3 = wrcoef('d',C,L,waveName,3);
D4 = wrcoef('d',C,L,waveName,4);
D5 = wrcoef('d',C,L,waveName,5);


Xa = D3 + D4 + D5;   % attck
Xb = D2 + D1 + A5;   % back

% ==============================figure(1) decomp=============begin========
x0 = (2^0:2^0:length(sig))/fs;
x1 = (2^1:2^1:length(sig))/fs;
x2 = (2^2:2^2:length(sig))/fs;
x3 = (2^3:2^3:length(sig))/fs;
x4 = (2^4:2^4:length(sig))/fs;
x5 = (2^5:2^5:length(sig))/fs;

figure(1);
subplot(6,1,1,'position', [lef_cor_x lef_cor_y+5*fig_hei fig_wei fig_hei]);
stem(x0,sig,'.');
axis([1/fs,N/fs,-inf,inf]);
ylabel('origin sig');
set(gca, 'XTick',[]);
title('wavelet decomposation');

subplot(6,1,2,'position', [lef_cor_x lef_cor_y+4*fig_hei fig_wei fig_hei]);
stem(x1,cD1,'.');
axis([1/fs,N/fs,1,inf]);
ylabel('level 1');
set(gca, 'XTick',[]);

subplot(6,1,3,'position', [lef_cor_x lef_cor_y+3*fig_hei fig_wei fig_hei]);
stem(x2,cD2,'.');
axis([1/fs,N/fs,1,inf]);
ylabel('level 2');
set(gca, 'XTick', []);

subplot(6,1,4,'position', [lef_cor_x lef_cor_y+2*fig_hei fig_wei fig_hei]);
stem(x3,cD3,'.');
axis([1/fs,N/fs,1,inf]);
ylabel('level 3');
set(gca, 'XTick', []);

subplot(6,1,5,'position', [lef_cor_x lef_cor_y+1*fig_hei fig_wei fig_hei]);
stem(x4,cD4,'.');
axis([1/fs,N/fs,1,inf]);
ylabel('level 4');
set(gca, 'XTick', []);

subplot(6,1,6,'position', [lef_cor_x lef_cor_y+0*fig_hei fig_wei fig_hei]);
stem(x5,cD5,'.');
axis([1/fs,N/fs,1,inf]);
ylabel('level 5');
% ==============================figure(1) decomp=============end========

% ==============================figure(2) reconstr==========begin========
figure(2);
subplot(3,1,1,'position', [lef_cor_x lef_cor_y+ 2*0.27 fig_wei 0.27]);
plot(x0,sig);
axis([1/fs,N/fs,-inf,inf]);
ylabel('origin sig');
set(gca, 'XTick', []);
title('wavelet reconstruction');

subplot(3,1,2,'position', [lef_cor_x lef_cor_y+ 1*0.27 fig_wei 0.27]);
plot(x0,Xa);
axis([1/fs,N/fs,-inf,inf]);
ylabel('attack Xa');
set(gca, 'XTick', []);

subplot(3,1,3,'position', [lef_cor_x lef_cor_y  fig_wei 0.27]);
plot(x0,Xb);
axis([1/fs,N/fs,-inf,inf]);
ylabel('background Xb');
% ==============================figure(2) reconstr=============end=======

% ==============================figure(3) statistics Xa(n)======end=======
% N = 10s; winlen = N*fs = 100
winlen = NTime*fs;
Column = floor(length(sig)/winlen);
Xa = abs(Xa(1:Column*winlen));
Xa = reshape(Xa, winlen, Column);

avgXa = mean(Xa);
Xtime = NTime:NTime:length(sig)/fs;
figure(3);
plot(Xtime, avgXa, '-*');
title('avg |Xa(n)| NT = 10s');
% ==============================figure(3) statistics Xa(n)======end=======

% Perform a single-level decomposition of the signal===begin========
% sig = [9 7 3 5];
% l_s = length(sig);
% [cA1, cD1] = dwt (sig,'haar');
% A1 = upcoef('a',cA1,'haar',1,l_s);
% D1 = upcoef('d',cD1,'haar',1,l_s);
% A0 = idwt(cA1,cD1,'haar',l_s);
% Perform a single-level decomposition of the signal===end===========
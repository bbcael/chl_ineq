clear all; close all; clc; load chl25.mat; % clear workspace & load data

%% STEP 1 create 15 versions of dataset: 1°/2°/4°, Southern/Atlantic/Indo-Pacific, Winter/Spring/Summer/Fall, 0-60N/0-60S/0-40N/0-40S/40-60N

for i = 1:240; % convert structure to array
    c(:,:,i) = data(i+6).chl; % grab each month as a new entry in 3rd dimension of chl array
end
clear i data; % clear extra variables
for i = 1:20; % convert 1° monthly to 1°/2°/4° annual
    C1(:,:,i) = nanmean(c(:,:,((i-1).*12+1):(12.*i)),3); % take NaN-ignoring mean of each month for each year
    dwi = nanmean(c(:,:,((i-1).*12+1):(12.*i-9)),3); % winter only
    dsp = nanmean(c(:,:,((i-1).*12+4):(12.*i-6)),3); % spring only
    dsu = nanmean(c(:,:,((i-1).*12+7):(12.*i-3)),3); % summer only
    dfa = nanmean(c(:,:,((i-1).*12+10):(12.*i)),3); % fall only
    d = squeeze(C1(:,:,i)); % take average of 2x2 squares to get to 2°
    d = (d(1:2:end,:)+d(2:2:end,:))./2;
    d = (d(:,1:2:end)+d(:,2:2:end))./2;
    dwi = (dwi(1:2:end,:)+dwi(2:2:end,:))./2;
    dwi = (dwi(:,1:2:end)+dwi(:,2:2:end))./2;
    dsp = (dsp(1:2:end,:)+dsp(2:2:end,:))./2;
    dsp = (dsp(:,1:2:end)+dsp(:,2:2:end))./2;
    dsu = (dsu(1:2:end,:)+dsu(2:2:end,:))./2;
    dsu = (dsu(:,1:2:end)+dsu(:,2:2:end))./2;
    dfa = (dfa(1:2:end,:)+dfa(2:2:end,:))./2;
    dfa = (dfa(:,1:2:end)+dfa(:,2:2:end))./2;
    C2(:,:,i) = d;
    Cwi(:,:,i) = dwi;
    Csp(:,:,i) = dsp;
    Csu(:,:,i) = dsu;
    Cfa(:,:,i) = dfa;
    d = squeeze(C2(:,:,i)); % take average of 2x2 squares to get to 4°
    d = (d(1:2:end,:)+d(2:2:end,:))./2;
    d = (d(:,1:2:end)+d(:,2:2:end))./2;
    C4(:,:,i) = d;
end
clear i d dwi dsp dsu dfa; % clear extra variables
N = sum(~isnan(C4),3); % find how many years of data for each grid cell
N(N<20) = 0; N(N==20) = 1; % set mask: yes if & only if value for all 20 years
N(40,12) = 0; N(4:8,23:29) = 0; N(13,79) = 0; N(22,76) = 0; N(13:14,47:50) = 0; N(4:5,52:61) = 0;  % remove inland waters, coasts, and cryosphere -- values chosen by eye based on disconnection from open ocean
N2 = imresize(N,2); N2(N2>.5) = 1; N2(N2<.5) = 0; % get 2° and 1° versions
N1 = imresize(N2,2); N1(N1>.5) = 1; N1(N1<.5) = 0;
for i = 1:20; % now mask out land areas from each product
    d = squeeze(C1(:,:,i));
    d(N1==0) = NaN;
    C1(:,:,i) = d;
    d = squeeze(Cwi(:,:,i));
    d(N2==0) = NaN;
    Cwi(:,:,i) = d;
    d = squeeze(Csp(:,:,i));
    d(N2==0) = NaN;
    Csp(:,:,i) = d;
    d = squeeze(Csu(:,:,i));
    d(N2==0) = NaN;
    Csu(:,:,i) = d;
    d = squeeze(Cfa(:,:,i));
    d(N2==0) = NaN;
    Cfa(:,:,i) = d;
    d = squeeze(C2(:,:,i));
    d(N2==0) = NaN;
    C2(:,:,i) = d;    
    d = squeeze(C4(:,:,i));
    d(N==0) = NaN;
    C4(:,:,i) = d;
end
clear d i c;
Cp = C2; Cp(65:end,:,:) = NaN; Cp(:,56:100,:) = NaN; Cp(27:40,45:55,:) = NaN;  Cp(33:34,43:44,:) = NaN; % remove southern and atlantic oceans to get pacific-only
Ca = C2; Ca(65:end,:,:) = NaN; Ca(:,101:end,:) = NaN; Ca(:,1:42,:) = NaN; Ca(41:64,43:54,:) = NaN; Ca(39:40,43:44,:) = NaN; % remove indo-pacific and southern oceans to get atlantic-only
Cs = C2; Cs(1:64,:,:) = NaN; %  remove indo-pacific and atlantic oceans to get southern-only
C6n = C2; C6n([1:15 46:end],:,:) = NaN; % remove 0-90S and 60-90N to get 0-60N-only
C6s = C2; C6s([1:45 76:end],:,:) = NaN; % remove 0-90N and 60-90S to get 0-60S-only
Cm = C2; Cm([1:15 26:end],:,:) = NaN; % remove 0-40 and 60-90 to get 40-60-only
C4n = C2; C4n([1:25 46:end],:,:) = NaN; % 0-40N
C4s = C2; C4s([1:45 66:end],:,:) = NaN; % 0-40S

cwi = Cwi; csp = Csp; csu = Csu;  cfa = Cfa; % swap S-hemisphere seasons
cwi(46:90,:,:) = Csu(46:90,:,:); Cwi = cwi;
csp(46:90,:,:) = Cfa(46:90,:,:); Csp = csp;
csu(46:90,:,:) = Cwi(46:90,:,:); Csu = csu;
cfa(46:90,:,:) = Csp(46:90,:,:); Cfa = cfa;
clear cwi csp csu cfa;

%% STEP 2 get SST product to match

T = ncread('sst.mon.mean.nc','sst'); % load SST data
T = T(:,:,1837:2076); % keep only relevant months
T = permute(T,[2 1 3]); % change dimension to match Chl
T(T>40) = NaN; % make missing value flag a NaN
t = T; t(:,181:360,:) = T(:,1:180,:); t(:,1:180,:) = T(:,181:360,:); T = t; clear t; % translate to match Chl
for i = 1:240; % same calculation as above to make 2°
    d = squeeze(T(:,:,i));
    d = (d(1:2:end,:)+d(2:2:end,:))./2;
    d = (d(:,1:2:end)+d(:,2:2:end))./2;
    t(:,:,i) = d;
end
clear i d; T = t; clear t;
for i = 1:20; % csame calculation as above to make annual
    t(:,:,i) = nanmean(T(:,:,((i-1).*12+1):(12.*i)),3);
end
clear i; T = t; clear t;

NT = sum(~isnan(T),3); % same calculation as above for SST land mask
NT(NT>0) = 1;
C2([13 18],[59 86],:) = NaN; % remove Chl data for which there aren't SST data
C6n(18,86,:) = NaN; % same for relevant alternate Chl's
Ca([13 18],[59 86],:) = NaN;
Cm(18,86,:) = NaN;
clear NT;

for i = 1:20; % same calculation as above to remove land etc from SST
    d = squeeze(T(:,:,i));
    d(N2==0) = NaN;
    T(:,:,i) = d;    
end
clear i d;

%% STEP 3 compute trends and means

t = -9.5:9.5; % year vector

% here replace C2 with C1 or C4/Cwi or Csp or Csu or Cfa/C4n or C4s/C6n or
 % C6s/Cs or Cp or Ca or Cm for sensitivity tests; 
 % change i+j for C1 or C4 & comment out T regression if irrelevant
    %N2 = N; % 4°
    %N2 = N1; % 1°
for i = 1:90; % loop through the regressions
    for j = 1:180;  % to get trends, means, & uncertainties
        if N2(i,j)==1;
        mdl = fitlm(t,squeeze(Cm(i,j,:)));
        t2(i,j) = table2array(mdl.Coefficients(2,1));
        m2(i,j) = table2array(mdl.Coefficients(1,1));
        t2u(i,j) = table2array(mdl.Coefficients(2,2));
        m2u(i,j) = table2array(mdl.Coefficients(1,2));
        mdl = fitlm(t,squeeze(T(i,j,:)));
        tT(i,j) = table2array(mdl.Coefficients(2,1));
        mT(i,j) = table2array(mdl.Coefficients(1,1));
        tTu(i,j) = table2array(mdl.Coefficients(2,2));
        mTu(i,j) = table2array(mdl.Coefficients(1,2));
        end
    end
    i % print progress
end
clear i j mdl t;
N2 = N2(1:80,:); % no ocean in Antarctica so no regressions
    %N2 = N2(1:40,:); % 4°
    %N2 = N2(1:160,:); % 1°
t2(N2==0) = NaN; t2u(N2==0) = NaN; m2(N2==0) = NaN; m2u(N2==0) = NaN; % NaNs where it's land
tT(N2==0) = NaN; tTu(N2==0) = NaN; mT(N2==0) = NaN; mTu(N2==0) = NaN;
clear N2;
t2 = t2./m2; % put Chl slopes in relative change per year
t2u = abs(t2./m2).*sqrt((t2u./t2).^2+(m2u./m2).^2); % convert and propagate uncertainty
m2u = abs(m2u./m2); % convert and propagate uncertainty in mean too
m2 = log(m2); % put Chl means in log terms

%% STEP 4 Make some maps

x = -179:2:179; x(1) = -180; x(end) = 180; x = repmat(x',1,80)'; y = 89:-2:-69; y = repmat(y',1,180);
figure;
axesm robinson
axesm('robinson','MapLatLimit',[-90 90],'MapLonLimit',[20 380])
framem off
axis off
h = contourfm(y,x,m2./log(10),log10([.01 .0173 .03 .0548 .1 .173 .3 .548 1 1.73 3 5.48 10]));
hold on;
colormap(cbrewer2('BuGn'))
c = colorbar;
set(c,'fontsize',24,'ytick',log10([.01 .03 .1 .3 1 3 10]),'yticklabel',{'0.01','0.03','0.1','0.3','1','3','10'},'ticklabelinterpreter','latex')
caxis([-1.81 .52])
hold on;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,'FaceColor',[.9 .85 .8]);
title('Mean Chl [(mg/m$^{3}$)]','fontsize',24,'interpreter','latex')

figure;
axesm robinson
axesm('robinson','MapLatLimit',[-90 90],'MapLonLimit',[20 380])
framem off
axis off
h = contourfm(y,x,1200.*nanconv(t2,fspecial('gaussian',5,2.5)),[-20:5:30]);
hold on;
colormap(flipud(cbrewer2('RdBu')))
c = colorbar;
set(c,'fontsize',24,'ytick',[-20:10:30],'ticklabelinterpreter','latex')
caxis([-35 35])
hold on;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,'FaceColor',[.9 .85 .8]);
title('Chl Trend [\%/decade]','fontsize',24,'interpreter','latex')

figure;
axesm robinson
axesm('robinson','MapLatLimit',[-90 90],'MapLonLimit',[20 380])
framem off
axis off
h = contourfm(y,x,10.*double(tT),[-1:.1:1]);
hold on;
colormap(cbrewer2('PuOr'))
c = colorbar;
set(c,'fontsize',24,'ytick',[-.5:.3:.9],'ticklabelinterpreter','latex')
caxis([-.9 .9])
hold on;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,'FaceColor',[.9 .85 .8]);
title('SST Trend [$^\circ$C/decade]','fontsize',24,'interpreter','latex')

%% STEP 5 test green-greener-blue-bluer hypothesis

w = 89:-2:-69; w = repmat(w,180,1)'; w(isnan(m2)) = NaN; w = cosd(w); w = w.*sum(~isnan(w(:)))./nansum(w(:)); % latitude weighting
    %w = 88:-4:-68; w = repmat(w,90,1)'; w(isnan(m2)) = NaN; w = cosd(w); w = w.*sum(~isnan(w(:)))./nansum(w(:)); % 4°
    %w = 89.5:-1:-69.5; w = repmat(w,360,1)'; w(isnan(m2)) = NaN; w = cosd(w); w = w.*sum(~isnan(w(:)))./nansum(w(:)); % 1°

t2u(t2u==0) = NaN; m2u(m2u==0) = NaN; % remove points with miscalculated uncertainties
W = w(:); W = W(~isnan(t2u(:))); % inputs for York regression
uX = m2u(:); uX = uX(~isnan(t2u(:))); uY = t2u(:); uY = uY(~isnan(t2u(:)));
X = m2(:); X = X(~isnan(t2u(:))); Y = t2(:); Y = Y(~isnan(t2u(:))); 

[INT, SLP, INTu, SLPu] = york_fit(X',Y',uX'./sqrt(W'),uY'./sqrt(W'));

%% STEP 6 find effect of temperature on green-greener-blue-bluer

w = 89:-2:-69; w = repmat(w,180,1)'; w(isnan(m2)) = NaN; w = cosd(w); w = w.*sum(~isnan(w(:)))./nansum(w(:)); % latitude weighting
tTu(tTu==0) = NaN; t2u(isnan(tTu)) = NaN; % remove points with miscalculated uncertainties
W = w(:); W = W(~isnan(t2u(:))); % inputs for York regression
A = t2-(SLP.*m2+INT); % anomalies in Chl trends
Y = A(:); Y = Y(~isnan(t2u(:))); X = tT(:); X = X(~isnan(t2u(:))); % inputs for York regression
uY = t2u(:); uY = uY(~isnan(t2u(:))); uX = tTu(:); uX = uX(~isnan(t2u(:)));

[INT2, SLP2, INT2u, SLP2u] = york_fit(X',Y',uX'./sqrt(W'),uY'./sqrt(W'));

%% STEP 7 Gini

clearvars -EXCEPT C2;
C2 = squeeze(nanmean(C2,2)); % bin into latitude bands
C2 = C2(16:45,:); % keep only 0-60N
C = nanmean(C2,2); % get the time-mean version
[C,i] = sort(C,'descend'); % sort it
N = sum((1:30)*C)/sum(C); % calculate the number of bands positively related to Gini index
twobox = max(i(1:N))-N; % Gini is a ~two-box calculatin if this is <0
for i = 1:20; % now find for each year
    [c,J(i,:)] = sort(C2(:,i),'descend');
    n(i) = sum((1:30)*c)/sum(c);
    tb(i) = max(J(i,1:n(i)))-n(i);
end

figure;
subplot(121)
plot(1:30,cumsum(flipud(C))./sum(C),'k','linewidth',4)
box on;
axis square;
hold on;
set(gca,'ticklabelinterpreter','latex','fontsize',18)
ylabel('Cumulative Share of Mean Chl (0-60$^\circ$N)','interpreter','latex')
xlabel('Rank of Latitudinal Band (2$^\circ$)','interpreter','latex')
axis([.5 30.5 0 1.01])
p3 = plot(linspace(21.5,21.5),linspace(-1,2),'-.','color',[.95 .675 .675],'linewidth',4)
d = flipud(C);
p2 = bar(1:21,cumsum(d(1:21))./sum(d),'facecolor',[.05 .15 .5],'edgecolor',[.05 .15 .5])
d = cumsum(d);
p1 = bar(22:30,d(22:30)./d(end),'facecolor',[.05 .5 .15],'edgecolor',[.05 .5 .15])
lgnd = legend([p1 p2 p3],'Bands 42-60$^\circ$N','Bands 0-40$^\circ$N','$J = \sum(C^\downarrow (i)_{i=1}^N)/\sum C$')
set(gca,'ytick',0:.2:1,'xtick',5:5:30)
set(lgnd,'interpreter','latex','fontsize',18,'location','northwest');
clear lgnd;

subplot(122)
t = 2003:2022;
tt = repmat(t,30,1)';
nn = round(31-n)-.5;
t9 = tt(:,1:9); j9 = J(:,1:9);
scatter(t9(:),31-j9(:),200,[.05 .5 .15],'s','filled')
hold on;
t21 = tt(:,10:end); j21 = J(:,10:end);
scatter(t21(:),31-j21(:),200,[.05 .15 .5],'s','filled')
scatter(t,nn,100,[.95 .675 .675],"_",'linewidth',4)
axis([2002.5 2022.5 0 31])
box on;
axis square;
set(gca,'ticklabelinterpreter','latex','fontsize',18)
ylabel('Rank of Latitudinal Band (2$^\circ$)','interpreter','latex')
xlabel('Year','interpreter','latex')
set(gca,'ytick',5:5:30,'xtick',2005:5:2020)
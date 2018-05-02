%----------------------------------------------------------------------------
%STOchastic Rainfall Model (STORM). A rainstorm generator in this case based on empirical data (Walnut Gulch, AZ)

%Code name: WG_storms_v3_01  %this version allows for simulations that are much longer than the input series in order to compare the distributions of storm characteristics to the original.
%It also includes the 'mean gauge approach' to determine when a simulation year stops. This involves summing storm totals at each gauge for each year until the mean for all gauges exceeds the
%selected annual total precip value. It also allows for fuzzy selection of intensity at the storm center based on a fixed value of duration and
%incorporates orographic effects, wherein there are separate intensity-duration curves derived for three intervals of basin elevation (1200-1350m, 1351-1500m, 1501-1650m)
%Current version also includes interarrival times between storms, allowing for output to drive other model frameworks (rainfall-ruonff, water balance,LEMs)
%This version will also include output at each grid location, rather than only at gauge locations.
%Author: Michael Singer 2017
%Date created: 2015-6
%----------------------------------------------------------------------------

clear;
t0 = clock;
t1 = [datestr(floor(now)) '_' datestr(rem(now,1))];
t2=regexprep(t1,'[^\w'']',''); %for naming output directories and files by current date/time
mkdir('C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2)

%Initialize variables for annual rainfall total (mm/h) storm center location (RG1-RG85), etc.

%add variable for number of simulations of simyears
simyears = 1; %number of years to simulate
numgauges = length(Xin); %number of rain gauges in the basin. NOTE: In this version this produces output on a grid, rather than at real gauge locations.
numcurves = 11; %number of intensity-duration curves (see below for curve equations)

storm_scaling = 1.00; %No storm scaling, as problem appears to be fixed with smaller grid spacing.
%storm_scaling = 1.15; %This scales the storm center intensity upward, so the values at each gauge are realistic once the gradient is applied.

load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Ptot_pdf % This is the pdf fitted to all available station precip data (normal dist). It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Duration_pdf % This is the pdf fitted to all available station duration data (GEV dist). It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Area_pdf % This is the pdf fitted to all available station area data (EV dist). It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Int_arr_pdf % This is the pdf fitted to all available station area data (GEV dist). It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Recess_pdf % This is the pdf of storm gradient recession coefficiencts from Morin et al, 2005 (normal dist). It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Easting % This is the Longitudinal data for each gauge.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Northing % This is the Latitudinal data for each gauge. It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\gauges % This is the list of gauge numbers. It will be sampled below.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\X % This is the Longitudinal data for each grid point. It will be sampled below to determine storm center location.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Y % This is the Latitudinal data for each grid point. It will be sampled below to determine storm center location.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Storm_depth_data % This is the storm depth data for use in model evaluation.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Intensity_data % This is the intensity data for use in model evaluation.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Duration_data % This is the duration data for use in model evaluation.
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Gauge_GR1 % This is the lowest elevation gauge grouping used for orography (1200-1350m).
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Gauge_GR2 % This is the middle elevation gauge grouping used for orography (1351-1500m).
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Gauge_GR3 % This is the highest elevation gauge grouping used for orography (1501-1650m).
load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\fuzz % This is a vector of fuzzy tolerace values for intensity selection.
a = shaperead('boundary/boundary.shp'); %This is shapefile of the watershed boundary
Yy = extractfield(a,'Y');
Xx = extractfield(a,'X');
[X1,Y1]=meshgrid(linspace(min(X),max(X),round((max(X)-min(X))/1000)),linspace(min(Y),max(Y),round((max(Y)-min(Y))/1000))); %creates a mesh with 1000m spacings 
isin=inpolygon(X1(:),Y1(:),Xx,Yy); 
Yin = Y1(isin);
Xin = X1(isin); %creates a list of points, Xin, Yin on the grid that are within the watershed boundary

%figure

%initialize matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Storm matrix = [Storm# StormArea StormDuration Int_DurCurve# StormIntensity #GaugesHit RecessionVal StormTotal Longitude Latitude Year]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Storm_matrix = zeros(500*simyears,11); %based on presumed maximum number of storms per year = 300
Ptot_ann_global = zeros(1,simyears);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gauge matrix = [Year Storm# StormIntensity StormDuration StormTotal Ann_Cum_PTotal]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:numgauges %initialize matrices of storm output for each gauge
    eval(['Gauge_matrix_',num2str(ii),' = zeros(500*simyears,6);']) %based on presumed maximum number of storms per year = 300
end

%Intensity_ave_local = zeros(simyears,numgauges); %these matrices will be analyzed for summary statistics and plotted in ArcGIS
% Duration_ave_local = zeros(simyears,numgauges);
% Ptot_sum_local = zeros(simyears,numgauges);

% Intensity_all = sparse(zeros(max_vals,1)); %these vectors will contain all data generated at all gauges for statistical comparison with input pdfs.
% Duration_all = sparse(zeros(max_vals,1));
% Ptot_all = sparse(zeros(max_vals,1));
%Area_all = sparse(zeros(max_vals,1));

Intensity_local_all = 0; %initialize all variables (concatenated matrices of generated output)
%Intensity_local_all = zeros(85*simyears,1); %initialize local_all variables (concatenated vector of generated output at each gauge location)
Storm_totals_all = 0;
%Ptot_ann_global = 0;
Duration_all = 0;
Duration_local_all = 0;
Gauges_hit_all = 0;

%int_dur_curve_num = (1:numcurves)'; % these were empirically determined based on data from WG (monsoon rainfall only).
gauge_locs = [Easting Northing]; %these are the UTM coords of the gauges
storm_count = 0;
master_storm_count = 0;
for year = 1:simyears
    Ptotal = 0;
    Ptot_ann_global(year) = random(Ptot_pdf,1,true); %samples from normal distribution and saves global value of Ptot (that must be equalled or exceeded) for each year
    clear Storm_total_local_year
    ann_cum_Ptot_gauge = zeros(1,numgauges);
    for storm = 1:5000
        clear cx cy r mask_name aa bb cc dd* North_hit North East East_hit rain_int_gauge intensity_val duration_val recess_val Ptotal* gdist center_val* int_dur_curve_num*
        rain_int_gauge = zeros(1,numgauges);
        storm_count = storm_count+1;
        master_storm_count = master_storm_count+1;
        center_val_X = datasample(X,1); %sample uniformly from storm center matrix from grid with 10 m spacings covering basin).
        center_val_Y = datasample(Y,1);
        North = center_val_Y;
        East = center_val_X;
          
        area_val = random(Area_pdf,1,true); %Samples from distribution of storm areas
        Storm_matrix(master_storm_count,1) = master_storm_count; %storm #
        Storm_matrix(master_storm_count,2) = area_val;
        cx = East; %value of coord should be set to storm center selected (same below)
        cy = North;
        Storm_matrix(master_storm_count,9) = cx;
        Storm_matrix(master_storm_count,10) = cy;
        Storm_matrix(master_storm_count,11) = year;
        r = sqrt(area_val/pi); %value here should be selected based on area above in meters to match the UTM values in North and East vectors.
        mask_name = (((Easting(:)-cx).^2 + (Northing(:)-cy).^2) <= r^2); %determine which gauges are hit by Euclidean geometry
        aa = find(mask_name == 1);
        gauges_hit = gauges(aa);
        num_gauges_hit = length(gauges_hit);
        %this routine below allows for orography in precip by first determining the closest gauge and then determining its orographic grouping
        gdist = (Easting(:)-cx).^2 + (Northing(:)-cy).^2;
        bb = min(gdist);
        cc = find(gdist == bb);
        closest_gauge = gauges(cc); %this will be compared against orographic gauge groupings to determine the appropriate set of intensity-duration curves
        %%%%%%
        Storm_matrix(master_storm_count,6) = num_gauges_hit;
        %gauges_hit_all(:,year) = [gauges_hit_all(:,year); gauges_hit];
        Gauges_hit_all = [Gauges_hit_all; gauges_hit]; %#ok<AGROW>
        North_hit = Northing(aa);
        East_hit = Easting(aa);
               
%                 if ~isempty(gauges_hit)
%                     viscircles([cx cy],r); %draw a circle with a particular radius around each storm center location
%                 end
        %this routine below determines to which orographic group the closest gauge to the storm center belongs to, and censors the number of curves accordingly
        %missing top curve in GR1, top and bottom curves for GR2, and bottom curve for GR3
        dd = ismember(Gauge_GR1,closest_gauge);
        dd1 = find(dd>0, 1);
        if ~isempty(dd1)
            %int_dur_curve_num = (2:numcurves)'; % these were empirically determined based on data from WG (monsoon rainfall only)-lowest orographic group.
            baa = 'a';
        end
        dd = ismember(Gauge_GR2,closest_gauge);
        dd1 = find(dd>0, 1);
        if ~isempty(dd1)
            %int_dur_curve_num = (2:numcurves-1)'; % these were empirically determined based on data from WG (monsoon rainfall only)-middle orographic group.
            baa = 'b';
        end
        dd = ismember(Gauge_GR3,closest_gauge);
        dd1 = find(dd>0);
        if ~isempty(dd1)
            %int_dur_curve_num = (1:numcurves-1)'; % these were empirically determined based on data from WG (monsoon rainfall only)-highest orographic group.
            baa = 'c';
        end
        int_dur_curve_num = (1:numcurves)';
        duration_val = random(Duration_pdf,1,true);
        duration_val = round(duration_val); %round to nearest minute for consistency with measured data
        %Duration_global(storm,year) = duration_val;
        Storm_matrix(master_storm_count,3) = duration_val;
        %Duration_all = [Duration_all; duration_val];
        
        %         int_dur_curve_numy = int_dur_curve_num;
        %         for poo = 1:100
        %             int_dur_curve_numy = [int_dur_curve_numy;int_dur_curve_num]; %this just repeats the sequence many times to allow datasample to function better.
        %         end
        %         int_dur_curve_num = int_dur_curve_numy;
        
        % original curve probs for 21-14-7%: [0.0718 0.0782 0.0845 0.0909 0.0909 0.0909 0.0909 0.0909 0.0973 0.1036 0.1100]
        % original curve probs for 24-16-8%: [0.0691 0.0764 0.0836 0.0909 0.0909 0.0909 0.0909 0.0909 0.0982 0.1055 0.1127]
        % original curve probs for 27-18-9%: [0.0664 0.0745 0.0827 0.0909 0.0909 0.0909 0.0909 0.0909 0.0991 0.1073 0.1155]
        % original curve probs for 30-20-10%: [0.0636 0.0727 0.0819 0.0909 0.0909 0.0909 0.0909 0.0909 0.1001 0.1090 0.1182]
        
        if baa == 'a'
            int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0318 0.0759 0.0851 0.0941 0.0941 0.0941 0.0941 0.0941 0.1033 0.1121 0.1213]); %add weights to reflect reasonable probabilities that favor lower curves.
        elseif baa == 'b'
            int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0478 0.0778 0.0869 0.0959 0.0959 0.0959 0.0959 0.0959 0.1051 0.1141 0.0888]); %add weights to reflect reasonable probabilities that favor lower curves.
        else
            int_dur_curve_val = datasample(int_dur_curve_num,1,'weights',[0.0696 0.0786 0.0878 0.0968 0.0968 0.0968 0.0968 0.0968 0.1060 0.1149 0.0591]); %add weights to reflect reasonable probabilities that favor lower curves.
            %Curve_num_global(storm,year) = int_dur_curve_val;
        end
        Storm_matrix(master_storm_count,4) = int_dur_curve_val;
        if int_dur_curve_val == 1
            intensity_val = 642.2*exp(-0.508*duration_val)+93.1*exp(-0.008*duration_val)+ 4.5;  %these curves are based on empirical data from WG.
        end
        if int_dur_curve_val == 2
            intensity_val = 578.0*exp(-0.508*duration_val)+83.8*exp(-0.008*duration_val)+ 4;
        end
        if int_dur_curve_val == 3
            intensity_val = 513.8*exp(-0.508*duration_val)+74.5*exp(-0.008*duration_val)+ 3.5;
        end
        if int_dur_curve_val == 4
            intensity_val = 449.5*exp(-0.508*duration_val)+65.2*exp(-0.008*duration_val)+ 3;
        end
        if int_dur_curve_val == 5
            intensity_val = 385.3*exp(-0.508*duration_val)+55.9*exp(-0.008*duration_val)+ 2.5;
        end
        if int_dur_curve_val == 6
            intensity_val = 321.1*exp(-0.508*duration_val)+46.6*exp(-0.008*duration_val)+ 2;
        end
        if int_dur_curve_val == 7
            intensity_val = 256.9*exp(-0.508*duration_val)+37.2*exp(-0.008*duration_val)+ 1.5;
        end
        if int_dur_curve_val == 8
            intensity_val = 192.7*exp(-0.508*duration_val)+27.9*exp(-0.008*duration_val)+ 1;
        end
        if int_dur_curve_val == 9
            intensity_val = 128.4*exp(-0.508*duration_val)+18.6*exp(-0.008*duration_val)+ 0.5;
        end
        if int_dur_curve_val == 10
            intensity_val = 64.1*exp(-0.508*duration_val)+9.3*exp(-0.008*duration_val)+ 0.25;
        end
        if int_dur_curve_val == 11
            intensity_val = 21*exp(-0.508*duration_val)+0.9*exp(-0.008*duration_val)+ 0.05;
        end
        
        fuzz_int_val = datasample(fuzz,1);
        intensity_val2 = intensity_val+fuzz_int_val; %allowing for +/-5 mm/hr fuzzy tolerance around selected intensity
        if intensity_val2 < 1 %cannot have zero or negative intensity
            intensity_val = 1;
        else
            intensity_val = intensity_val2;
        end
        intensity_val = intensity_val*storm_scaling; %This scales the storm center intensity upward, so the values at each gauge are realistic once the gradient is applied.
        Storm_matrix(master_storm_count,5) = intensity_val;
        %Intensity_global(storm,year) = intensity_val;
        %Intensity_all = [Intensity_all; intensity_val]; %selected storm center intensities
        
        %area to determine which gauges are hit
        
        
        recess_val = random(Recess_pdf,1,true); %this pdf of recession coefficients determines how intensity declines with distance from storm center (see below)
        Storm_matrix(master_storm_count,7) = recess_val;
        for j = 1:numgauges %determine cartesian distances to all hit gauges and associated intensity values at each gauge hit by the storm
            if ismember(j,gauges_hit)
                gauge_dist = sqrt(((Easting(j)-cx).^2 + (Northing(j)-cy).^2));
                gauge_dist_km = gauge_dist./1000; %change to km
                %rain_int_gauge(j) = intensity_val(storm)*exp(-2*(0.39^2)*gauge_dist_km.^2); %Rodriguez-Iturbe et al., 1986; Morin et al., 2005
                rain_int_gauge(j) = intensity_val*exp(-2*(recess_val^2)*gauge_dist_km.^2); %Rodriguez-Iturbe et al., 1986; Morin et al., 2005 but sampled from a distribution
            elseif ~ismember(j,gauges_hit)
                %rain_int_gauge(j) = NaN; %give NaN to gauges not hit
                rain_int_gauge(j) = 0; %give zero intensity to gauges not hit
            end
            ann_cum_Ptot_gauge(j) = ann_cum_Ptot_gauge(j)+rain_int_gauge(j)*duration_val/60;
        end
        
        for jj = 1:numgauges
            eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,1) = year;']) %year
            eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,2) = master_storm_count;']) %storm #
            eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,3) = rain_int_gauge(jj);'])
            eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,4) = duration_val;'])
            eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,5) = rain_int_gauge(jj)*duration_val/60;']) %storm total
            eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,6) = ann_cum_Ptot_gauge(jj);']) %ann cum total (Ptot)
        end
        %FIX THIS PART
        Intensity_local_all = [Intensity_local_all,rain_int_gauge]; %#ok<AGROW> %collect into vector of all simulated intensities (at all gauges)
        dur_step(1:85) = duration_val;
        Duration_local_all = [Duration_local_all,dur_step]; %#ok<AGROW>
        Storm_total_local_year(storm,1:numgauges) = rain_int_gauge.*duration_val/60; %#ok<SAGROW> %collect storm total data for all gauges into rows by storm
        Storm_totals_all = [Storm_totals_all,rain_int_gauge.*duration_val/60]; %#ok<AGROW>
        Storm_matrix(master_storm_count,8) = intensity_val*duration_val/60;
        for k = 1:numgauges
            Ptotal(k) = nansum(Storm_total_local_year(:,k)); %#ok<SAGROW> %sum the annual storm total at each gauge
        end
        %Ptotal_test = find((nanmean(Ptotal) + nanstd(Ptotal)/sqrt(85)) > Ptot_ann_global(year));    %once the mean of all gauges exceeds the selected annual storm total, a new simulation year begins
        Ptotal_test = find(nanmedian(Ptotal) > Ptot_ann_global(year));    %once the median of all gauges exceeds the selected annual storm total, a new simulation year beginsPtotal_test = find(Ptotal > Ptot_ann_global(year));    %once the first gauge exceeds the selected annual storm total, a new simulation year begins
        %Ptotal = Ptotal + duration_val(storm)/60*intensity_val(storm);
        if ~isempty(Ptotal_test)
            %eval(['Ptotal_local_',num2str(year),'(1:numgauges) = Ptotal;'])
            %Ptotal_local(year,1:numgauges) = Ptotal;
            for l = 1:numgauges
                Ptotal_local(year,l) = Ptotal(l); %#ok<SAGROW>
            end
            break %start a new simulation year
        end
        
        %local and concatenated vector output
        %concatenate data into single vectors
        %         a = find(Intensity_ave_local == 0);
        %         Intensity_ave_local(a) = NaN;
        
        
        %         a = find(intensity_val == 0); %FIX
        %         intensity_val(a) = NaN;
        %         a = find(duration_val == 0);
        %         duration_val(a) = NaN;
        %         a = find(area_val == 0);
        %         area_val(a) = NaN;
        %eval(['ann_gauges_',num2str(year),'(:,storm) = gauges_hit(:,storm);'])
        %         eval(['ann_intensity_',num2str(year),'(storm) = intensity_val;'])
        %         eval(['ann_duration_',num2str(year),'(storm) = duration_val./60;']) %convert duration values to hours.
        %         eval(['ann_area_',num2str(year),'(storm) = area_val;'])
        eval(['Storm_total_local_year_',num2str(year),'(storm,1:numgauges) = Storm_total_local_year(storm,1:numgauges);']) %collect all local annual storm totals for each gauge.
        int_arr_val = random(Int_arr_pdf,1,true); %Samples from distribution of interarrival times (hours). This can be used to develop STORM output for use in rainfall-runoff models or any water balance application.
    end %storm loop
end %year loop
%hold on
%plot(Easting,Northing,'o')
%grid on
AA = find(Storm_matrix(:,9)>0); %gets rid of trailing zeros from initialized matrix
Storm_matrix = Storm_matrix(AA,:);
AB = find(Gauge_matrix_1(:,2)>0);
for CC = 1:numgauges
    eval(['Gauge_matrix_',num2str(CC),' = Gauge_matrix_',num2str(CC),'(AB,:);'])
end
% for ii = 1:numgauges
%     Storm_total_mean_local(ii) = nanmean(Storm_total_local(:,ii));
%     Storm_total_med_local(ii) = nanmedian(Storm_total_local(:,ii));
%     Storm_total_max_local(ii) = nanmax(Storm_total_local(:,ii));
%     Storm_total_in_local(ii) = nanmin(Storm_total_local(:,ii));
% end

%concatenate data into single vectors
% for jj = 1:simmyears %years of data
% Gauges_hit_all = horzcat(gauges_hit_all(1,:),gauges_hit_all(jj,:));
% end

%global output
% a = find(Area_global == 0);
% Area_global(a) = NaN;
% a = find(Duration_global == 0);
% Duration_global(a) = NaN;
% a = find(Curve_num_global == 0);
% Curve_num_global(a) = NaN;
%a = find(Storm_center_global == 0);
%Storm_center_global(a) = NaN;
a = find(Gauges_hit_all == 0);
Gauges_hit_all(a) = NaN;

GZ = find(Intensity_local_all>0);
Intensity_all = Intensity_local_all(GZ);
Duration_all = Duration_local_all(GZ);
Storm_totals_all = Storm_totals_all(GZ);
Ptot_ann_global = Ptot_ann_global(2:length(Ptot_ann_global));
Gauges_hit_all = Gauges_hit_all(2:length(Gauges_hit_all));

eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Ptot_ann_',num2str(simyears),'y_global_',t2,' Ptot_ann_global'])
eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Storm_matrix_',num2str(simyears),'y_',t2,' Storm_matrix'])
eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Gauges_hit_',num2str(simyears),'y_all_',t2,' Gauges_hit_all'])
eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Storm_totals_',num2str(simyears),'y_all_',t2,' Storm_totals_all'])
eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Intensity_',num2str(simyears),'y_selected_',t2,' Intensity_all'])
eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Duration_',num2str(simyears),'y_selected_',t2,' Duration_all'])

for kk = 1:numgauges
    eval(['save C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2,'\Gauge_matrix_',num2str(kk),'_y_',t2,' Gauge_matrix_',num2str(kk);])
end

boo = etime(clock,t0);
runtime_seconds = boo %#ok<NOPTS>
runtime_minutes = boo/60 %#ok<NOPTS>










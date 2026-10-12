function summary = mseed_cross_corr_2_branch(alignedH5, targetCSV, outputDir, repeatCycles, slideCycles, taperName, options)
% Studio 2.5.3: the uploaded MATLAB Goertzel/MAD/sign-branch framework.
% options defaults: stabilization_s=600, frequency_method='matlab_bin',
% amplitude_rejection=true, split_branches=true, screening_mode='curve'.
% screening_mode='matlab' uses each branch's preliminary phase +/-90 degrees.
% Corrected aligned H5 preserves gaps; raw discontinuous MiniSEED is not joined.
% The unchanged uploaded source is reference/mseed_cross_corr_2.m.
if nargin<1 || isempty(alignedH5)
    [f,p]=uigetfile('*.h5','Select aligned_waveforms.h5');
    if isequal(f,0), summary=table(); return; end
    alignedH5=fullfile(p,f);
end
if nargin<2 || isempty(targetCSV)
    [f,p]=uigetfile('*.csv','Select selected_phase_targets.csv');
    if isequal(f,0), summary=table(); return; end
    targetCSV=fullfile(p,f);
end
if nargin<3 || isempty(outputDir), outputDir=fullfile(fileparts(alignedH5),'matlab_phase_results'); end
if nargin<4, repeatCycles=20; end
if nargin<5, slideCycles=1; end
if nargin<6, taperName='none'; end
if nargin<7,options=struct();end
settings=struct('stabilization_s',600,'frequency_method','matlab_bin','amplitude_rejection',true,'split_branches',true,'screening_mode','curve');
fields=fieldnames(settings);for i=1:numel(fields),if ~isfield(options,fields{i}),options.(fields{i})=settings.(fields{i});end,end
assert(isfinite(options.stabilization_s) && options.stabilization_s>=0,'Invalid stabilization time');
assert(any(strcmp(options.frequency_method,{'matlab_bin','exact'})),'Invalid frequency method');
assert(any(strcmp(options.screening_mode,{'curve','matlab'})),'Invalid screening mode');
assert(isfinite(repeatCycles) && repeatCycles>0 && isfinite(slideCycles) && slideCycles>0,'Invalid cycles');
assert(any(strcmp(taperName,{'none','hann'})),'taperName must be none or hann');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
T=readtable(targetCSV,'TextType','string');
required={'curve_id','segment_id','pair','branch','frequency_hz','predicted_phase_deg','phase_offset_low_deg','phase_offset_high_deg'};
assert(all(ismember(required,T.Properties.VariableNames)),'Missing required target columns');
assert(all(isfinite(T.frequency_hz) & T.frequency_hz>0),'Invalid frequencies');
assert(all(T.cross_convention=="Gtarget*conj(Greference)"),'Cross-spectrum convention mismatch');
fs=double(h5readatt(alignedH5,'/','fs_hz'));
startUTC=string(h5readatt(alignedH5,'/','start_utc'));
assert(logical(h5readatt(alignedH5,'/','response_corrected')),'Response corrected data required');
assert(all(T.frequency_hz<fs/2),'Requested frequency exceeds Nyquist');
allStations=unique(split(join(unique(T.pair),"--"),"--"));
D=containers.Map('KeyType','char','ValueType','any');
N=[];B=containers.Map('KeyType','char','ValueType','any');
for i=1:numel(allStations)
    name=char(allStations(i));
    x=double(h5read(alignedH5,['/stations/' name '/data']));x=x(:);
    valid=logical(h5read(alignedH5,['/stations/' name '/valid']));valid=valid(:);
    assert(numel(x)==numel(valid),'Data/mask length mismatch');
    x(~valid)=NaN;
    if isempty(N),N=numel(x);else,assert(N==numel(x),'Station lengths not UTC aligned');end
    D(name)=x;bounds=[0,numel(x)];
    try,bounds(1)=double(h5readatt(alignedH5,['/stations/' name],'record_start_sample'));catch,end
    try,bounds(2)=double(h5readatt(alignedH5,['/stations/' name],'record_end_sample_exclusive'));catch,end
    B(name)=bounds;
end
rows={};gateNo=0;
frequencies=unique(T.frequency_hz);
for fi=1:numel(frequencies)
    f=frequencies(fi);n=round(repeatCycles*fs/f);step=max(1,round(slideCycles*fs/f));
    assert(n>=2,'Invalid window length');
    actualF=f;if strcmp(options.frequency_method,'matlab_bin'),actualF=round(f*n/fs)*fs/n;end
    assert(actualF>0 && actualF<fs/2,'Invalid Goertzel frequency');
    w=ones(n,1);
    if strcmp(taperName,'hann'),w=.5-.5*cos(2*pi*(0:n-1)'/(n-1));end
    Q=T(abs(T.frequency_hz-f)<1e-10,:);
    cache=containers.Map('KeyType','char','ValueType','any');
    segments=unique(Q.segment_id,'stable');
    for si=1:numel(segments)
        R=Q(Q.segment_id==segments(si),:);r=R(1,:);pair=split(r.pair,"--");
        assert(numel(pair)==2 && all(R.pair==r.pair),'Inconsistent pair');
        key=char(r.pair);
        if ~isKey(cache,key)
            refBounds=B(char(pair(1)));targetBounds=B(char(pair(2)));
            first=max(1,max(refBounds(1),targetBounds(1))+round(options.stabilization_s*fs)+1);
            last=min([N,refBounds(2),targetBounds(2)]);starts=(first:step:last-n+1)';
            Gref=exactGoertzelWindows(D(char(pair(1))),starts,n,actualF,fs,w);
            Gtarget=exactGoertzelWindows(D(char(pair(2))),starts,n,actualF,fs,w);
            cache(key)=struct('starts',starts,'cross',Gtarget.*conj(Gref));
        end
        entry=cache(key);starts=entry.starts;C=entry.cross;A=abs(C);phase=angle(C);
        center=deg2rad(r.predicted_phase_deg);lo=deg2rad(r.phase_offset_low_deg);hi=deg2rad(r.phase_offset_high_deg);
        assert(isfinite(center) && isfinite(lo) && isfinite(hi) && hi>=lo,'Invalid phase bounds');
        assert(all(abs(angle(exp(1i*deg2rad(R.predicted_phase_deg-r.predicted_phase_deg))))<1e-8),'Aliases have inconsistent centers');
        assert(all(abs(R.phase_offset_low_deg-r.phase_offset_low_deg)<1e-8) && all(abs(R.phase_offset_high_deg-r.phase_offset_high_deg)<1e-8),'Aliases have inconsistent ranges');
        valid=isfinite(real(C)) & isfinite(imag(C)) & A>0;
        if hi-lo>=2*pi-1e-12, inRange=true(size(C));
        else,inRange=abs(angle(exp(1i*(phase-center-(lo+hi)/2))))<=(hi-lo)/2+1e-12;end
        [amplitudeLarge,amplitudeJump]=amplitudeFlags(A);amplitudeValid=valid;
        if options.amplitude_rejection,amplitudeValid=valid & ~(amplitudeLarge|amplitudeJump);end
        positive=amplitudeValid & phase>=0;negative=amplitudeValid & phase<0;
        if strcmp(options.screening_mode,'matlab')
            preliminaryPositive=sum(C(positive));preliminaryNegative=sum(C(negative));cp=NaN;cn=NaN;
            if abs(preliminaryPositive)>0,cp=angle(preliminaryPositive);end
            if abs(preliminaryNegative)>0,cn=angle(preliminaryNegative);end
            acceptedPositive=positive & abs(angle(C.*exp(-1i*cp)))<=pi/2;
            acceptedNegative=negative & abs(angle(C.*exp(-1i*cn)))<=pi/2;
        else
            acceptedPositive=positive & inRange;acceptedNegative=negative & inRange;
        end
        if options.split_branches || strcmp(options.screening_mode,'matlab')
            if center>=0,accepted=acceptedPositive;else,accepted=acceptedNegative;end
        else,accepted=acceptedPositive | acceptedNegative;end
        sumPositive=sum(C(acceptedPositive));sumNegative=sum(C(acceptedNegative));
        z=C(accepted);rawSum=sum(z);unitSum=sum(z./abs(z));count=sum(accepted);
        rawPhase=NaN;unitPhase=NaN;unitR=NaN;meanAmplitude=NaN;
        if count>0
            if abs(rawSum)>0,rawPhase=rad2deg(angle(rawSum));end
            if abs(unitSum)>0,unitPhase=rad2deg(angle(unitSum));end
            unitR=abs(unitSum)/count;meanAmplitude=abs(rawSum)/count;
        end
        gateNo=gateNo+1;gateFile=sprintf('gate_%05d.mat',gateNo);
        W=struct('start_sample_zero_based',starts-1,'start_offset_s',(starts-1)/fs,'start_utc',startUTC,'frequency_hz',f, ...
            'actual_frequency_hz',actualF,'options',options,'accepted_positive',acceptedPositive,'accepted_negative',acceptedNegative, ...
            'total_cross_positive',sumPositive,'total_cross_negative',sumNegative,'window_samples',n,'slide_samples',step,'fs_hz',fs,'taper',taperName,'pair',r.pair,'segment_id',r.segment_id, ...
            'curve_ids',R.curve_id,'cross',C,'observed_phase_deg',rad2deg(phase),'valid',valid,'accepted',accepted, ...
            'amplitude_large',amplitudeLarge,'amplitude_jump',amplitudeJump,'predicted_phase_deg',r.predicted_phase_deg, ...
            'phase_offset_low_deg',r.phase_offset_low_deg,'phase_offset_high_deg',r.phase_offset_high_deg, ...
            'raw_cumulative',cumsum(z),'unit_cumulative',cumsum(z./abs(z)));
        save(fullfile(outputDir,gateFile),'W','-v7');
        for ci=1:height(R)
            rows(end+1,:)={R.curve_id(ci),r.segment_id,r.pair,r.branch,f,r.predicted_phase_deg,r.phase_offset_low_deg,r.phase_offset_high_deg, ...
                numel(C),sum(valid),count,count/max(1,sum(valid)),sum((amplitudeLarge|amplitudeJump)&accepted), ...
                real(rawSum),imag(rawSum),rawPhase,meanAmplitude,unitPhase,unitR,string(gateFile),actualF,sum(acceptedPositive),sum(acceptedNegative), ...
                real(sumPositive),imag(sumPositive),phaseOrNaN(sumPositive),real(sumNegative),imag(sumNegative),phaseOrNaN(sumNegative)}; %#ok<AGROW>
        end
        if gateNo==1
            figure('Color','w');subplot(2,1,1);residual=rad2deg(angle(exp(1i*(phase-center))));
            histogram(residual(valid),-180:6:180);hold on;histogram(residual(accepted),-180:6:180);xline(rad2deg(lo));xline(rad2deg(hi));
            xlabel('Circular phase residual (deg)');ylabel('Windows');legend('Before','Accepted');title(sprintf('%s %.2f Hz',r.segment_id,f));
            subplot(2,1,2);plot(rad2deg(angle(cumsum(z))));hold on;plot(rad2deg(angle(cumsum(z./abs(z)))));xlabel('Accepted windows');ylabel('Stack phase (deg)');legend('Raw','Unit phase');
            saveas(gcf,fullfile(outputDir,'first_gate.png'));
        end
    end
end
columns={'curve_id','segment_id','pair','branch','frequency_hz','predicted_phase_deg','phase_offset_low_deg','phase_offset_high_deg', ...
    'total_windows','valid_windows','accepted_windows','accepted_fraction','amplitude_marked_but_accepted','raw_real','raw_imag', ...
    'raw_phase_deg','raw_mean_amplitude','unit_phase_deg','unit_resultant_length','windows_file','actual_frequency_hz', ...
    'accepted_positive_windows','accepted_negative_windows','positive_raw_real','positive_raw_imag','positive_raw_phase_deg', ...
    'negative_raw_real','negative_raw_imag','negative_raw_phase_deg'};
summary=cell2table(rows,'VariableNames',columns);writetable(summary,fullfile(outputDir,'phase_stack_summary.csv'));
save(fullfile(outputDir,'phase_stack_summary.mat'),'summary','repeatCycles','slideCycles','taperName','startUTC','fs');
assignin('base','curveGuidedSummary',summary);
fprintf('Finished %d gates, %d curve/frequency rows. Results: %s\n',gateNo,height(summary),outputDir);
end

function G=exactGoertzelWindows(x,starts,N,f,fs,taper)
% Standard sum(x[n]*exp(-i*w*n)) at the exact frequency; no nearest-bin rounding.
% Vectorize recurrence across windows. Correct common recurrence phase rotation.
w=2*pi*f/fs;coeff=2*cos(w);mu=zeros(size(starts));valid=true(size(starts));
for j=0:N-1
    v=x(starts+j);valid=valid & isfinite(v);v(~isfinite(v))=0;mu=mu+v;
end
mu=mu/N;s1=zeros(size(starts));s2=s1;
for j=0:N-1
    v=x(starts+j);v(~isfinite(v))=0;
    s0=(v-mu)*taper(j+1)+coeff*s1-s2;s2=s1;s1=s0;
end
G=(s1-exp(-1i*w)*s2)*exp(-1i*w*(N-1));G(~valid)=complex(NaN,NaN);
end

function [large,jump]=amplitudeFlags(A)
vals=A(isfinite(A));large=false(size(A));jump=large;
if isempty(vals),return;end
m=median(vals);s=1.4826*median(abs(vals-m));if s>0,limit=m+5*s;else,limit=5*m;end
large=A>limit;d=[0;abs(diff(A))];v=d(2:end);v=v(isfinite(v));if isempty(v),v=0;end;m=median(v);s=1.4826*median(abs(v-m));if s>0,limit=m+5*s;else,limit=5*m;end
jump=d>limit;
end

function phase=phaseOrNaN(value)
phase=NaN;if abs(value)>0,phase=rad2deg(angle(value));end
end

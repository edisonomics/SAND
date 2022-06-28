close all;
clear all;
%% Set your toolbox paths; functions imported from these directories:
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox='/Users/yuewu/Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/';
% NMR decompositon program found @ https://github.com/edisonomics/SAND
localPaths.nmrdecomp_path='/Users/yuewu/Documents/GitHub/SAND/';
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.nmrdecomp_path));
pause(1),clc
% the path should be changed accordingly in hte users' computer
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/'
projdir=[paredir 'result_reproduce/urine_spikin/'];
datadir=[projdir 'data/'];
libdir=[projdir 'archive/data_exap/1/test_trans.fid'];% a template fid file containing useful header information. https://www.dropbox.com/s/1i0dixw4vasctwu/test_trans.fid?dl=0
preresdirpath=[projdir 'res/res/'];
cd([projdir]);
load([projdir 'res/saved_preprocessing.mat']);
specppm=ppm_r;
%
sampleseq=1:nsample;
ppmrange_dss=[-0.1 0.1];
deltapm_threshold=0.01;%distance threshold for peak matching
% vis check of the deconv
for i=sampleseq
  foldpath=[preresdirpath num2str(i) '/'];
  dircont=dir([foldpath '*ft*.fig'])
  figpath=[foldpath dircont.name];
  uiopen(figpath,1);
end
close all;
obj_refine=[];
time_cost=[];
for runid=sampleseq
  load([preresdirpath num2str(runid) '/runid' num2str(runid) '_refine_res.mat']);
  obj_refine=[obj_refine obj_scaled];
  time_cost=[time_cost timecost];
end
tab_eval=table(sampleseq',obj_refine',time_cost','VariableNames',{'id','obj_refine','time_cost'});
save('performance_eval.mat','tab_eval');
% mean performance
mean(log10(tab_eval{:,'obj_refine'}))
% mean time
mean(tab_eval{:,'time_cost'})
% load raw FT of all samples
ft_mat=[];
for runid=sampleseq
  runid
  sampi_str=num2str(runid);
  load([preresdirpath sampi_str '/runid' sampi_str '_trainingdata.mat']);
  ft_mat=[ft_mat;ft_ori_tab{:,2}'];
end
fig=figure();
plotr(ppm,ft_mat);
savefig(fig,['ft_all_sample_stack.fig']);
close all;
%
stackSpectra(ft_mat,ppm,0,50,'all stack')
fig=gcf;
saveas(fig,['stack_all_sample.fig']);
close all;
% load decomposation estimation of different spectra
namelist={};
est_tab=[];%PPM, lambda, A, simulation_ind
ppmstr={};
ftstr={};
for sampi=sampleseq
  samplestr=num2str(sampi);
  load([preresdirpath samplestr '/runid' samplestr '_env_final.mat']);
  runtab=array2table(tabsumm_refine,'VariableNames',{'PPM','lambda','A','phase'});%f, lambda, A, phi
  nfeature=size(runtab,1);
  runtab{:,'PPM'}=runtab{:,'PPM'}/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  simu=repmat(sampi,[nfeature,1]);
  runtab=[runtab table(simu)];
  namelist=[namelist; {repmat({'unknown'},[1,nfeature])}];
  est_tab=[est_tab; runtab];
  ppmstr{sampi}=ppm;
  ftstr{sampi}=ft_ori_tab;
end
est_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation'};
%intensity and integral based estimation
est_other_tab=[];
for sampi=sampleseq
  spec_here=specmat(sampi,:);
  runtab=[];
  samplestr=num2str(sampi);
  load([preresdirpath samplestr '/runid' samplestr '_env_final.mat']);
  for bini=1:size(binrange,1)
    ppmrange=binrange(bini,:);
    [indrang]=sort(matchPPMs(ppmrange,specppm));
    indseq=indrang(1):indrang(2);
    spec_reg_shift=spec_here(indseq);
    baseval=min(spec_reg_shift);
    % baseval=0;
    [est_inten maxind]=max(spec_reg_shift-baseval);
    est_auc=trapz(spec_reg_shift-baseval);
    est_temp_tab=table(specppm(indrang(1)+maxind-1),nan,est_inten,est_auc,0,sampi);
    runtab=[runtab; est_temp_tab];
  end
  est_other_tab=[est_other_tab; runtab];
end
est_other_tab.Properties.VariableNames={'PPM','lambda','intensity','integral','phase','simulation'};
%
quan_str=struct();
quan_str.deconv=est_tab;
quan_str.intensity=est_other_tab(:,{'PPM','lambda','intensity','phase','simulation'});
quan_str.integral=est_other_tab(:,{'PPM','lambda','integral','phase','simulation'});
%
save('data_temp.mat')
% match peaks from different samples
% example visualization regions
exampreg=[0.965 1.047; 1.457 1.495; 1.73 1.774; 4.61 4.67; 5.2 5.3; 6.887 6.917; 7.339 7.353; 7.851 7.869];
% peak match for spike in peaks [ref bin, mixture bin]
ppmref_match=[0.97 0.98 0.97 0.98; 0.983 0.99 0.983 0.99; 1.02 1.03 1.02 1.03; 1.035 1.042 1.035 1.042; 1.468 1.473 1.472 1.477; 1.48 1.485 1.484 1.489; 1.733 1.746 1.733 1.746; 1.746 1.761 1.746 1.761; 1.761 1.773 1.761 1.773; 1.895 1.902 1.897 1.904; 3.373 3.38 3.383 3.388; 3.389 3.396 3.399 3.405; 3.593 3.599 3.605 3.609; 3.601 3.607 3.612 3.616; 4.624 4.637 4.633 4.645; 4.637 4.649 4.648 4.658; 5.218 5.224 5.228 5.234; 5.224 5.231 5.234 5.241; 6.89 6.902 6.927 6.939; 6.904 6.913 6.943 6.953; 7.34 7.351 7.371 7.382; 7.853 7.868 7.87 7.883; 3.478 3.484 3.487 3.494];
% ppmref_hard=[3.218 3.226 3.228 3.232; 3.237 3.24 3.244 3.249; 3.248 3.254 3.257 3.26; 3.406 3.411 3.414 3.419];
% ppmref_match=[ppmref_match; ppmref_hard];
% peak match for urine peaks
ppmurine_match=[1.466 1.471 1.466 1.471; 1.479 1.483 1.479 1.483; 1.913 1.917 1.912 1.915];%; 4.043 4.057 4.043 4.057
est_tab_renorm=est_tab;
vol_add=[0,20,40,60,80,20];%the last one is reference
sampnames={'0 (urine)','20','40','60','80','ref'};
nsample=6;
sampi_interest=[1:nsample];
sampi_spike=[1:(nsample-1)];
dssrange=[-0.5 0.5];
refspec=ft_mat(nsample,:);
normpeak=[3.034 3.043];
% use one urine peak to normalize samples with urines and then use dss peak to normalize the reference spectra
for sampi=sampi_interest
  sampind=find(est_tab_renorm{:,'simulation'}==sampi);
  loctab=est_tab_renorm(sampind,:);
  if sampi~=6
    refrangind=find(loctab{:,1}>normpeak(1)&loctab{:,1}<normpeak(2));
    a_ref=sum(loctab{refrangind,3});
    loctab{:,3}=loctab{:,3}/a_ref;
  else
    dssrangind=find(loctab{:,1}>dssrange(1)&loctab{:,1}<dssrange(2));
    dss_here=max(loctab{dssrangind,3});
    loctab{:,3}=loctab{:,3}/dss_here*dss_ref;
  end
  if sampi==2
    dssrangind=find(loctab{:,1}>dssrange(1)&loctab{:,1}<dssrange(2));
    dss_ref=max(loctab{dssrangind,3});
  end
  est_tab_renorm(sampind,:)=loctab;
end
% plot spike spectra quantification
quan_methods={'max','sum'};
ppmpairtab=[ppmref_match; ppmurine_match];
pairsize=[size(ppmref_match,1) size(ppmurine_match,1)];
labels=[repmat([1],[pairsize(1),1]); repmat([0],[pairsize(2),1])];
npairs=size(ppmpairtab,1);
specvec=[0 10 10 10 10 50 50];% space settings for plot
quan_record=struct();
for quanmethod=quan_methods
  quanmethod=quanmethod{1};
  quan_coll=[];%peak * samples
  % fetch quantifications
  for ppair_i=1:npairs
    locppm=ppmpairtab(ppair_i,:);
    if labels(ppair_i)==1
      searchppm_reg=[repmat(locppm(3:4),[nsample-1,1]); locppm(1:2)];
    else
      searchppm_reg=[locppm(1:2); repmat(locppm(3:4),[4,1]); nan nan];
    end
    tempvec=[];
    for sampi=1:size(searchppm_reg,1)
      locreg=searchppm_reg(sampi,:);
      if isnan(locreg(1))
        tempvec=[tempvec 0];
      else
        mask_samp=est_tab_renorm{:,'simulation'}==sampi;
        mask_ppm=est_tab_renorm{:,'PPM'}>locreg(1) & est_tab_renorm{:,'PPM'}<locreg(2);
        Avec=est_tab_renorm{mask_samp&mask_ppm,'A'};
        if length(Avec)==0
          tempvec=[tempvec 0];
        elseif strcmp(quanmethod,'max')
          tempvec=[tempvec max(Avec)];
        else
          tempvec=[tempvec sum(Avec)];
        end
      end
    end
    quan_coll=[quan_coll; tempvec];
    % scatter for each pair
    locxvec=sampi_interest;
    locyvec=tempvec;
    h=figure();
    spikeind=sampi_spike;
    dlm=fitlm(locxvec(spikeind),locyvec(spikeind));
    lineestpara=dlm.Coefficients.Estimate';
    gscatter(locxvec,locyvec,[],[],[],[20]);
    refline(lineestpara(2),lineestpara(1));
    if ppair_i>pairsize(1)
      titadd='urine';
    else
      titadd='spikein';
    end
    titlefig=[titadd ' ' num2str(locppm(1)) ' ' num2str(locppm(2))];
    xlabel('spike volumn');
    ylabel('estimation');
    title([titlefig ' quantificaiton estimation']);
    set(gca,'XTick',1:6,'XTickLabel',sampnames);
    saveas(h,['scatter_eval.' quanmethod '_' num2str(ppair_i)  '_'  titlefig '.fig']);
    close(h);
  end
  for regi=1:size(exampreg,1)
    ppm_sele_range=exampreg(regi,:);
    ppmbounds=sort(matchPPMs(ppm_sele_range,ppm));
    ind_sele=ppmbounds(1):ppmbounds(2);
    loc_para_tab=est_tab_renorm(est_tab_renorm{:,'PPM'}>ppm_sele_range(1)&est_tab_renorm{:,'PPM'}<ppm_sele_range(2),:);
    stackmat=[];
    for sampi=[nsample:-1:1]
      tab_reg=loc_para_tab(loc_para_tab{:,'simulation'}==sampi,:);
      tab_reg{:,1}=(tab_reg{:,1}-para_add_list.conv_f(1))*para_add_list.conv_f(2);
      sumsig=sin_mixture_simu(tab_reg{:,[1 2 3 4]},timevec_sub_front',0.0,'complex');
      scalfactor=0.5;
      sumsig(1)=sumsig(1)*scalfactor;
      sumsig=[zeros([1,shifttimeadd]) sumsig];
      spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig)',imag(sumsig)'),libdir,'temp');
      stackmat=[stackmat; spec_new_sum{ind_sele,2}'];
    end
    refft=refspec(ind_sele);
    stackmat=[refft/max(refft)*max(stackmat(1,:)); stackmat];
    stackmat=flip(stackmat,1);
    stackspec_time(stackmat,ppm(ind_sele)',0.0,5,['example_stack_plot'],'timeVect',specvec);
    fig=gcf;
    saveas(fig,['stack_example_region_sele' quanmethod num2str(regi) '.fig']);
    close all;
    % line plot for the regions
    ind_reg=find(all(ppmpairtab(:,1:2)>ppm_sele_range(1) & ppmpairtab(:,1:2)<ppm_sele_range(2),2));
    quan_coll_loc=quan_coll(ind_reg,:);
    %
    locxvec=repelem(sampi_interest,length(ind_reg));
    locyvec=quan_coll_loc(:);
    pairsize_loc=length(find(ind_reg<=pairsize(1)));
    pairsize_loc=[pairsize_loc length(ind_reg)-pairsize_loc(1)];
    colorvec=repmat([repmat({'spike'},[pairsize_loc(1),1]);repmat({'urine'},[pairsize_loc(2),1])],[nsample,1]);
    % plot
    h=figure();
    linecoll=[];
    x_spikein=repelem(sampi_spike,length(ind_reg));
    y_spikein=quan_coll_loc(:,sampi_spike);
    y_spikein=y_spikein(:);
    color_spikein=repmat([repmat({'spike'},[pairsize_loc(1),1]);repmat({'urine'},[pairsize_loc(2),1])],[nsample-1,1])';
    % lines by spike-in or urine
    % for color_group=unique(color_spikein)
    %   groupind=find(cellfun(@(x) strcmp(x,color_group),color_spikein));
    %   % linear line
    %   dlm=fitlm(x_spikein(groupind),y_spikein(groupind));
    %   % % plot the line
    %   % subh=plot(dlm);
    %   % delete(subh([1,3,4]));
    %   linecoll=[linecoll; dlm.Coefficients.Estimate'];
    % end
    % lines by each peak
    for rowi=1:size(quan_coll_loc,1)
      seqind=sampi_spike;
      dlm=fitlm(seqind,quan_coll_loc(rowi,seqind));
      linecoll=[linecoll; dlm.Coefficients.Estimate'];
    end
    % plot points
    gscatter(locxvec,locyvec,colorvec,[],[],[20]);
    for linei=1:size(linecoll,1)
      hline=refline(linecoll(linei,2),linecoll(linei,1));
    end
    xlabel('spike volumn');
    ylabel('estimation');
    title(['quantificaiton estimation' num2str(ppm_sele_range(1)) ' ' num2str(ppm_sele_range(2))]);
    set(gca,'XTick',1:6,'XTickLabel',sampnames);
    saveas(h,['scatter_eval.' quanmethod num2str(regi) '_linepoints.fig']);
    close(h);
  end
  % scale each features (sample 3 spike 2 set as 2)
  quan_coll_s2=quan_coll./quan_coll(:,3)*2;
  % scatter plot
  locxvec=repelem(sampi_interest,npairs);
  locyvec=quan_coll_s2(:);
  colorvec=repmat([repmat({'spike'},[pairsize(1),1]);repmat({'urine'},[pairsize(2),1])],[nsample,1]);
  %
  h=figure();
  gscatter(locxvec,locyvec,colorvec,[],[],[20]);
  xlabel('spike volumn');
  ylabel('estimation');
  title(['quantificaiton estimation']);
  set(gca,'XTick',1:6,'XTickLabel',sampnames);
  saveas(h,['scatter_eval.' quanmethod '.fig']);
  close(h);
  % scale by range
  quan_coll_shift=quan_coll;
  quan_coll_shift(:,sampi_spike)=quan_coll_shift(:,sampi_spike)-quan_coll_shift(:,1);
  quan_coll_s_scale=quan_coll_shift./quan_coll_shift(:,(nsample-1))*4;
  quan_coll_s_scale=quan_coll_s_scale(1:pairsize(1),:);
  % scatter plot
  locxvec=repelem(sampi_interest,pairsize(1));
  locyvec=quan_coll_s_scale(:);
  colorvec=repmat([repmat({'spike'},[pairsize(1),1])],[nsample,1]);
  %
  h=figure();
  gscatter(locxvec,locyvec,colorvec,[],[],[20]);
  xlabel('spike volumn');
  ylabel('estimation');
  title(['quantificaiton estimation']);
  set(gca,'XTick',1:6,'XTickLabel',sampnames);
  saveas(h,['scatter_eval.' quanmethod 'scale2.fig']);
  close(h);
  % plot with linear lines
  seleseq=1:pairsize(1);
  [intensityscale ordind]=sort(quan_coll(seleseq,nsample-1));
  cmap=jet(length(intensityscale));
  cinfo=customColormap(intensityscale,'colors',cmap(1:end,:));
  cmap=cinfo.rgb;
  %
  bvec=[];
  fig=figure();
  hold on;
  for rowindi=1:length(ordind)
    rowind=ordind(rowindi);
    locxvec=sampi_interest;
    locyvec=quan_coll(rowind,:);
    gscatter(locxvec,locyvec,[],cmap(rowindi,:),[],[20]);
    groupind=sampi_spike;
    % linear line
    dlm=fitlm(locxvec(groupind),locyvec(groupind));
    para=dlm.Coefficients.Estimate'
    hline=refline(para(2),para(1));
    hline.Color=cmap(rowindi,:);
    bvec=[bvec para(1)];
  end
  xlabel('spike volumn');
  ylabel('estimation');
  title(['quantificaiton estimation']);
  set(gca,'XTick',sampi_interest,'XTickLabel',sampnames);
  % make a separate colorbar
  colormap(cinfo.cmap);
  t=colorbar;
  ticks=intensityscale;
  set(t,'Ticks',(ticks-min(ticks))/range(ticks));
  set(t,'TickLabels',num2cell(round(ticks,2)));
  set(get(t,'ylabel'),'String','intensity');
  set(gca,'FontSize',20);
  saveas(fig,['scatter_eval.' quanmethod 'all_linepoints.fig']);
  close(fig);
  quan_record.(quanmethod)=quan_coll;
end
% spaced whole stack plot
specvec=[0 10 10 10 10 50 50];
ppm_sele_range=[0.5 9.0];
remreg=[4.681 4.887];%the water signals
ppmbounds=sort(matchPPMs(ppm_sele_range,ppm));
ind_sele=ppmbounds(1):ppmbounds(2);
loc_para_tab=est_tab_renorm(est_tab_renorm{:,'PPM'}>ppm_sele_range(1)&est_tab_renorm{:,'PPM'}<ppm_sele_range(2),:);
loc_para_tab=loc_para_tab(loc_para_tab{:,1}<remreg(1)|loc_para_tab{:,1}>remreg(2),:);
stackmat=[];
for sampi=[nsample:-1:1]
  tab_reg=loc_para_tab(loc_para_tab{:,'simulation'}==sampi,:);
  tab_reg{:,1}=(tab_reg{:,1}-para_add_list.conv_f(1))*para_add_list.conv_f(2);
  sumsig=sin_mixture_simu(tab_reg{:,[1 2 3 4]},timevec_sub_front',0.0,'complex');
  scalfactor=0.5;
  sumsig(1)=sumsig(1)*scalfactor;
  sumsig=[zeros([1,shifttimeadd]) sumsig];
  spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig)',imag(sumsig)'),libdir,'temp');
  stackmat=[stackmat; spec_new_sum{ind_sele,2}'];
end
refft=refspec(ind_sele);
stackmat=[refft/max(refft)*max(stackmat(1,:)); stackmat];
stackmat=flip(stackmat,1);
stackspec_time(stackmat,ppm(ind_sele)',0.0,5,['example_stack_plot'],'timeVect',specvec);
fig=gcf;
saveas(fig,['stack_example_region_all.ps']);
close all;

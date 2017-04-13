function [nidmfile, prov] = spm_results_nidm(SPM,xSPM,TabDat,opts)
% Export SPM stats results using the Neuroimaging Data Model (NIDM)
% FORMAT [nidmfile, prov] = spm_results_nidm(SPM,xSPM,TabDat,opts)
% SPM      - structure containing analysis details (see spm_spm.m)
% xSPM     - structure containing inference details (see spm_getSPM.m)
% TabDat   - structure containing results details (see spm_list.m)
% opts     - structure containing extra information about:
%   .group - subject/group(s) under study
%   .mod   - data modality
%   .space - reference space
%
% nidmfile - output NIDM zip archive filename
% prov     - provenance object (see spm_provenance.m)
%__________________________________________________________________________
% References:
% 
% Neuroimaging Data Model (NIDM):
%   http://nidm.nidash.org/
%
% PROV-DM: The PROV Data Model:
%   http://www.w3.org/TR/prov-dm/
%__________________________________________________________________________
% Copyright (C) 2013-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_results_nidm.m 6903 2016-10-12 11:36:41Z guillaume $


%-Get input parameters, interactively if needed
%==========================================================================
if nargin && ischar(SPM) && strcmpi(SPM,'upload')
    if nargin > 2
        upload_to_neurovault(xSPM,TabDat); % token
    else
        upload_to_neurovault(xSPM);
    end
    return;
end

if nargin < 1
    [SPM,xSPM] = spm_getSPM;
elseif nargin < 2
    if isstruct(SPM)
        xSPM = struct('swd',SPM.swd);
    else
        xSPM = struct('swd',spm_file('cpath',SPM));
    end
    [SPM,xSPM] = spm_getSPM(xSPM);
end
if nargin < 3
    % Consider Inf local maxima more than 0mm apart (i.e. all)
    TabDat = spm_list('Table',xSPM,Inf,0);
end
if nargin < 4
    opts = struct;
end

nidm_json = containers.Map();

%-Options
%==========================================================================

%-General options
%--------------------------------------------------------------------------
gz           = '.gz';                        %-Compressed NIfTI {'.gz', ''}
NIDMversion  = '1.3.0';
SVNrev       = '$Rev: 6903 $';

%-Reference space
%--------------------------------------------------------------------------
if ~isfield(opts,'space')
    s1 = {'Subject space','Normalised space (using Segment)',...
        'Normalised space (using Old Segment)','Custom space',...
        'Other normalised MNI space','Other normalised Talairach space'};
    s2 = {'subject','ixi','icbm','custom','mni','talairach'};
    opts.space = spm_input('Reference space :',1,'m',s1);
    opts.space = s2{opts.space};
end
switch opts.space
    case 'subject'
        coordsys = 'nidm_SubjectCoordinateSystem';
    case 'ixi'
        coordsys = 'nidm_Ixi549CoordinateSystem';
    case 'icbm'
        coordsys = 'nidm_IcbmMni152LinearCoordinateSystem';
    case 'custom'
        coordsys = 'nidm_CustomCoordinateSystem';        
    case 'mni'
        coordsys = 'nidm_MNICoordinateSystem';        
    case 'talairach'
        coordsys = 'nidm_TalairachCoordinateSystem';  
    otherwise
        error('Unknown reference space.');
end
nidm_json('nidm_CoordinateSpace/nidm_inWorldCoordinateSystem') = coordsys;

%-Data modality
%--------------------------------------------------------------------------
MRIProtocol  = '';
if ~isfield(opts,'mod')
    m1 = {'Anatomical MRI','functional MRI','Diffusion MRI','PET','SPECT','EEG','MEG'};
    m2 = {'AMRI','FMRI','DMRI','PET','SPECT','EEG','MEG'};
    opts.mod = spm_input('Data modality :','+1','m',m1);
    opts.mod = m2{opts.mod};
end
switch opts.mod
    case 'AMRI'
        ImagingInstrument      = 'nlx_Magneticresonanceimagingscanner';
        
        MRIProtocol            = 'nlx_AnatomicalMRIprotocol';    
    case 'FMRI'
        ImagingInstrument      = 'nlx_Magneticresonanceimagingscanner';
        
        MRIProtocol            = 'nlx_FunctionalMRIprotocol';
    case 'DMRI'
        ImagingInstrument      = 'nlx_Magneticresonanceimagingscanner';
       
        MRIProtocol            = 'nlx_DiffusionMRIprotocol';
    case 'PET'
        ImagingInstrument      = 'nlx_Positronemissiontomographyscanner';
    case 'SPECT'
        ImagingInstrument      = 'nlx_Singlephotonemissioncomputedtomographyscanner';
    case 'EEG'
        ImagingInstrument      = 'nlx_Electroencephalographymachine';
    case 'MEG'
        ImagingInstrument      = 'nlx_Magnetoencephalographymachine';
    otherwise
        error('Unknown modality.');
end

%-Subject/Group(s)
%--------------------------------------------------------------------------
if ~isfield(opts,'group')
    opts.group.N = spm_input('Number of subjects per group :','+1','e');
    if isequal(opts.group.N,1)
        opts.group.name = {'single subject'};
    else
        for i=1:numel(opts.group.N)
            opts.group.name{i} = spm_input(...
                sprintf('Name of group %d :',i),'+1','s');
        end
    end
end
groups = opts.group;

%-Units
%--------------------------------------------------------------------------
try
    units = xSPM.units;
catch
    try
        units = SPM.xVol.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
end

%==========================================================================
%-Populate output directory
%==========================================================================
if ~exist(SPM.swd,'dir'), SPM.swd = pwd; end

%-Design Matrix values
%--------------------------------------------------------------------------
nidm_json('nidm_DesignMatrix/prov:value') = SPM.xX.xKXs.X;
nidm_json('nidm_DesignMatrix/nidm_regressorNames') = SPM.xX.name;

nidm_json('nidm_CoordinateSpace/units') = units;

%-Maximum Intensity Projection image (as png)
%--------------------------------------------------------------------------
temp_img = {};
files.mip    = fullfile(SPM.swd,'MaximumIntensityProjection.png');
MIP          = spm_mip(xSPM.Z,xSPM.XYZmm,xSPM.M,units);
imwrite(MIP,gray(64),files.mip,'png');

temp_img{end+1} = files.mip;

%-Beta images (as NIfTI)
%--------------------------------------------------------------------------
for i=1:numel(SPM.Vbeta)
%     files.beta{i} = fullfile(outdir,[sprintf('ParameterEstimate_%04d',i) '.nii' gz]);
    files.beta{i} = fullfile(xSPM.swd,SPM.Vbeta(i).fname);
%     img2nii(fullfile(xSPM.swd,SPM.Vbeta(i).fname), files.beta{i});
end

%-SPM{.}, contrast, contrast standard error, and contrast explained mean square images (as NIfTI)
%-----------------------x---------------------------------------------------
for i=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%04d',i); end
%     files.spm{i}  = fullfile(outdir,[xSPM.STAT 'Statistic' postfix '.nii' gz]);
%     img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vspm.fname), files.spm{i}, xSPM);
    files.spm{i} = fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vspm.fname);
    if xSPM.STAT == 'T'
%         files.con{i} = fullfile(outdir,['Contrast' postfix '.nii' gz]);
%         img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname), files.con{i},...
%             struct('STAT','con'));
        files.con{i} = fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname);
        
        % Create conse map
        files.conse{i} = fullfile(xSPM.swd,['ContrastStandardError' postfix '.nii' gz]);
        Vc = SPM.xCon(xSPM.Ic(i)).c' * SPM.xX.Bcov * SPM.xCon(xSPM.Ic(i)).c;
        img2nii(fullfile(xSPM.swd,SPM.VResMS.fname), files.conse{i}, struct('fcn',@(x) sqrt(x*Vc)));
        
        temp_img{end+1} = files.conse{i};
        
    elseif xSPM.STAT == 'F'
        files.effms{i} = fullfile(xSPM.swd,['ContrastExplainedMeanSquare' postfix '.nii' gz]);
        eidf = SPM.xCon(xSPM.Ic(i)).eidf;
        img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname), files.effms{i}, struct('fcn',@(x) x/eidf));
        
        temp_img{end+1} = files.effms{i};
    end
end

%-Thresholded SPM{.} image (as NIfTI)
%--------------------------------------------------------------------------
files.tspm = fullfile(xSPM.swd,['ExcursionSet.nii' gz]);
if ~isempty(gz), files.tspm = spm_file(files.tspm,'ext',''); end
spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'',files.tspm);
if ~isempty(gz), gzip(files.tspm); spm_unlink(files.tspm); files.tspm = [files.tspm gz]; end

temp_img{end+1} = files.tspm;

%-Residual Mean Squares image (as NIfTI)
%--------------------------------------------------------------------------
files.resms = fullfile(xSPM.swd,SPM.VResMS.fname);
% files.resms = fullfile(outdir,['ResidualMeanSquares.nii' gz]);
% img2nii(fullfile(xSPM.swd,SPM.VResMS.fname), files.resms);

%-Resels per Voxel image (as NIfTI)
%--------------------------------------------------------------------------
files.rpv = fullfile(xSPM.swd,SPM.xVol.VRpv.fname);
% files.rpv = fullfile(outdir,['ReselsPerVoxel.nii' gz]);
% img2nii(fullfile(xSPM.swd,SPM.xVol.VRpv.fname), files.rpv);

%-Analysis mask image (as NIfTI)
%--------------------------------------------------------------------------
files.mask = fullfile(xSPM.swd,SPM.VM.fname);
% files.mask = fullfile(outdir,['Mask.nii' gz]);
% img2nii(fullfile(xSPM.swd,SPM.VM.fname), files.mask);

%-Grand mean image (as NIfTI)
%--------------------------------------------------------------------------
files.grandmean = fullfile(xSPM.swd,'GrandMean.nii');

sf  = mean(SPM.xX.xKXs.X,1);
Vb  = SPM.Vbeta;
for i=1:numel(Vb), Vb(i).pinfo(1:2,:) = Vb(i).pinfo(1:2,:) * sf(i); end
Vgm = struct(...
    'fname',   files.grandmean,...
    'dim',     Vb(1).dim,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     Vb(1).mat,...
    'pinfo',   [1 0 0]',...
    'descrip', 'Grand Mean');
Vgm = spm_create_vol(Vgm);
Vgm.pinfo(1,1) = spm_add(Vb,Vgm);
Vgm = spm_create_vol(Vgm);
if ~isempty(gz), gzip(files.grandmean); spm_unlink(files.grandmean); files.grandmean = [files.grandmean gz]; end
temp_img{end+1} = files.grandmean;

%-Explicit mask image (as NIfTI)
%--------------------------------------------------------------------------
if ~isempty(SPM.xM.VM)
    files.emask = fullfile(xSPM.swd,['CustomMask.nii' gz]);
    temp_img{end+1} = files.emask;
    
    if isempty(spm_file(SPM.xM.VM.fname,'path'))
        Vem = fullfile(xSPM.swd,SPM.xM.VM.fname);
    else
        Vem = SPM.xM.VM.fname;
    end
    img2nii(Vem, files.emask);
else
    files.emask = '';
end

%-Clusters n-ary image (as NIfTI)
%--------------------------------------------------------------------------
files.clust = fullfile(xSPM.swd,['ClusterLabels.nii' gz]);
temp_img{end+1} = files.clust;
if ~isempty(gz), files.clust = spm_file(files.clust,'ext',''); end
Z   = spm_clusters(xSPM.XYZ);
idx = find(~cellfun(@isempty,{TabDat.dat{:,5}}));
n   = zeros(1,numel(idx));
for i=1:numel(idx)
    [unused,j] = spm_XYZreg('NearestXYZ',TabDat.dat{idx(i),12}',xSPM.XYZmm);
    n(i) = Z(j);
end
n(n) = 1:numel(n);
if max(Z) ~= numel(idx)
    warning('Small Volume Correction not handled yet.');
    n(numel(idx)+1:max(Z)) = 0;
end
Z    = n(Z);
spm_write_filtered(Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'',files.clust);
if ~isempty(gz), gzip(files.clust); spm_unlink(files.clust); files.clust = [files.clust gz]; end

%-Display mask images (as NIfTI)
%--------------------------------------------------------------------------
for i=1:numel(xSPM.Im)
    files.dmask{i} = fullfile(xSPM.swd,[sprintf('DisplayMask_%04d.nii',i) gz]);
    temp_img{end+1} = files.dmask{i};
    if isnumeric(xSPM.Im)
        um = spm_u(xSPM.pm,[SPM.xCon(xSPM.Im(i)).eidf,SPM.xX.erdf],...
            SPM.xCon(xSPM.Im(i)).STAT);
        if ~xSPM.Ex, fcn = @(x) x > um;
        else         fcn = @(x) x <= um; end
        img2nii(SPM.xCon(xSPM.Im(i)).Vspm.fname, files.dmask{i}, struct('fcn',fcn));
    else
        if ~xSPM.Ex, fcn = @(x) x~=0 & ~isnan(x);
        else         fcn = @(x) ~(x~=0 & ~isnan(x)); end
        img2nii(xSPM.Im{i}, files.dmask{i}, struct('fcn',fcn));
    end
end
if numel(xSPM.Im) == 0, files.dmask = {}; end

%-SVC Mask (as NIfTI)
%--------------------------------------------------------------------------
if strcmp(TabDat.tit,'p-values adjusted for search volume')
    files.svcmask = '';
elseif strncmp(TabDat.tit,'search volume: ',15)
    warning('Small Volume Correction not handled yet.'); % see spm_VOI.m
    % '%0.1fmm sphere at [%.0f,%.0f,%.0f]'
    % '%0.1f x %0.1f x %0.1f mm box at [%.0f,%.0f,%.0f]'
    % 'image mask: %s'
    files.svcmask = '';
else
    warning('Unable to retrieve search volume details: assuming whole brain search.');
    files.svcmask = '';
end

%-Search Space mask image (as NIfTI)
%--------------------------------------------------------------------------
files.searchspace = fullfile(xSPM.swd,SPM.VM.fname);
% files.searchspace = fullfile(outdir,['SearchSpaceMask.nii' gz]);
% img2nii(fullfile(xSPM.swd,SPM.VM.fname), files.searchspace);


%==========================================================================
%-                          D A T A   M O D E L
%==========================================================================

clear coordspace originalfile isHumanReadable

niifmt = {'image/nifti','xsd:string'};
isHumanReadable(false);

%-Provenance
%--------------------------------------------------------------------------
[V,R] = spm('Ver');

nidm_json('nidm_NIDMResultsExporter/prov:type') = 'spm_results_nidm';
nidm_json('nidm_NIDMResultsExporter/nidm_softwareVersion') = [V(4:end) '.' char(regexp(SVNrev,'\$Rev: (\w.*?) \$','tokens','once'))];
nidm_json('nidm_NIDMResults/nidm_version') = NIDMversion;

%-Agent: SPM
%--------------------------------------------------------------------------
nidm_json('nidm_NeuroimagingAnalysisSoftware/prov:type') = 'SPM';
nidm_json('nidm_NeuroimagingAnalysisSoftware/nidm_softwareVersion') = [V(4:end) '.' R];

%-Entity: Coordinate Space
%--------------------------------------------------------------------------
% id_data_coordspace = coordspace(p,xSPM.M,xSPM.DIM,units,coordsys,1);

%-Agent: Scanner
%--------------------------------------------------------------------------
nidm_json('nlx_Imaginginstrument/prov:type') = ImagingInstrument;

%-Agent: Person
%--------------------------------------------------------------------------
if ~isequal(groups.N,1)
    %-Agent: Group
    %----------------------------------------------------------------------
    nidm_groups = containers.Map();
    for i=1:numel(groups.N)
        nidm_groups(groups.name{i}) = groups.N(i);
    end
    nidm_json('Groups') = nidm_groups;
end

%-Entity: Image Data
%--------------------------------------------------------------------------
if isfield(SPM,'Sess')
    nidm_json('nidm_Data/nidm_grandMeanScaling') = true;
    nidm_json('nidm_Data/nidm_targetIntensity') = SPM.xGX.GM;
else
    nidm_json('nidm_Data/nidm_grandMeanScaling') = false;
end
if ~isempty(MRIProtocol)
    nidm_json('nidm_Data/nidm_hasMRIProtocol') = MRIProtocol;
end

%-Entity: Drift Model
%--------------------------------------------------------------------------
if isfield(SPM,'Sess') && isfield(SPM.xX,'K')
    nidm_json('nidm_DesignMatrix/nidm_hasDriftModel') = 'spm_DiscreteCosineTransformbasisDriftModel';
    nidm_json('nidm_DesignMatrix/spm_SPMsDriftCutoffPeriod') = SPM.xX.K(1).HParam;
end

%-Entity: Design Matrix
%--------------------------------------------------------------------------
if isfield(SPM,'xBF')
    switch SPM.xBF.name
        case 'hrf'
            nidm_json('nidm_DesignMatrix/nidm_hasHRFBasis') = 'spm_SPMsCanonicalHRF';
        case 'hrf (with time derivative)'
            nidm_json('nidm_DesignMatrix/nidm_hasHRFBasis') = ['spm_SPMsCanonicalHRF', 'spm_SPMsTemporalDerivative'];
        case 'hrf (with time and dispersion derivatives)'
            nidm_json('nidm_DesignMatrix/nidm_hasHRFBasis') = ['spm_SPMsCanonicalHRF', 'spm_SPMsTemporalDerivative', 'spm_SPMsDispersionDerivative'];
        case 'Finite Impulse Response'
            nidm_json('nidm_DesignMatrix/nidm_hasHRFBasis') = ['nidm_FiniteImpulseResponseBasisSet'];
        case 'Fourier set'
            nidm_json('nidm_DesignMatrix/nidm_hasHRFBasis') = ['nidm_FourierBasisSet'];
        case 'Gamma functions'
            nidm_json('nidm_DesignMatrix/nidm_hasHRFBasis') = ['nidm_GammaBasisSet'];
        case {'Fourier set (Hanning)'}
            warning('Not implemented "%s".',SPM.xBF.name);
        otherwise
            warning('Unknown basis set.');
    end
end

%-Entity: Explicit Mask
%--------------------------------------------------------------------------
if ~isempty(SPM.xM.VM)
    nidm_json('nidm_CustomMap/prov:atLocation') = files.emask;
end

%-Entity: Error Model
%--------------------------------------------------------------------------
if isfield(SPM.xVi,'form')
    if strcmp(SPM.xVi.form,'i.i.d')
        nidm_json('nidm_ErrorModel/nidm_errorVarianceHomogeneous') = true;
        nidm_json('nidm_ErrorModel/nidm_hasErrorDependence') = 'nidm_IndependentError';
        nidm_json('nidm_ModelParameterEstimation/nidm_withEstimationMethod') = 'obo_ordinaryleastsquaresestimation';
    else
        nidm_json('nidm_ErrorModel/nidm_errorVarianceHomogeneous') = true;
        nidm_json('nidm_ErrorModel/nidm_hasErrorDependence') = 'obo_Toeplitzcovariancestructure';
        nidm_json('nidm_ErrorModel/nidm_dependenceMapWiseDependence') = 'nidm_ConstantParameter';
        nidm_json('nidm_ErrorModel/nidm_varianceMapWiseDependence') = 'nidm_IndependentParameter';
        nidm_json('nidm_ModelParameterEstimation/nidm_withEstimationMethod') = 'obo_generalizedleastsquaresestimation';
    end
else
    if ~isfield(SPM.xVi,'Vi') || numel(SPM.xVi.Vi) == 1 % assume it's identity
        nidm_json('nidm_ErrorModel/nidm_errorVarianceHomogeneous') = true;
        nidm_json('nidm_ErrorModel/nidm_hasErrorDependence') = 'nidm_IndependentError';
        nidm_json('nidm_ErrorModel/nidm_varianceMapWiseDependence') = 'nidm_IndependentParameter';
        nidm_json('nidm_ModelParameterEstimation/nidm_withEstimationMethod') = 'obo_ordinaryleastsquaresestimation';
    else
        nidm_json('nidm_ErrorModel/nidm_errorVarianceHomogeneous') = false;
        nidm_json('nidm_ErrorModel/nidm_hasErrorDependence') = 'obo_unstructuredcovariancestructure';
        nidm_json('nidm_ErrorModel/nidm_dependenceMapWiseDependence') = 'nidm_ConstantParameter';
        nidm_json('nidm_ErrorModel/nidm_varianceMapWiseDependence') = 'nidm_IndependentParameter';
        nidm_json('nidm_ModelParameterEstimation/nidm_withEstimationMethod') = 'obo_generalizedleastsquaresestimation';
    end
end

nidm_json('nidm_ErrorModel/nidm_hasErrorDistribution') = 'obo_normaldistribution';

%-Activity: Model Parameters Estimation
%==========================================================================

%-Entity: Mask Map
%--------------------------------------------------------------------------
nidm_json('nidm_MaskMap/prov:atLocation') = uri(spm_file(files.mask,'cpath'));

%-Entity: Grand Mean Map
%--------------------------------------------------------------------------
nidm_json('nidm_GrandMeanMap/prov:atLocation') = uri(spm_file(files.grandmean,'cpath'));

%-Entity: Parameter Estimate (Beta) Maps
%--------------------------------------------------------------------------
nidm_betas = containers.Map();
for i=1:numel(SPM.Vbeta)
    nidm_betas(nidm_esc(SPM.xX.name{i})) = uri(files.beta{i});
end
nidm_json('ParameterEstimateMaps') = nidm_betas;

%-Entity: ResMS Map
%--------------------------------------------------------------------------
nidm_json('nidm_ResidualMeanSquaresMap/prov:atLocation') = uri(spm_file(files.resms,'cpath'));

%-Entity: RPV Map
%--------------------------------------------------------------------------
nidm_json('nidm_ReselsPerVoxelMap/prov:atLocation') = uri(spm_file(files.rpv,'cpath'));

%-Activity: Contrast Estimation
%==========================================================================
STAT = xSPM.STAT;
if STAT == 'T', STAT = lower(STAT); end

nidm_contrasts = containers.Map();
contrast_names = {};

for c=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%d',c); end
       
    if xSPM.STAT == 'T'
        con_name = nidm_esc(SPM.xCon(xSPM.Ic(c)).name);
        contrast_names{end+1} = con_name;
        nidm_contrasts(con_name) = ...
            containers.Map({...
                'obo_contrastweightmatrix/prov:value', ...
                'nidm_StatisticMap/nidm_statisticType', ...
                'nidm_StatisticMap/nidm_errorDegreesOfFreedom', ...
                'nidm_StatisticMap/prov:atLocation', ...
                'nidm_ContrastMap/prov:atLocation', ...
                'nidm_ContrastStandardErrorMap/prov:atLocation' ...
            }, ...
            { ...
                SPM.xCon(xSPM.Ic(c)).c', ...
                ['obo_' STAT 'statistic'], ...
                xSPM.df(2), ...
                uri(spm_file(files.spm{c},'cpath')), ...
                uri(spm_file(files.con{c},'cpath')), ...
                uri(spm_file(files.conse{c},'cpath')) ...
            });
    end
    if xSPM.STAT == 'F'
        con_name = nidm_esc(SPM.xCon(xSPM.Ic(c)).name);
        contrast_names{end+1} = con_name;
        nidm_contrasts(con_name) = ...
            containers.Map({...
                'obo_contrastweightmatrix/prov:value', ...
                'nidm_StatisticMap/nidm_statisticType', ...
                'nidm_StatisticMap/nidm_errorDegreesOfFreedom', ...
                'nidm_StatisticMap/nidm_effectDegreesOfFreedom', ...
                'nidm_StatisticMap/prov:atLocation', ...
                'nidm_ContrastExplainedMeanSquareMap/prov:atLocation' ...
            }, ...
            { ...
                SPM.xCon(xSPM.Ic(c)).c', ...
                ['obo_' STAT 'statistic'], ...
                xSPM.df(2), ...
                xSPM.df(1), ...
                uri(spm_file(files.spm{c},'cpath')), ...
                uri(uri(spm_file(files.effms{c},'cpath'))) ...
            });
    end
end
nidm_json('Contrasts') = nidm_contrasts;

%-Entity: Height Threshold
%--------------------------------------------------------------------------
thresh(1).type  = 'obo_statistic';
thresh(1).label = 'Height Threshold';
thresh(1).value = xSPM.u; % TabDat.ftr{1,2}(1)
thresh(2).type  = 'nidm_PValueUncorrected';
thresh(2).label = 'Height Threshold';
thresh(2).value = TabDat.ftr{1,2}(2);
thresh(3).type  = 'obo_FWERadjustedpvalue';
thresh(3).label = 'Height Threshold';
thresh(3).value = TabDat.ftr{1,2}(3);
td = regexp(xSPM.thresDesc,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
if isempty(td)
    td = regexp(xSPM.thresDesc,'\w=(?<u>[\.\d]+)','names');
    if ~isempty(td)
        thresh_order = 1:3; % Statistic
        thresh_desc  = sprintf(': %s=%f)',xSPM.STAT,xSPM.u);
    else
        warning('Unkwnown threshold type.');
        thresh_order = 1:3; % unknown
        thresh_desc  = '';
    end
else
    switch td.thresDesc
        case 'FWE'
            thresh_order = [3 1 2]; % FWE
            thresh_desc  = sprintf(': p<%f (FWE)',TabDat.ftr{1,2}(3));
        case 'unc.'
            thresh_order = [2 1 3]; % uncorrected
            thresh_desc  = sprintf(': p<%f (unc.)',TabDat.ftr{1,2}(2));
            % Set uncorrected p-value threshold to the user-defined value
            % (to avoid possible floating point approximations)            
            %thresh(2).value = str2double(td.u);
        case 'FDR'
            thresh(3).type  = 'obo_qvalue';
            thresh(3).label = 'Height Threshold';
            thresh(3).value = str2double(td.u);
            thresh_order = [3 1 2]; % FDR
            thresh_desc  = sprintf(': p<%s (FDR)',td.u);
        otherwise
            warning('Unkwnown threshold type.');
            thresh_order = 1:3; % unknown
            thresh_desc  = '';
    end
end
thresh = thresh(thresh_order);
thresh(1).label = [thresh(1).label thresh_desc];

%-Entity: Extent Threshold
%--------------------------------------------------------------------------
if spm_get_defaults('stats.rft.nonstat')
    warning('Non-stationary RFT results not handled yet.');
end
V2R = 1 / prod(xSPM.FWHM);

if TabDat.ftr{2,2}(1) > 0
    kk = [TabDat.ftr{2,2}(2) TabDat.ftr{2,2}(3)];
else
    kk = [1 1];
end

nidm_inference = containers.Map();
nidm_inference('nidm_HeightThreshold/prov:type') = thresh(1).type;
nidm_inference('nidm_HeightThreshold/prov:value') = thresh(1).value;
nidm_inference('nidm_ExtentThreshold/prov:type') = 'obo_statistic';
nidm_inference('nidm_ExtentThreshold/nidm_clusterSizeInVoxels') = TabDat.ftr{2,2}(1);
nidm_inference('nidm_ExtentThreshold/nidm_clusterSizeInResels') = TabDat.ftr{2,2}(1)*V2R;

nidm_height_equivthresh = containers.Map();
nidm_extent_equivthresh = containers.Map();
ex_equiv_types = {'obo_FWERadjustedpvalue', 'nidm_PValueUncorrected'};
for i = 2:3
    nidm_height_equivthresh(nidm_esc([thresh(i).label ' ' num2str(i)])) = containers.Map(...
        {...
            'nidm_HeightThreshold/prov:type', ...
            'nidm_HeightThreshold/prov:value'
        }, {...
            thresh(i).type, ...
            thresh(i).value, ...
        });
    nidm_extent_equivthresh(nidm_esc([thresh(i).label ' ' num2str(i)])) = containers.Map(...
        {...
            'nidm_ExtentThreshold/prov:type', ...
            'nidm_ExtentThreshold/prov:value' % SPM equiv thresholds are always p-values (i.e. not cluster sizes)
        }, {...
            ex_equiv_types{i-1}, ...
            kk(i-1)
        });
end
nidm_inference('nidm_HeightThreshold/nidm_equivalentThreshold') = nidm_height_equivthresh;
nidm_inference('nidm_ExtentThreshold/nidm_equivalentThreshold') = nidm_extent_equivthresh;

nidm_inferences = containers.Map();
%-Entity: Peak & Cluster Definition Criteria
%--------------------------------------------------------------------------
% TabDat.str = 'table shows %d local maxima more than %.1fmm apart'
maxNumberOfPeaksPerCluster = spm_get_defaults('stats.results.volume.nbmax');
minDistanceBetweenPeaks = spm_get_defaults('stats.results.volume.distmin');
clusterConnectivityCriterion = 18; % see spm_max.m

nidm_json('nidm_ClusterDefinitionCriteria/nidm_hasConnectivityCriterion') = ...
    sprintf('nidm_voxel%dconnected',clusterConnectivityCriterion);
nidm_json('nidm_PeakDefinitionCriteria/nidm_minDistanceBetweenPeaks') = ...
    minDistanceBetweenPeaks;
nidm_json('nidm_PeakDefinitionCriteria/nidm_maxNumberOfPeaksPerCluster') = ...
    maxNumberOfPeaksPerCluster;
    
%-Activity: Inference
%==========================================================================
nidm_inference('nidm_Inference/nidm_hasAlternativeHypothesis') = ...
    'nidm_OneTailedTest';
if numel(xSPM.Ic) == 1
    conj = false;
else
    conj = true;
    nidm_inference('nidm_contrastName') = contrast_names;
    if xSPM.n == 1
        nidm_inference('nidm_Inference/prov:type') = 'nidm_ConjunctionInference';
% 
%         st = {'prov:type',nidm_conv('nidm_ConjunctionInference',p), ...
%               nidm_conv('nidm_hasAlternativeHypothesis',p),nidm_conv('nidm_OneTailedTest',p),...
%               'prov:label','Conjunction Inference'};
    else
        nidm_inference('nidm_Inference/prov:type') = 'spm_PartialConjunctionInference';
        nidm_inference('spm_PartialConjunctionInference/spm_partialConjunctionDegree') = xSPM.n;
%         
%         st = {'prov:type',nidm_conv('spm_PartialConjunctionInference',p), ...
%               nidm_conv('nidm_hasAlternativeHypothesis',p),nidm_conv('nidm_OneTailedTest',p),...
%               'prov:label','Partial Conjunction Inference', ...
%               nidm_conv('spm_partialConjunctionDegree',p),{xSPM.n,'xsd:int'}};
    end
end

%-Entity: Display Mask Maps
%--------------------------------------------------------------------------
for i=1:numel(files.dmask)
    % TODO: why more than 1??
    
    nidm_inference('nidm_DisplayMaskMap/prov:atLocation') = ...
        uri(spm_file(files.dmask{i},'cpath'));
end

%-Entity: SVC Mask Map
%--------------------------------------------------------------------------
if ~isempty(files.svcmask)
    nidm_inference('nidm:SubVolumeMap/prov:atLocation') = ...
        uri(spm_file(files.svcmask,'cpath'));    
end

%-Entity: Search Space
%--------------------------------------------------------------------------
nidm_inference('nidm_SearchSpaceMaskMap/prov:atLocation') = ...
        uri(spm_file(files.searchspace,'cpath'));
nidm_inference('nidm_SearchSpaceMaskMap/nidm_searchVolumeInVoxels') = xSPM.S;    
nidm_inference('nidm_SearchSpaceMaskMap/nidm_searchVolumeInUnits') = TabDat.ftr{8,2}(1);    
nidm_inference('nidm_SearchSpaceMaskMap/nidm_reselSizeInVoxels') = TabDat.ftr{9,2}(end);
nidm_inference('nidm_SearchSpaceMaskMap/nidm_searchVolumeInResels') = TabDat.ftr{8,2}(3);
nidm_inference('nidm_SearchSpaceMaskMap/spm_searchVolumeReselsGeometry') = xSPM.R;
nidm_inference('nidm_SearchSpaceMaskMap/nidm_noiseFWHMInVoxels') = xSPM.FWHM;
nidm_inference('nidm_SearchSpaceMaskMap/nidm_noiseFWHMInUnits') = TabDat.ftr{7,2}(1:3);
nidm_inference('nidm_SearchSpaceMaskMap/nidm_randomFieldStationarity') = ...
    ~spm_get_defaults('stats.rft.nonstat');
nidm_inference('nidm_SearchSpaceMaskMap/nidm_expectedNumberOfVoxelsPerCluster') = TabDat.ftr{3,2};
nidm_inference('nidm_SearchSpaceMaskMap/nidm_expectedNumberOfClusters') = TabDat.ftr{4,2};
nidm_inference('nidm_SearchSpaceMaskMap/nidm_heightCriticalThresholdFWE05') = xSPM.uc(1);
nidm_inference('nidm_SearchSpaceMaskMap/nidm_heightCriticalThresholdFDR05') = xSPM.uc(2);

if isfinite(xSPM.uc(3))
    nidm_inference('nidm_SearchSpaceMaskMap/spm_smallestSignificantClusterSizeInVoxelsFWE05') = xSPM.uc(3);
end
if isfinite(xSPM.uc(4))
    nidm_inference('nidm_SearchSpaceMaskMap/spm_smallestSignificantClusterSizeInVoxelsFDR05') = xSPM.uc(4);
end

%-Entity: Excursion Set
%--------------------------------------------------------------------------
if size(TabDat.dat,1) > 0
    c  = TabDat.dat{1,2};
    pc = TabDat.dat{1,1};
else
    c  = 0;
    pc = NaN;
end

nidm_inference('nidm_ExcursionSetMap/prov:atLocation') = uri(spm_file(files.tspm,'cpath'));
nidm_inference('nidm_ExcursionSetMap/nidm_numberOfSupraThresholdClusters') = c;
nidm_inference('nidm_ExcursionSetMap/nidm_pValue') = pc;

nidm_inference('nidm_ClusterLabelsMap/prov:atLocation') = uri(spm_file(files.clust,'cpath'));
nidm_inference('nidm_ExcursionSetMap/nidm_hasMaximumIntensityProjection') = uri(spm_file(files.mip,'cpath'));


%-Entity: Clusters
%--------------------------------------------------------------------------
nidm_clusters = containers.Map();

idx = find(~cellfun(@isempty,{TabDat.dat{:,5}}));

clustidx = cumsum(~cellfun(@isempty,{TabDat.dat{:,5}}));
num_digits = numel(num2str(max(clustidx)));

for i=1:numel(idx)
    clust_key = num2str(i, ['%0' num2str(num_digits) 'd']);
    
    clust_keys = {...
            'nidm_SupraThresholdCluster/nidm_clusterSizeInVoxels',...
            'nidm_SupraThresholdCluster/nidm_clusterSizeInResels',...
            'nidm_SupraThresholdCluster/nidm_pValueUncorrected',...
            'nidm_SupraThresholdCluster/nidm_pValueFWER',...
            'nidm_SupraThresholdCluster/nidm_qValueFDR'...
        };
    clust_values = {...
            TabDat.dat{idx(i),5}, ...
            TabDat.dat{idx(i),5}*V2R, ...
            TabDat.dat{idx(i),6}, ...
            TabDat.dat{idx(i),3}, ...
            TabDat.dat{idx(i),4} ...   
        };
    
%-Entity: Peaks
%--------------------------------------------------------------------------
    nidm_peaks = containers.Map();
    for j=1:size(TabDat.dat,1)
        if clustidx(j) == i
            %     iPeak  = sprintf('%04d',i);   
            nidm_peaks(num2str(j)) = containers.Map(...
                { ...
                    'nidm_Peak/prov:value', ...
                    'nidm_Coordinate/nidm_coordinateVector', ...
                    'nidm_Peak/nidm_pValueUncorrected', ...
                    'nidm_Peak/nidm_equivalentZStatistic', ...
                    'nidm_Peak/nidm_pValueFWER', ...
                    'nidm_Peak/nidm_qValueFDR'
                }, ...
                { ...
                    TabDat.dat{j,9}, ...
                    TabDat.dat{j,12}(1:3), ...
                    TabDat.dat{j,11}, ...
                    xsdfloat(TabDat.dat{j,10}), ...
                    TabDat.dat{j,7}, ...
                    TabDat.dat{j,8} 
                });
        end
    end
    clust_keys{end+1} = 'Peaks';
    clust_values{end+1} = nidm_peaks;
    
    nidm_clusters(clust_key) = containers.Map(clust_keys, clust_values);
end

nidm_inference('Clusters') = nidm_clusters;

nidm_inferences([contrast_names{:}]) = nidm_inference;

if ~conj
    nidm_json('Inferences') = nidm_inferences;
else
    nidm_json('ConjunctionInferences') = nidm_inferences;    
end

[nidmfile, prov] = spm_nidmresults(nidm_json, SPM.swd);

for i = 1:numel(temp_img)
    spm_unlink(temp_img{i});
end

% % pp.bundle(idResults,p);
% 
% %==========================================================================
% %-                  P R O V   S E R I A L I Z A T I O N
% %==========================================================================
% %serialize(pp,fullfile(outdir,'nidm.provn'));
% serialize(pp,fullfile(outdir,'nidm.ttl'));
% try, serialize(pp,fullfile(outdir,'nidm.jsonld')); end
% %serialize(pp,fullfile(outdir,'nidm.json'));
% %serialize(pp,fullfile(outdir,'nidm.pdf'));
% 
% i = 1;
% while true
%     nidmfile = fullfile(SPM.swd,sprintf('spm_%04d.nidm.zip',i));
%     if spm_existfile(nidmfile), i = i + 1; else break; end
% end
% f = zip(nidmfile,'*',outdir);
% for i=1:numel(f)
%     spm_unlink(fullfile(outdir,f{i}));
% end
% rmdir(outdir);
% 
% prov = pp;


%==========================================================================
% function v = xsdfloat(v)
%==========================================================================
function v = xsdfloat(v)
% See http://books.xmlschemata.org/relaxng/ch19-77095.html
if numel(v) == 1 && isinf(v) && v > 0, v = 'INF';  end
if numel(v) == 1 && isinf(v) && v < 0, v = '-INF'; end
if numel(v) == 1 && isnan(v),          v = 'NaN';  end


%==========================================================================
% function str = html_esc(str)
%==========================================================================
function str = html_esc(str)
%-Escape
% See http://www.w3.org/TR/html4/charset.html#h-5.3.2
str = strrep(str,'&','&amp;');
str = strrep(str,'<','&lt;');
str = strrep(str,'>','&gt;');
str = strrep(str,'"','&quot;');


%==========================================================================
% function u = uri(u)
%==========================================================================
function u = uri(u)
%-File URI scheme
%if ispc, s='/'; else s=''; end
%u = ['file://' s strrep(spm_file(u,'cpath'),'\','/')];
e = ' ';
for i=1:length(e)
    u = strrep(u,e(i),['%' dec2hex(e(i))]);
end
u = spm_file(u,'filename');


%==========================================================================
% function checksum = sha512sum(file)
%==========================================================================
function checksum = sha512sum(file)
md   = java.security.MessageDigest.getInstance('SHA-512');
file = spm_file(file,'cpath');
fid  = fopen(file,'rb');
if fid == -1, error('Cannot open "%s".',file); end
md.update(fread(fid,Inf,'*uint8'));
fclose(fid);
checksum = typecast(md.digest,'uint8');
checksum = lower(reshape(dec2hex(checksum)',1,[]));


%==========================================================================
% function checksum = md5sum(data)
%==========================================================================
function checksum = md5sum(data)
if ~nargin
    data = char(java.util.UUID.randomUUID);
end
md   = java.security.MessageDigest.getInstance('MD5');
if ischar(data)
    md.update(uint8(data));
else
    md.update(typecast(data,'uint8'));
end
checksum = typecast(md.digest,'uint8');
checksum = lower(reshape(dec2hex(checksum)',1,[]));


%==========================================================================
% function img2nii(img,nii,xSPM)
%==========================================================================
function img2nii(img,nii,xSPM)
if nargin == 2, xSPM = struct; end
if ~isfield(xSPM,'STAT'), xSPM.STAT = ''; end
if ~isfield(xSPM,'fcn'), xSPM.fcn = @(x) x; end
if nargin == 1, nii = spm_file(img,'ext','.nii'); end
gz = strcmp(spm_file(nii,'ext'),'gz');
if gz, nii = spm_file(nii,'ext',''); end
ni     = nifti(img);
no     = nifti;
no.dat = file_array(nii,...
                    ni.dat.dim,...
                    ni.dat.dtype,...
                    0,...
                    ni.dat.scl_slope,...
                    ni.dat.scl_inter);
no.mat  = ni.mat;
no.mat_intent = ni.mat_intent;
no.mat0 = ni.mat0;
no.mat0_intent = ni.mat0_intent;
no.descrip = ni.descrip;
switch xSPM.STAT
    case 'T'
        no.intent.name  = ['spm' xSPM.STATstr];
        no.intent.code  = 3;
        no.intent.param = xSPM.df(2);
    case 'F'
        no.intent.name  = ['spm' xSPM.STATstr];
        no.intent.code  = 4;
        no.intent.param = xSPM.df;
    case 'con'
        no.intent.name  = 'SPM contrast';
        no.intent.code  = 1001;
end

create(no);
no.dat(:,:,:) = xSPM.fcn(ni.dat(:,:,:));
if gz
    gzip(nii);
    spm_unlink(nii);
end


%==========================================================================
% function make_ROI(fname,DIM,M,xY)
%==========================================================================
function make_ROI(fname,DIM,M,xY)
gz = strcmp(spm_file(fname,'ext'),'gz');
if gz, fname = spm_file(fname,'ext',''); end
R = struct(...
    'fname',  fname,...
    'dim',    DIM,...
    'dt',     [spm_type('uint8'), spm_platform('bigend')],...
    'mat',    M,...
    'pinfo',  [1,0,0]',...
    'descrip','ROI');
Q    = zeros(DIM);
[xY, XYZmm, j] = spm_ROI(xY, struct('dim',DIM,'mat',M));
Q(j) = 1;
R    = spm_write_vol(R,Q);
if gz
    gzip(R.fname);
    spm_unlink(R.fname);
end


%==========================================================================
% function id = coordspace(p,M,DIM,units,coordsys,idx)
%==========================================================================
function id = coordspace(p,M,DIM,units,coordsys,idx)
persistent index
if nargin == 6
    index = idx;
else
    if isempty(index)
        index = 1;
    else
        index = index + 1;
    end
end
% Convert from first voxel at [1,1,1] to first voxel at [0,0,0]
v2wm = M * [eye(4,3) [1 1 1 1]'];
M    = M(1:3,1:3);
id = getid(['niiri:coordinate_space_id_' num2str(index)],isHumanReadable);
p.entity(id,{...
    'prov:type',nidm_conv('nidm_CoordinateSpace',p),...
    'prov:label',{['Coordinate space ' num2str(index)],'xsd:string'},...
    nidm_conv('nidm_voxelToWorldMapping',p),{v2wm,'xsd:string'},...
    nidm_conv('nidm_voxelUnits',p),{units,'xsd:string'},...
    nidm_conv('nidm_voxelSize',p),{sqrt(diag(M'*M))','xsd:string'},...
    nidm_conv('nidm_inWorldCoordinateSystem',p),coordsys,...
    nidm_conv('nidm_numberOfDimensions',p),{numel(DIM),'xsd:int'},...
    nidm_conv('nidm_dimensionsInVoxels',p),{DIM,'xsd:string'}
    });

%==========================================================================
% function id = originalfile(p,file,idx,typ)
%==========================================================================
function id = originalfile(p,file,idx,typ)
id = getid([idx '_der'],isHumanReadable);
p.entity(id,{...
    'prov:type',typ,...
    'nfo:fileName',{spm_file(file,'filename'),'xsd:string'},...
    'dct:format',{'image/nifti','xsd:string'},...
    'crypto:sha512',{sha512sum(spm_file(file,'cpath')),'xsd:string'},...
    });

%==========================================================================
% function id = getid(id,humanReadable,checksum)
%==========================================================================
function id = getid(id,humanReadable,checksum)
if ~humanReadable
    if nargin == 2
        id = md5sum;
    else
        id = md5sum(checksum);
    end
    id = ['niiri:' id];
end

%==========================================================================
% function i = isHumanReadable(i)
%==========================================================================
function i = isHumanReadable(i)
persistent isHR
if nargin, isHR = i; end
if isempty(isHR), error('Default not set.'); end
i = isHR;

%==========================================================================
% function out = nidm_conv(in,p)
%==========================================================================
function out = nidm_conv(in,p)
persistent C
if isempty(C), C = nidm_constants; end

i = find(ismember(C(:,2),in));
if ~isempty(i)
    out = [C{i,2} ':'];
    if nargin == 2
        prefix = '';
        qname = C{i,1};
        j = find(qname == ':');
        if ~isempty(j)
            prefix = qname(1:j(end)-1);
            qname = qname(j(end)+1:end);
        end
        % should instead use ns = p.get_namespace;
        switch prefix
            case 'nidm'
                url = 'http://purl.org/nidash/nidm#';
            case 'spm'
                url = 'http://purl.org/nidash/spm#';
            case 'obo'
                url = 'http://purl.obolibrary.org/obo/';
            case 'nlx'
                url = 'http://uri.neuinfo.org/nif/nifstd/';
            case 'src'
                url = 'http://scicrunch.org/resolver/';
            otherwise
                warning('Unknown prefix "%s".',prefix);
                url = '';
        end
        p.add_namespace(C{i,2},[url qname]);
    end
else
    warning('Unknown element ''%s''.',in);
    out = in;
end

%==========================================================================
% function S = nidm_esc(S)
%==========================================================================
function S = nidm_esc(S)
S = strrep(S, sprintf('\n'), '\n');

%==========================================================================
% function nidm_store_constants
%==========================================================================
function nidm_store_constants
urlwrite('https://raw.githubusercontent.com/incf-nidash/nidm/master/nidm/nidm-results/terms/prefixes.csv','prefixes.csv');
C = reshape(textread('prefixes.csv','%s','delimiter',',','headerlines',1),2,[])';
fprintf('C = {...\n');
for i=1:size(C,1)
    fprintf('''%s'', ''%s'';...\n',C{i,1},strrep(C{i,2},':',''));
end
fprintf('};\n');

%==========================================================================
% function C = nidm_constants
%==========================================================================
function C = nidm_constants
% automatically generated by nidm_store_constants
C = {...
'obo:BFO_0000179', 'obo_BFOOWLspecificationlabel';...
'obo:BFO_0000180', 'obo_BFOCLIFspecificationlabel';...
'obo:IAO_0000002', 'obo_exampletobeeventuallyremoved';...
'obo:IAO_0000111', 'obo_editorpreferredterm';...
'obo:IAO_0000111', 'obo_editorpreferredterm';...
'obo:IAO_0000111', 'obo_editorpreferredterm';...
'obo:IAO_0000112', 'obo_exampleofusage';...
'obo:IAO_0000114', 'obo_hascurationstatus';...
'obo:IAO_0000115', 'obo_definition';...
'obo:IAO_0000115', 'obo_definition';...
'obo:IAO_0000115', 'obo_definition';...
'obo:IAO_0000115', 'obo_definition';...
'obo:IAO_0000116', 'obo_editornote';...
'obo:IAO_0000117', 'obo_termeditor';...
'obo:IAO_0000118', 'obo_alternativeterm';...
'obo:IAO_0000119', 'obo_definitionsource';...
'obo:IAO_0000120', 'obo_metadatacomplete';...
'obo:IAO_0000121', 'obo_organizationalterm';...
'obo:IAO_0000122', 'obo_readyforrelease';...
'obo:IAO_0000123', 'obo_metadataincomplete';...
'obo:IAO_0000124', 'obo_uncurated';...
'obo:IAO_0000125', 'obo_pendingfinalvetting';...
'obo:IAO_0000136', 'obo_isabout';...
'obo:IAO_0000232', 'obo_curatornote';...
'obo:IAO_0000412', 'obo_importedfrom';...
'obo:IAO_0000423', 'obo_tobereplacedwithexternalontologyterm';...
'obo:IAO_0000428', 'obo_requiresdiscussion';...
'obo:IAO_0000600', 'obo_elucidation';...
'obo:OBI_0000251', 'obo_cluster';...
'obo:OBI_0001265', 'obo_FWERadjustedpvalue';...
'obo:OBI_0001442', 'obo_qvalue';...
'obo:STATO_0000030', 'obo_ChiSquaredstatistic';...
'obo:STATO_0000039', 'obo_statistic';...
'obo:STATO_0000051', 'obo_Poissondistribution';...
'obo:STATO_0000067', 'obo_continuousprobabilitydistribution';...
'obo:STATO_0000117', 'obo_discreteprobabilitydistribution';...
'obo:STATO_0000119', 'obo_modelparameterestimation';...
'obo:STATO_0000129', 'obo_hasvalue';...
'obo:STATO_0000176', 'obo_tstatistic';...
'obo:STATO_0000193', 'obo_studygrouppopulation';...
'obo:STATO_0000225', 'obo_probabilitydistribution';...
'obo:STATO_0000227', 'obo_normaldistribution';...
'obo:STATO_0000276', 'obo_binomialdistribution';...
'obo:STATO_0000282', 'obo_Fstatistic';...
'obo:STATO_0000323', 'obo_contrastweightmatrix';...
'obo:STATO_0000346', 'obo_covariancestructure';...
'obo:STATO_0000357', 'obo_Toeplitzcovariancestructure';...
'obo:STATO_0000362', 'obo_compoundsymmetrycovariancestructure';...
'obo:STATO_0000370', 'obo_ordinaryleastsquaresestimation';...
'obo:STATO_0000371', 'obo_weightedleastsquaresestimation';...
'obo:STATO_0000372', 'obo_generalizedleastsquaresestimation';...
'obo:STATO_0000373', 'obo_iterativelyreweightedleastsquaresestimation';...
'obo:STATO_0000374', 'obo_feasiblegeneralizedleastsquaresestimation';...
'obo:STATO_0000376', 'obo_Zstatistic';...
'obo:STATO_0000405', 'obo_unstructuredcovariancestructure';...
'iao:iao.owl', 'iao_IAORelease20150223';...
'dc:contributor', 'dc_Contributor';...
'dc:creator', 'dc_Creator';...
'dc:date', 'dc_Date';...
'dc:date', 'dc_Date';...
'dc:description', 'dc_Description';...
'dc:title', 'dc_Title';...
'fsl:FSL_0000001', 'fsl_FSLsGammaDifferenceHRF';...
'fsl:FSL_0000002', 'fsl_GaussianRunningLineDriftModel';...
'fsl:FSL_0000003', 'fsl_FSLsTemporalDerivative';...
'fsl:FSL_0000004', 'fsl_driftCutoffPeriod';...
'fsl:FSL_0000005', 'fsl_featVersion';...
'nidm:NIDM_0000001', 'nidm_ContrastEstimation';...
'nidm:NIDM_0000002', 'nidm_ContrastMap';...
'nidm:NIDM_0000004', 'nidm_BinaryMap';...
'nidm:NIDM_0000007', 'nidm_ClusterDefinitionCriteria';...
'nidm:NIDM_0000008', 'nidm_ClusterLabelsMap';...
'nidm:NIDM_0000009', 'nidm_Colin27CoordinateSystem';...
'nidm:NIDM_0000011', 'nidm_ConjunctionInference';...
'nidm:NIDM_0000012', 'nidm_ConnectivityCriterion';...
'nidm:NIDM_0000013', 'nidm_ContrastStandardErrorMap';...
'nidm:NIDM_0000015', 'nidm_Coordinate';...
'nidm:NIDM_0000016', 'nidm_CoordinateSpace';...
'nidm:NIDM_0000017', 'nidm_CustomCoordinateSystem';...
'nidm:NIDM_0000019', 'nidm_DesignMatrix';...
'nidm:NIDM_0000020', 'nidm_DisplayMaskMap';...
'nidm:NIDM_0000021', 'nidm_regressorNames';...
'nidm:NIDM_0000023', 'nidm_ErrorModel';...
'nidm:NIDM_0000024', 'nidm_ExchangeableError';...
'nidm:NIDM_0000025', 'nidm_ExcursionSetMap';...
'nidm:NIDM_0000026', 'nidm_ExtentThreshold';...
'nidm:NIDM_0000027', 'nidm_NIDMResults';...
'nidm:NIDM_0000028', 'nidm_FiniteImpulseResponseBasisSet';...
'nidm:NIDM_0000029', 'nidm_GammaDifferenceHRF';...
'nidm:NIDM_0000030', 'nidm_GammaBasisSet';...
'nidm:NIDM_0000031', 'nidm_GammaHRF';...
'nidm:NIDM_0000033', 'nidm_GrandMeanMap';...
'nidm:NIDM_0000034', 'nidm_HeightThreshold';...
'nidm:NIDM_0000035', 'nidm_HemodynamicResponseFunction';...
'nidm:NIDM_0000036', 'nidm_ConvolutionBasisSet';...
'nidm:NIDM_0000037', 'nidm_HemodynamicResponseFunctionDerivative';...
'nidm:NIDM_0000038', 'nidm_Icbm452AirCoordinateSystem';...
'nidm:NIDM_0000039', 'nidm_Icbm452Warp5CoordinateSystem';...
'nidm:NIDM_0000040', 'nidm_IcbmMni152LinearCoordinateSystem';...
'nidm:NIDM_0000041', 'nidm_IcbmMni152NonLinear2009aAsymmetricCoordinateSystem';...
'nidm:NIDM_0000042', 'nidm_IcbmMni152NonLinear2009aSymmetricCoordinateSystem';...
'nidm:NIDM_0000043', 'nidm_IcbmMni152NonLinear2009bAsymmetricCoordinateSystem';...
'nidm:NIDM_0000044', 'nidm_IcbmMni152NonLinear2009bSymmetricCoordinateSystem';...
'nidm:NIDM_0000045', 'nidm_IcbmMni152NonLinear2009cAsymmetricCoordinateSystem';...
'nidm:NIDM_0000046', 'nidm_IcbmMni152NonLinear2009cSymmetricCoordinateSystem';...
'nidm:NIDM_0000047', 'nidm_IcbmMni152NonLinear6thGenerationCoordinateSystem';...
'nidm:NIDM_0000048', 'nidm_IndependentError';...
'nidm:NIDM_0000049', 'nidm_Inference';...
'nidm:NIDM_0000050', 'nidm_Ixi549CoordinateSystem';...
'nidm:NIDM_0000051', 'nidm_MNICoordinateSystem';...
'nidm:NIDM_0000052', 'nidm_Map';...
'nidm:NIDM_0000053', 'nidm_MapHeader';...
'nidm:NIDM_0000054', 'nidm_MaskMap';...
'nidm:NIDM_0000055', 'nidm_Mni305CoordinateSystem';...
'nidm:NIDM_0000056', 'nidm_ModelParametersEstimation';...
'nidm:NIDM_0000057', 'nidm_NIDMObjectModel';...
'nidm:NIDM_0000059', 'nidm_NonParametricSymmetricDistribution';...
'nidm:NIDM_0000060', 'nidm_OneTailedTest';...
'nidm:NIDM_0000061', 'nidm_ParameterEstimateMap';...
'nidm:NIDM_0000062', 'nidm_Peak';...
'nidm:NIDM_0000063', 'nidm_PeakDefinitionCriteria';...
'nidm:NIDM_0000064', 'nidm_PixelConnectivityCriterion';...
'nidm:NIDM_0000066', 'nidm_ResidualMeanSquaresMap';...
'nidm:NIDM_0000067', 'nidm_CustomBasisSet';...
'nidm:NIDM_0000068', 'nidm_SearchSpaceMaskMap';...
'nidm:NIDM_0000069', 'nidm_FourierBasisSet';...
'nidm:NIDM_0000070', 'nidm_SupraThresholdCluster';...
'nidm:NIDM_0000071', 'nidm_ErrorParameterMapWiseDependence';...
'nidm:NIDM_0000072', 'nidm_ConstantParameter';...
'nidm:NIDM_0000073', 'nidm_IndependentParameter';...
'nidm:NIDM_0000074', 'nidm_RegularizedParameter';...
'nidm:NIDM_0000075', 'nidm_StandardizedCoordinateSystem';...
'nidm:NIDM_0000076', 'nidm_StatisticMap';...
'nidm:NIDM_0000077', 'nidm_SubjectCoordinateSystem';...
'nidm:NIDM_0000078', 'nidm_TalairachCoordinateSystem';...
'nidm:NIDM_0000079', 'nidm_TwoTailedTest';...
'nidm:NIDM_0000080', 'nidm_VoxelConnectivityCriterion';...
'nidm:NIDM_0000081', 'nidm_WorldCoordinateSystem';...
'nidm:NIDM_0000082', 'nidm_clusterLabelId';...
'nidm:NIDM_0000083', 'nidm_clusterSizeInVertices';...
'nidm:NIDM_0000084', 'nidm_clusterSizeInVoxels';...
'nidm:NIDM_0000085', 'nidm_contrastName';...
'nidm:NIDM_0000086', 'nidm_coordinateVector';...
'nidm:NIDM_0000087', 'nidm_DriftModel';...
'nidm:NIDM_0000088', 'nidm_hasDriftModel';...
'nidm:NIDM_0000089', 'nidm_dependenceMapWiseDependence';...
'nidm:NIDM_0000090', 'nidm_dimensionsInVoxels';...
'nidm:NIDM_0000091', 'nidm_effectDegreesOfFreedom';...
'nidm:NIDM_0000092', 'nidm_equivalentZStatistic';...
'nidm:NIDM_0000093', 'nidm_errorDegreesOfFreedom';...
'nidm:NIDM_0000094', 'nidm_errorVarianceHomogeneous';...
'nidm:NIDM_0000096', 'nidm_grandMeanScaling';...
'nidm:NIDM_0000097', 'nidm_hasAlternativeHypothesis';...
'nidm:NIDM_0000098', 'nidm_hasClusterLabelsMap';...
'nidm:NIDM_0000099', 'nidm_hasConnectivityCriterion';...
'nidm:NIDM_0000100', 'nidm_hasErrorDependence';...
'nidm:NIDM_0000101', 'nidm_hasErrorDistribution';...
'nidm:NIDM_0000102', 'nidm_hasHRFBasis';...
'nidm:NIDM_0000103', 'nidm_hasMapHeader';...
'nidm:NIDM_0000104', 'nidm_inCoordinateSpace';...
'nidm:NIDM_0000105', 'nidm_inWorldCoordinateSystem';...
'nidm:NIDM_0000106', 'nidm_isUserDefined';...
'nidm:NIDM_0000107', 'nidm_maskedMedian';...
'nidm:NIDM_0000108', 'nidm_maxNumberOfPeaksPerCluster';...
'nidm:NIDM_0000109', 'nidm_minDistanceBetweenPeaks';...
'nidm:NIDM_0000110', 'nidm_GaussianHRF';...
'nidm:NIDM_0000111', 'nidm_numberOfSupraThresholdClusters';...
'nidm:NIDM_0000112', 'nidm_numberOfDimensions';...
'nidm:NIDM_0000113', 'nidm_objectModel';...
'nidm:NIDM_0000114', 'nidm_pValue';...
'nidm:NIDM_0000115', 'nidm_pValueFWER';...
'nidm:NIDM_0000116', 'nidm_pValueUncorrected';...
'nidm:NIDM_0000117', 'nidm_pixel4connected';...
'nidm:NIDM_0000118', 'nidm_pixel8connected';...
'nidm:NIDM_0000119', 'nidm_qValueFDR';...
'nidm:NIDM_0000120', 'nidm_randomFieldStationarity';...
'nidm:NIDM_0000121', 'nidm_searchVolumeInVoxels';...
'nidm:NIDM_0000122', 'nidm_softwareVersion';...
'nidm:NIDM_0000123', 'nidm_statisticType';...
'nidm:NIDM_0000124', 'nidm_targetIntensity';...
'nidm:NIDM_0000126', 'nidm_varianceMapWiseDependence';...
'nidm:NIDM_0000127', 'nidm_version';...
'nidm:NIDM_0000128', 'nidm_voxel18connected';...
'nidm:NIDM_0000129', 'nidm_voxel26connected';...
'nidm:NIDM_0000130', 'nidm_voxel6connected';...
'nidm:NIDM_0000131', 'nidm_voxelSize';...
'nidm:NIDM_0000132', 'nidm_voxelToWorldMapping';...
'nidm:NIDM_0000133', 'nidm_voxelUnits';...
'nidm:NIDM_0000134', 'nidm_withEstimationMethod';...
'nidm:NIDM_0000135', 'nidm_ContrastVarianceMap';...
'nidm:NIDM_0000136', 'nidm_searchVolumeInUnits';...
'nidm:NIDM_0000137', 'nidm_searchVolumeInVertices';...
'nidm:NIDM_0000138', 'nidm_hasMaximumIntensityProjection';...
'nidm:NIDM_0000139', 'nidm_coordinateVectorInVoxels';...
'nidm:NIDM_0000140', 'nidm_ClusterCenterOfGravity';...
'nidm:NIDM_0000141', 'nidm_expectedNumberOfClusters';...
'nidm:NIDM_0000142', 'nidm_expectedNumberOfVerticesPerCluster';...
'nidm:NIDM_0000143', 'nidm_expectedNumberOfVoxelsPerCluster';...
'nidm:NIDM_0000144', 'nidm_ReselsPerVoxelMap';...
'nidm:NIDM_0000145', 'nidm_noiseRoughnessInVoxels';...
'nidm:NIDM_0000146', 'nidm_heightCriticalThresholdFDR05';...
'nidm:NIDM_0000147', 'nidm_heightCriticalThresholdFWE05';...
'nidm:NIDM_0000148', 'nidm_reselSizeInVoxels';...
'nidm:NIDM_0000149', 'nidm_searchVolumeInResels';...
'nidm:NIDM_0000150', 'nidm_LinearSplineBasisSet';...
'nidm:NIDM_0000151', 'nidm_SineBasisSet';...
'nidm:NIDM_0000156', 'nidm_clusterSizeInResels';...
'nidm:NIDM_0000157', 'nidm_noiseFWHMInUnits';...
'nidm:NIDM_0000158', 'nidm_noiseFWHMInVertices';...
'nidm:NIDM_0000159', 'nidm_noiseFWHMInVoxels';...
'nidm:NIDM_0000160', 'nidm_PValueUncorrected';...
'nidm:NIDM_0000161', 'nidm_equivalentThreshold';...
'nidm:NIDM_0000162', 'nidm_Threshold';...
'nidm:NIDM_0000163', 'nidm_ContrastExplainedMeanSquareMap';...
'nidm:NIDM_0000164', 'nidm_NeuroimagingAnalysisSoftware';...
'nidm:NIDM_0000165', 'nidm_NIDMResultsExporter';...
'nidm:NIDM_0000166', 'nidm_NIDMResultsExport';...
'nidm:NIDM_0000167', 'nidm_nidmfsl';...
'nidm:NIDM_0000168', 'nidm_spm_results_nidm';...
'nidm:NIDM_0000169', 'nidm_Data';...
'nidm:NIDM_0000170', 'nidm_groupName';...
'nidm:NIDM_0000171', 'nidm_numberOfSubjects';...
'nidm:NIDM_0000172', 'nidm_hasMRIProtocol';...
'spm:SPM_0000001', 'spm_SPMsDriftCutoffPeriod';...
'spm:SPM_0000002', 'spm_DCTDriftModel';...
'spm:SPM_0000003', 'spm_SPMsDispersionDerivative';...
'spm:SPM_0000004', 'spm_SPMsCanonicalHRF';...
'spm:SPM_0000005', 'spm_PartialConjunctionInference';...
'spm:SPM_0000006', 'spm_SPMsTemporalDerivative';...
'spm:SPM_0000010', 'spm_searchVolumeReselsGeometry';...
'spm:SPM_0000011', 'spm_smallestSignificantClusterSizeInVerticesFDR05';...
'spm:SPM_0000012', 'spm_smallestSignificantClusterSizeInVerticesFWE05';...
'spm:SPM_0000013', 'spm_smallestSignificantClusterSizeInVoxelsFDR05';...
'spm:SPM_0000014', 'spm_smallestSignificantClusterSizeInVoxelsFWE05';...
'spm:SPM_0000015', 'spm_partialConjunctionDegree';...
'prv:PropertyReification', 'prv_PropertyReification';...
'prv:object_property', 'prv_hasobjectproperty';...
'prv:reification_class', 'prv_hasreificationclass';...
'prv:shortcut', 'prv_hasshortcut';...
'prv:shortcut_property', 'prv_hasshortcutproperty';...
'prv:subject_property', 'prv_hassubjectproperty';...
'src:SCR_002823', 'src_FSL';...
'src:SCR_007037', 'src_SPM';...
'nlx:birnlex_2094', 'nlx_Imaginginstrument';...
'nlx:birnlex_2100', 'nlx_Magneticresonanceimagingscanner';...
'nlx:birnlex_2177', 'nlx_MRIprotocol';...
'nlx:birnlex_2250', 'nlx_FunctionalMRIprotocol';...
'nlx:birnlex_2251', 'nlx_StructuralMRIprotocol';...
'nlx:ixl_0050000', 'nlx_Positronemissiontomographyscanner';...
'nlx:ixl_0050001', 'nlx_Singlephotonemissioncomputedtomographyscanner';...
'nlx:ixl_0050002', 'nlx_Magnetoencephalographymachine';...
'nlx:ixl_0050003', 'nlx_Electroencephalographymachine';...
'nlx:ixl_0050004', 'nlx_AnatomicalMRIprotocol';...
'nlx:nlx_inv_20090249', 'nlx_Diffusionweightedimagingprotocol';...
};

%==========================================================================
% Upload to NeuroVault
%==========================================================================
function upload_to_neurovault(nidmfile,token)

neurovault = 'http://neurovault.org';

% Get token
if nargin > 1
    addpref('neurovault','token',token);
elseif ispref('neurovault','token')
    token = getpref('neurovault','token');
else
    if spm('CmdLine')
        error('Upload to NeuroVault requires a token-based authentication.');
    end
    fprintf([...
        'To set up integration with NeuroVault you need to obtain a\n',...
        'secret token first. Please go to http://neurovault.org/ and\n',...
        'generate a new token. To complete the procedure you will need\n',...
        'to log into NeuroVault. If you don''t have an account you can\n',...
        'create one or simply log in using Facebook or Google. When you\n',...
        'generated a token (a string of 40 random characters), copy it\n',...
        'into the dialog box. If you ever lose access to this machine\n',...
        'you can delete the token on NeuroVault.org thus preventing\n',...
        'unauthorized parties to access your NeuroVault data\n']);
    token = inputdlg('Token','Enter NeuroVault token',1,{''},'on');
    if isempty(token) || isempty(token{1}), return; else token = char(token); end
    addpref('neurovault','token',token);
end
auth = ['Bearer ' token];

% Check token
url = [neurovault '/api/my_collections/'];
statusCode = http_request('head', url, 'Authorization', auth);
if isnan(statusCode)
    error('Cannot connect to NeuroVault.');
elseif statusCode ~= 200
    warning('Failed authentication with NeuroVault: invalid token.');
    rmpref('neurovault','token');
    upload_to_neurovault(nidmfile);
    return;
end

% Get user name
url = [neurovault '/api/user/'];
[statusCode, responseBody] = http_request('jsonGet', url, 'Authorization', auth);
if statusCode == 200
    responseBody = spm_jsonread(responseBody);
    owner_name = responseBody.username;
else
    err = spm_jsonread(responseBody);
    error(char(err.detail));
end
my_collection = sprintf('%s''s %s Collection',owner_name,spm('Ver'));

% Get the list of collections
url   = [neurovault '/api/my_collections/'];
[statusCode, responseBody] = http_request('jsonGet', url, 'Authorization', auth);
if statusCode == 200
    collections = spm_jsonread(responseBody);
else
    err = spm_jsonread(responseBody);
    error(char(err.detail));
end

% Create a new collection if needed
id = NaN;
for i=1:numel(collections.results)
    if strcmp(collections.results(i).name,my_collection)
        id = collections.results(i).id;
        break;
    end
end
if isnan(id)
    url   = [neurovault '/api/collections/'];
    requestParts = [];
    requestParts.Type = 'string';
    requestParts.Name = 'name';
    requestParts.Body = my_collection;
    [statusCode, responseBody] = http_request('multipartPost', url, requestParts, 'Authorization', auth);
    if statusCode == 201
        collections = spm_jsonread(responseBody);
        id = collections.id;
    else
        err = spm_jsonread(responseBody);
        error(char(err.name));
    end
end

% Upload NIDM-Results in NeuroVault's collection
url = [neurovault sprintf('/api/collections/%d/nidm_results/',id)];
requestParts = [];
requestParts(1).Type = 'string';
requestParts(1).Name = 'name';
requestParts(1).Body = 'No name'; % NAME
requestParts(2).Type = 'string';
requestParts(2).Name = 'description';
requestParts(2).Body = 'No description'; % DESCRIPTION
requestParts(3).Type = 'file';
requestParts(3).Name = 'zip_file';
requestParts(3).Body = nidmfile;
[statusCode, responseBody] = http_request('multipartPost', url, requestParts, 'Authorization', auth);
if statusCode == 201
    responseBody = spm_jsonread(responseBody);
    cmd = 'web(''%s'',''-browser'')';
    fprintf('Uploaded to %s\n',spm_file(responseBody.url,'link',cmd));
else
    err = spm_jsonread(responseBody);
    error(char(err.name));
end

%==========================================================================
% HTTP requests with missing-http
%==========================================================================
function varargout = http_request(action, url, varargin)
% Use missing-http from Paul Sexton:
% https://github.com/psexton/missing-http/

persistent jar_added_to_path
if isempty(jar_added_to_path)
    try
        net.psexton.missinghttp.MatlabShim;
    catch
        err = lasterror;
        if strcmp(err.identifier, 'MATLAB:dispatcher:noMatchingConstructor')
            jar_added_to_path = true;
        else
            jarfile = fullfile(spm('Dir'),'external','missing-http','missing-http.jar');
            if spm_existfile(jarfile)
                javaaddpath(jarfile); % this clears global
                jar_added_to_path = true;
            else
                error('HTTP library missing.');
            end
        end
    end
end

switch lower(action)
    case 'head'
        try
            response = char(net.psexton.missinghttp.MatlabShim.head(url, varargin));
        catch
            response = 'NaN';
        end
        statusCode   = str2double(response);
        varargout    = { statusCode };
    case 'jsonget'
        try
            response = cell(net.psexton.missinghttp.MatlabShim.jsonGet(url, varargin));
        catch
            response = {'NaN',''};
        end
        statusCode   = str2double(response{1});
        responseBody = response{2};
        varargout    = { statusCode, responseBody };
    case 'multipartpost'
        requestParts = varargin{1};
        crp = cell(1, numel(requestParts));
        for k=1:numel(requestParts)
            crp{k}   = sprintf('%s\n%s\n%s', requestParts(k).Type, requestParts(k).Name, requestParts(k).Body);
        end
        try
            response = cell(net.psexton.missinghttp.MatlabShim.multipartPost(url, crp, varargin(2:end)));
        catch
            response = {'NaN',''};
        end
        statusCode   = str2double(response{1});
        responseBody = response{2};
        varargout    = { statusCode, responseBody };
    otherwise
        error('Unknown HTTP request.');
end

%==========================================================================
% HTTP requests with MATLAB R2016b
%==========================================================================
function varargout = http_request2(action, url, varargin)

switch lower(action)
    case 'head'
        opt = weboptions('HeaderFields',varargin,'ContentType','text');
        try
            webread(url,opt); % it's a GET, not a HEAD
            statusCode = 200;
        catch
            err = lasterror;
            s = regexp(err.identifier,'MATLAB:webservices:HTTP(\d+)','tokens');
            if isempty(s)
                statusCode = NaN;
            else
                statusCode = str2double(s{1}{1});
            end
        end
        varargout    = { statusCode };
    case 'jsonget'
        opt = weboptions('HeaderFields',varargin,'ContentType','text');
        try
            responseBody = webread(url,opt);
            statusCode = 200;
        catch
            err = lasterror;
            s = regexp(err.identifier,'MATLAB:webservices:HTTP(\d+)','tokens');
            if isempty(s)
                statusCode = NaN;
            else
                statusCode = str2double(s{1}{1});
            end
        end
        varargout    = { statusCode, responseBody };
    case 'multipartpost'
        requestParts = varargin{1};
        opt = weboptions('HeaderFields',varargin(2:end),'ContentType','text');
        try
            responseBody = webread(url,opt);
            statusCode = 200;
        catch
            err = lasterror;
            s = regexp(err.identifier,'MATLAB:webservices:HTTP(\d+)','tokens');
            if isempty(s)
                statusCode = NaN;
            else
                statusCode = str2double(s{1}{1});
            end
        end
        varargout    = { statusCode, responseBody };
    otherwise
        error('Unknown HTTP request.');
end

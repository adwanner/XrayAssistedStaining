basepath='/home/awanner/Documents/PostDoc/Manuscripts/2020 x-ray assisted staining/codeExample/ExampleData/';
in_filename=[basepath filesep 'IAC20220210_B20220208_20190118.txrm'];

TIFpath=[basepath filesep 'TIF' filesep ];
alignedTIFpath=[TIFpath '/aligned/'];
boundarypath=[alignedTIFpath filesep '/BoundaryCoords/'];

metadatafile=[basepath filesep 'metadata.mat'];
displacementfile=[basepath filesep 'pairwise_displacements.mat'];
densityfile=[basepath filesep 'DensityFile.mat'];
modelfile=[basepath filesep 'ModelFile.mat'];

%% Step01: Extract TXRM file
dooverwrite=1;
doextract_xrm=1;
read_metadata=1;
calc_offsets=1;
writeAlignedImage=1;
Step01_ExtractTXRM;

%% Step02: Detect sample surface
doRadialProjectionCorrection=true;
dooverwrite=true;
% doplot=1;
% dosave=0;
doplot=0;
dosave=1;
%the pixel range for detecting the surface of the sample might has to be
%adapted to the specific sample:
YStart=1; 
YEnd=800;
Step02_DetectSampleSurface;

%% Step03: Measure diffusion
dosave=1;
calc_diffusion=1; 
dooverwrite=1;
doExpansionCorrection=1;
doRadialProjectionCorrection=1;
Step03_MeasureDiffusion;

%% Step04: Model diffusion
dosavedata=0;
dooverwrite=1;

domodel=1;
doplotprogress=0;

dobaselinecorr=1;
donorm=1;
dogrowth=1;
domasking=1;
Step04_ModelDiffusion;

%% Step05: Plot diffusion data and model
dobaselinecorr=true;
donorm=true;
dounmask=1;
dogrowth=1;

dodetectpeaks=0;

Step05_PlotDiffusionAndModel;
function y = proj_forw(G, Gopt, x)
%--------------------------------------------------------------------------
% forward projection 
%
% Guobao Wang @ UC Davis (10-01-2012)
%

if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end
x = x(:);

switch Gopt.mtype
    case 'matlab'
        y = G * x;
        
    case 'handle'
        y = Gopt.forw(x);
        
    case 'fessler'
        y = G * x(Gopt.ig.mask,:);
        y = double(y);
        
    case 'usc'
        put_data('temp.image.usc', x, 'float32');
        setenv('NUMOFWORKERS', '32');
        unix(sprintf('%s %s temp.image.usc', Gopt.pforw, Gopt.cfile));
        y = get_data('temp.image.usc.forw');
        
    case 'zhou'
        put_data('temp.image', x, 'float32');
        setenv('OMP_NUM_THREADS', '32');
        unix(sprintf('%s -c %s -f temp.image -o temp.image.forw', Gopt.pforw, Gopt.cfile));
        y = get_data('temp.image.forw');
        
    case 'geprt'
        acqParams.nU                   = 249;
        acqParams.nV                   = 553;
        acqParams.nPhi                 = 210;
        acqParams.sU                   = 3.1966;
        acqParams.sV                   = 3.2646;
        reconParams.nx                 = 214;
        acqParams.nZ                   = 47;
        reconParams.FOV                = 700;
        acqParams.sV                   = 3.2646;
        scanner.numBlocksPerRing       = 70;
        scanner.radialCrystalsPerBlock = 6;
        scanner.radBlockSize           = 38.35;
        subsetList                     = 0;
        reconParams.numSubsets         = 1;
        reconParams.xOffset            = -0.17089;
        reconParams.yOffset            = 1.6078;
        reconParams.rotate             = -4.6979;
        scanner.effectiveRingDiameter  = 903;
        x = reshape(x,reconParams.nx, reconParams.nx, acqParams.nZ);
        y = FDD3D_ng(single(x), acqParams.nU, acqParams.nV, acqParams.nPhi, acqParams.sU, acqParams.sV, ...
                     reconParams.nx, acqParams.nZ, reconParams.FOV/reconParams.nx, acqParams.sV, ...
                     scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, scanner.radBlockSize, ...
                     subsetList, reconParams.numSubsets, reconParams.xOffset, reconParams.yOffset, ...
                     reconParams.rotate, scanner.effectiveRingDiameter);
        y = double(y(:));
       
end

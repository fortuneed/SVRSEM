function x = proj_back(G, Gopt, y)
%--------------------------------------------------------------------------
% back projection 
%
% Guobao Wang @ UC Davis (10-01-2012)
%

if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end
y = y(:);

switch Gopt.mtype
    case 'matlab'
        x = G' * y;
        
    case 'handle'
        x = Gopt.back(y);
        
    case 'fessler'
        a = double( G' * y ); 
        for m = 1:size(a,2)
            temp = Gopt.ig.embed(a(:,m));
            x(:,m) = temp(:);
        end
        
    case 'usc'
        put_data('temp.proj.usc', y, 'float32');
        setenv('NUMOFWORKERS', '32');
        unix(sprintf('%s %s temp.proj.usc', Gopt.pback, Gopt.cfile));
        x = get_data('temp.proj.usc.back');
        
    case 'zhou'
        put_data('temp.proj', y, 'float32');
        setenv('OMP_NUM_THREADS', '32');
        unix(sprintf('%s -c %s -b temp.proj -o temp.proj.back', Gopt.pback, Gopt.cfile));
        x = get_data('temp.proj.back');
    
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
        y = reshape(y, acqParams.nU, acqParams.nV, acqParams.nPhi);
        x = BDD3D_ng(single(y), acqParams.nU, acqParams.nV, acqParams.nPhi, acqParams.sU, acqParams.sV, ...
                     reconParams.nx, acqParams.nZ, reconParams.FOV/reconParams.nx, acqParams.sV, ...
                     scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, scanner.radBlockSize, ...
                     subsetList, reconParams.numSubsets, reconParams.xOffset, reconParams.yOffset, ...
                     reconParams.rotate, scanner.effectiveRingDiameter);
        x = double(x(:));
end   

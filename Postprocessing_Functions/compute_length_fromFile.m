function [dist_hist, time, Ln, Lw ] =  compute_length_fromFile(numFramesFileName, positionsFileName, dt, write_freq, minL, maxL  )

numFramesFile=fopen(numFramesFileName);
positionsFile=fopen(positionsFileName);
%read number of frames
numFrames=fscanf(numFramesFile,'%g',1);
numBins = (maxL-minL)/minL;
lengths = minL:(maxL-minL)/numBins :maxL;

time = (1:numFrames)*dt*write_freq;
dist_hist =zeros(numFrames,length(lengths));
Ln = zeros(1,numFrames);
Lw = zeros(1,numFrames);
for frame= 1 : numFrames
    segment = 0;
    % read the number of fibers in this frame
    numFibers=fscanf(positionsFile,'%g',1);
    for  fiber = 1: numFibers
        % read number of hinges
        numHinges=fscanf(positionsFile,'%g',1);
        length_ =0;
        %load x1, y1, z1
        pts(1:3)=fscanf(positionsFile,'%g',3);
        x1 = pts(1); %
        y1 = pts(2);
        z1 = pts(3);
        
        for hinge = 2 : numHinges
            
            %load x2, y2, z2
            pts(1:3)=fscanf(positionsFile,'%g',3);
            x2 = pts(1); %
            y2 = pts(2);
            z2 = pts(3);
            
            % find the components of the segments
            % do components have to be positive? I dont know
            p1 = (x2-x1);
            p2 = (y2-y1);
            p3 = (z2-z1);
            
            % find segment length and add it to the fiber length
            
            length_ = length_ + sqrt(p1*p1+p2*p2+p3*p3);
            
            % normalize components
                      
            x1 = x2;
            y1 = y2;
            z1 = z2;
            
            segment = segment +1;
            %
        end
        
        binIndex = floor( (length_ -minL/2)/minL   ) +1;
        
        dist_hist(frame, binIndex) = dist(frame, binIndex) +1;
    
    end
    
    Ln(frame) = sum( dist_hist(frame, :).*lengths )/sum(dist_hist(frame, :));
    Lw(frame) = sum( dist_hist(frame, :).*lengths.*lengths )/sum(dist_hist(frame, :).*lengths);
    
    
    if ( mod(frame,numFrames/10)==0 )
        fprintf('Loading Mechanistic model results %d percent completed \n',frame/numFrames*100 );
        
    end
    
    
end

fclose(numFramesFile);
fclose(positionsFile);


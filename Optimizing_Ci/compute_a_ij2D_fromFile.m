function [a_ij, totalDeformation ] =  compute_a_ij2D_fromFile(numFramesFileName, positionsFileName, dt, shear_rate, write_freq  )
% this I will leave as manual inputs.

numFramesFile=fopen(numFramesFileName)
positionsFile=fopen(positionsFileName)
%read number of frames
numFrames=fscanf(numFramesFile,'%g',1);


totalDeformation = (1:numFrames)*dt*shear_rate*write_freq;

for frame= 1 : numFrames
    segment = 0;
    % read the number of fibers in this frame
    numFibers=fscanf(positionsFile,'%g',1);
    for  fiber = 1: numFibers
        % read number of hinges
        numHinges=fscanf(positionsFile,'%g',1);
        
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
            %p3 = (z2-z1);
            
            % find segment length
            
            %length = sqrt(p1*p1+p2*p2+p3*p3);
            length = sqrt(p1*p1+p2*p2);
            
            % normalize components
            p1 = p1/ length;
            p2 = p2/ length;
            %p3 = p3/ length;
            
            segment = segment + 1;
            p1_array(segment) = p1;
            p2_array(segment) = p2;
            %p3_array(segment) = p3;
            
            x1 = x2;
            y1 = y2;
            %z1 = z2;
            
            %
        end
    end
    
    
    a11 = dot(p1_array,p1_array)/segment; % a11 is the average value of the product p1*p1
    a12 = dot(p1_array,p2_array)/segment; % the others are analogous
    %a13 = dot(p1_array,p3_array)/segment;
    
    a22 = dot(p2_array,p2_array)/segment;
    %a23 = dot(p2_array,p3_array)/segment;
    %a33 = dot(p3_array,p3_array)/segment;
    
    a_ij(frame,:) = [a11 a12 a22 ];
    
    if ( mod(frame,numFrames/10)==0 )
        fprintf('Loading Mechanistic model results %d percent completed \n',frame/numFrames*100 );
        
    end
    
    
end


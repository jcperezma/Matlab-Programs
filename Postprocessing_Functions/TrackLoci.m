%transforms cluster into single segment fibers
clear all
close all
clc

dt =1e-6;
write_freq = 1000;
positionsFile=fopen('positions.out','r');
framesFile=fopen('nbr_frames.txt','r');

numFrames=fscanf(framesFile,'%g ',1);
%read number of frames
pts=zeros(numFrames,3);
for frame =1 : numFrames
    % read the number of fibers in this frame
    numFibers=fscanf(positionsFile,'%g',1);
    
    for  fiber = 1: numFibers
        % read number of hinges
        numHinges=fscanf(positionsFile,'%g ',1);
        
        pts(frame,1:3)=fscanf(positionsFile,'%g',3);
        
        % ignore rest of the hinges
        for hinge = 2 : numHinges 
           % Skip lines
            %load x2, y2, z2
            fscanf(positionsFile,'%g',3);
            %
        end
    
    end
end 
    fclose(positionsFile);
    fclose(framesFile);
 t = (1:numFrames)*dt*write_freq;
plot(t,pts(:,2),'x')



writeVTKLine( pts(:,1),pts(:,2),pts(:,3), 'orbit.vtp' );



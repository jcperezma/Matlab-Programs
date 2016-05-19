% computes the length distribution of a simulation
folder_name = uigetdir; % select the folder that contains the output folder
numFramesFileName=[folder_name '\output\nbr_frames.txt'];
positionsFileName=[folder_name '\output\positions.out'];

% get this from the fibers.in file. 
dt = 1e-5;
write_freq =2000;

% define the minimum and maximum length in the simulation
minL = 0.25;
maxL = 12.5;


[dist_hist, time, Ln, Lw ] =  compute_length_fromFile(numFramesFileName, positionsFileName, dt, write_freq, minL, maxL  );

%Plot the average length by number
figure
plot(time,Ln*1000 )
ylabel('Average Length by number [mm]')
xlabel('Time[s]')


%Plot the average lenggth by length
figure
plot(time,Lw*1000)
ylabel('Average Length by number [mm]')
xlabel('Time[s]')
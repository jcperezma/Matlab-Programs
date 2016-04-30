% computes the 2D orientation tensor of a simulation
folder_name = uigetdir; % select the folder that contains the output folder
numFramesFileName=[folder_name '\output\nbr_frames.txt'];
positionsFileName=[folder_name '\output\positions.out'];

dt = 1e-5;
shear_rate = 10;
write_freq =2000;
[a_ij, totalDeformation ] =  compute_a_ij2D_fromFile(numFramesFileName, positionsFileName, dt, shear_rate, write_freq  );

figure
hold off
plot(totalDeformation, a_ij(:,1)) %a11
hold on
plot(totalDeformation, a_ij(:,3),'r') %a22
index = (floor(length(totalDeformation)/2));
%plot([totalDeformation(index) totalDeformation(index)], [0 1]   )
%plot([totalDeformation(end) totalDeformation(end)], [0 1]   )
legend('a_{11}','a_{22}')
set(gca,'FontSize', 18)
xlabel('Total Strain - \gamma ')
ylabel('Orientation Tensor - Component ')
ylim([0 1])


mean(a_ij(index:end,1))
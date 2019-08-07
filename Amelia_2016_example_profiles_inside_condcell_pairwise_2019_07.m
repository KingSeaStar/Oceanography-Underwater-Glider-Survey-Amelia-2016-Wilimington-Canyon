%% pair wise correction
% make a movie of each pair of profiles

%% movie for potential temperature comparison

n_iter = size(strictly_filtered_down_up_pair_indices, 1);

fig1 = figure('position',[200 200 1200 600]);
suptitle('Pair-wise correction: potential temperature')
hold on;
writerObj = VideoWriter('Potential_Temperature_Profile.avi');
writerObj.FrameRate = 5;
open(writerObj);

for ii = 1:n_iter
    
    subplot(1,2,1)
    plot(pair_wise_cor_downcast(ii).ptemp, pair_wise_cor_downcast(ii).z, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).ptemp, pair_wise_cor_upcast(ii).z, '.b')
    hold off
    grid on
    xlim([7, 16])
    ylim([-350, 0])
    xlabel('Potential Temperature (C)')
    ylabel('Depth (m)')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,2,2)
    plot(pair_wise_cor_downcast(ii).ptemp_cor_outside, pair_wise_cor_downcast(ii).z, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).ptemp_cor_outside, pair_wise_cor_upcast(ii).z, '.b')
    hold off
    grid on
    xlim([7, 16])
    ylim([-350, 0])
        xlabel('Potential Temperature (C)')
    ylabel('Depth (m)')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['after','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);


%% movie for practical salinity comparison

n_iter = size(strictly_filtered_down_up_pair_indices, 1);

fig1 = figure('position',[200 200 1200 600]);
suptitle('Pair-wise correction: Salinity')
hold on;
writerObj = VideoWriter('Salinity_Profile.avi');
writerObj.FrameRate = 5;
open(writerObj);

for ii = 1:n_iter
    
    subplot(1,2,1)
    plot(pair_wise_cor_downcast(ii).salt, pair_wise_cor_downcast(ii).z, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).salt, pair_wise_cor_upcast(ii).z, '.b')
    hold off
    grid on
    xlim([33 36])
    ylim([-350, 0])
    xlabel('Salinity')
    ylabel('Depth (m)')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,2,2)
    plot(pair_wise_cor_downcast(ii).salt_cor_inside, pair_wise_cor_downcast(ii).z, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).salt_cor_inside, pair_wise_cor_upcast(ii).z, '.b')
    hold off
    grid on
    xlim([33 36])
    ylim([-350, 0])
    xlabel('Salinity')
    ylabel('Depth (m)')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['after','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);



%% movie for potential density comparison

n_iter = size(strictly_filtered_down_up_pair_indices, 1);

fig1 = figure('position',[200 200 1200 600]);
suptitle('Pair-wise correction: Potential Density Anomaly')
hold on;
writerObj = VideoWriter('Density_Profile.avi');
writerObj.FrameRate = 5;
open(writerObj);

for ii = 1:n_iter
    
    subplot(1,2,1)
    plot(pair_wise_cor_downcast(ii).sigma0, pair_wise_cor_downcast(ii).z, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).sigma0, pair_wise_cor_upcast(ii).z, '.b')
    hold off
    grid on
    xlim([24.5 27.6])
    ylim([-350, 0])
    xlabel('Density Anomaly (Kg/m^{3})')
    ylabel('Depth (m)')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,2,2)
    plot(pair_wise_cor_downcast(ii).sigma0_cor_inside, pair_wise_cor_downcast(ii).z, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).sigma0_cor_inside, pair_wise_cor_upcast(ii).z, '.b')
    hold off
    grid on
     xlim([24.5 27.6])
    ylim([-350, 0])
    xlabel('Density Anomaly (Kg/m^{3})')
    ylabel('Depth (m)')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['after','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);

%% movie for TS diagram comparison

n_iter = size(strictly_filtered_down_up_pair_indices, 1);

fig1 = figure('position',[200 200 1200 600]);
suptitle('Pair-wise correction: TS diagram')
hold on;
writerObj = VideoWriter('TS_diagram.mov');
writerObj.FrameRate = 5;
open(writerObj);

for ii = 1:n_iter
    
    subplot(1,2,1)
    plot(pair_wise_cor_downcast(ii).salt, pair_wise_cor_downcast(ii).ptemp, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).salt, pair_wise_cor_upcast(ii).ptemp, '.b')
    hold off
    grid on
    xlim([33 36])
    ylim([7, 16])
    xlabel('Potential Temperature (C)')
    ylabel('Salinity')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,2,2)
    plot(pair_wise_cor_downcast(ii).salt_cor_inside, pair_wise_cor_downcast(ii).ptemp_cor_outside, '.r');
    hold on;
    plot(pair_wise_cor_upcast(ii).salt_cor_inside, pair_wise_cor_upcast(ii).ptemp_cor_outside, '.b')
    hold off
    grid on
    xlim([33 36])
    ylim([7, 16])
    xlabel('Potential Temperature (C)')
    ylabel('Salinity')
    legend({'downcast','upcast'}, 'location', 'southwest');
    title (sprintf(['after','\n',...
        'Pair:', num2str(ii)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);
clc;
clear;
close all;
format long;

%% Input Data
histflag = 1;
nvalsarr = [32,48,64,72,80,100]; ngraft = 64;
free_energy = zeros(10,2);
diff_energy = zeros(length(nvalsarr),4);
err_energy = zeros(length(nvalsarr),4);
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m',brown,green,'k','b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Plot Free Energy

fout = fopen('deltaF_all.dat','w');
fprintf(fout,'%s\t%s\t%s\t%s\n','Nfree','Arch','deltaF','Error');
for nvals = 1:length(nvalsarr)
    
    for i = 1:4
        
        if i == 1
            dirstr = 'block_block';
        elseif i == 2
            dirstr = 'block_alter';
        elseif i == 3
            dirstr = 'alter_block';
        else
            dirstr = 'alter_alter';
        end
        
        fylename = sprintf('./whamout_all/whamnew/n_%d_%s/whamout.txt',nvalsarr(nvals),dirstr);
        fid = fopen(fylename,'r');
        free_energy = zeros(10,3);
        header = fgetl(fid);
        k = 1;
        
        while ~feof(fid)
            
            tline = fgetl(fid);
            strarr = strsplit(tline);
            
            if strcmp(strarr{1},'#Window')
                break;
            else
                free_energy(k,1) = str2double(strarr{1});
                free_energy(k,2) = str2double(strarr{2});
                free_energy(k,3) = str2double(strarr{3});
                k = k + 1;
            end
            
        end
        
        [minfree,indmin] = min(free_energy(:,2));
        [maxfree,indmax] = max(free_energy(:,2));
        
        diff_energy(nvals,i) =  minfree - maxfree;
        err_energy(nvals,i) = sqrt(free_energy(indmin,3)^2 + free_energy(indmax,3)^2);        
        fprintf(fout,'%d\t%s\t%g\t%g\n',nvalsarr(nvals),dirstr,diff_energy(nvals,i),err_energy(nvals,i));
        
    end
    fclose(fid);
end



h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$\Delta F$ ($k_B T$)','FontSize',20,'Interpreter','Latex')

for i = 1:4
    errorbar(nvalsarr/ngraft,diff_energy(:,i),err_energy(:,i),'color',pclr{i},'LineWidth',2, ...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'delta_F','png');

%% Free_Energy Plots for n = 32

nvalsarr = [32];
h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$z$','FontSize',20,'Interpreter','Latex')
ylabel('$\Delta F(z)$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
for nvals = 1:length(nvalsarr)
    
    for i = 1:4
        
        if i == 1
            dirstr = 'block_block';
        elseif i == 2
            dirstr = 'block_alter';
        elseif i == 3
            dirstr = 'alter_block';
        else
            dirstr = 'alter_alter';
        end
        
        fylename = sprintf('./whamout_all/whamnew/n_%d_%s/whamout.txt',nvalsarr(nvals),dirstr);
        fid = fopen(fylename,'r');
        free_energy = zeros(10,3);
        header = fgetl(fid);
        
        while ~feof(fid)
            
            tline = fgetl(fid);
            strarr = strsplit(tline);
            
            if strcmp(strarr{1},'#Window')
                break;
            else
                free_energy(k,1) = str2double(strarr{1});
                free_energy(k,2) = str2double(strarr{2});
                free_energy(k,3) = str2double(strarr{3});
                k = k + 1;
            end
            
        end
        
        errorbar(free_energy(:,1),free_energy(:,2),free_energy(:,3),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{2},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
        
    end
    
    fclose(fid);
    
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';
xlim([3 60])

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','NorthWest')
legend boxoff
saveas(h1,'delta_Fz','png');


%% Histogram Plots

if histflag == 1
    nvalsarr = [32]; nvals = 1;
    winarr  = {'5.0';'8.0';'10.0';'12.0';'14.0';'18.0';'22.0';'26.0';'30.0';'32.0';'34.0';'38.0';'42.0';'46.0'};
    
    for i = 1:1
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$z$','FontSize',20,'Interpreter','Latex')
        ylabel('Probability','FontSize',20,'Interpreter','Latex')
        
        for wval = 1:length(winarr)
            
            if i == 1
                dirstr = 'block_block';
            elseif i == 2
                dirstr = 'block_alter';
            elseif i == 3
                dirstr = 'alter_block';
            else
                dirstr = 'alter_alter';
            end
            
            fylename = sprintf('./whamout_all/USbackup_colvar/n_%d/out.colvars.traj_%s_%s',...
                nvalsarr(nvals),dirstr,winarr{wval});
            fid = fopen(fylename,'r');
            trval = zeros(1000,2);
            tline = fgetl(fid);
            k = 1;
            
            while ~feof(fid)
                
                tline = fgetl(fid);
                strarr = strsplit(tline);
                
                if strcmp(strarr{1},'#')
                    continue;
                else
                    trval(k,1) = str2double(strarr{2});
                    trval(k,2) = str2double(strarr{3});
                    k = k + 1;
                end
                
            end
            
            histogram(trval(:,2),80)
            
        end
        
        fclose(fid);
        saveas(h1,'hist','png');
        
    end
    
    
    
    nvalsarr = [32]; nvals = 1;
    winarr  = {'5.0';'8.0';'10.0';'12.0';'14.0';'18.0';'22.0';'26.0';'30.0';'32.0';'34.0';'38.0';'39.0';'40.0';'41.0';'42.0';'46.0'};
    
    for i = 1:1
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$z$','FontSize',20,'Interpreter','Latex')
        ylabel('Probability','FontSize',20,'Interpreter','Latex')
        
        for wval = 1:length(winarr)
            
            if i == 1
                dirstr = 'block_block';
            elseif i == 2
                dirstr = 'block_alter';
            elseif i == 3
                dirstr = 'alter_block';
            else
                dirstr = 'alter_alter';
            end
            
            fylename = sprintf('./whamout_all/USbackup_colvar/n_%d/out.colvars.traj_%s_%s',...
                nvalsarr(nvals),dirstr,winarr{wval});
            fid = fopen(fylename,'r');
            trval = zeros(1000,2);
            tline = fgetl(fid);
            k = 1;
            
            while ~feof(fid)
                
                tline = fgetl(fid);
                strarr = strsplit(tline);
                
                if strcmp(strarr{1},'#')
                    continue;
                else
                    trval(k,1) = str2double(strarr{2});
                    trval(k,2) = str2double(strarr{3});
                    k = k + 1;
                end
                
            end
            
            histogram(trval(:,2),80)
            
        end
        
        fclose(fid);
        saveas(h1,'hist','png');
        
    end
    
end
    

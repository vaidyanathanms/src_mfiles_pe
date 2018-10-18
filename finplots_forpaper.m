%% Final Figure Plots

clear;
clc;
% close all;
format long;


%% Flags

fig2a = 0; %fads-Npa/Npc
fig2b = 0; %delF-Npa/Npc
fig3a = 1; %delU/delN-Npa/Npc
fig3b = 1; %delS-Npa/Npc
fig4a = 0; %rhopa,rhopc-z @n=32
fig4b = 0; %rhopa,rhopc-z @n=100
fig5a = 0; %q(z) @n=100
fig5b = 0; %Qb vs Npa/Npc

%% Inputs

nmonfree = 30; nmongraft = 30; ngraft = 64;nsalt=510;nbackmons=10;
nbase = nbackmons*ngraft;
nfreearr = [16;32;48;64;80;100;150];
cutoff = '1.50'; lz = 120; area=53^2;
rhofree = nfreearr*nmonfree/(lz*area);


green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m',brown,green,'k','b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

nadschain = zeros(length(nfreearr),4);

%% Fig 2a. Plot adsorbed fraction as a function of number of ADSORBED CHAINS
%  Ref:adsfrac.m
if fig2a == 1
    fout = fopen(sprintf('adsorbed_chain_ave_rcut_%s.dat',cutoff),'w');
    fprintf(fout,'%s\t%s\t%s\n','N_f','Arch','fraction');
    errvals = importdata('error_fvals.txt');
    for ncnt = 1:length(nfreearr)
        nval = nfreearr(ncnt);
        for i = 1:4
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            else
                dirstr = 'al_al';
            end
            
            filename = sprintf('./results_adsfrac/results_%d_%s/adsfracchain_rcut_%s.lammpstrj',...
                nval,dirstr,cutoff);
            data = importdata(filename);
            nadschain(ncnt,i) = mean(data(:,3));
            fprintf(fout,'%d\t%s\t%g\n',nval,dirstr,nadschain(ncnt,i));
        end
    end
    fclose(fout);
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{ads}$','FontSize',20,'Interpreter','Latex')
    
    errorbar(nfreearr/ngraft,nadschain(:,1),errvals.data(:,2),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
    errorbar(nfreearr/ngraft,nadschain(:,3),errvals.data(:,3),'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
    errorbar(nfreearr/ngraft,nadschain(:,2),errvals.data(:,4),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
    errorbar(nfreearr/ngraft,nadschain(:,4),errvals.data(:,5),'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Alter-Block';
    legendinfo{3} = 'Block-Alter';
    legendinfo{4} = 'Alter-Alter';
    
    
    %overlay y = x line
    
    x = 0:0.1:1.1; y = x;
    plot(x,y,'LineWidth',2,'Color',orange,'LineStyle','--')
    
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('adsorbchain_rcut_%s.png',cutoff));
    
    
    %overlay y = x line
    
    x = 0:0.1:1.1; y = x;
    plot(x,y,'LineWidth',2,'Color',orange,'LineStyle','--')
    
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('adsorbchain_rcut_%s.png',cutoff));
end

%% Fig2b. Plot Free Energy
% Ref: whamplots.m
if fig2b == 1
    nvalsarr = [32;48;64;72;80;100];
    diff_energy = zeros(length(nvalsarr),4);
    err_energy = zeros(length(nvalsarr),4);
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
    
end

%% Fig 3a. Intensive (delU/delN) with errorbar
% Ref: analyzeenergy.m

if fig3a == 1
    
    nvalsarr = [32,48,64,72];   
    
    all_Emeans = fopen('./biascalc/AllEner.dat','r');
    errUbyN = importdata('./biascalc/deludelN.txt');
    
    varUbyN  = zeros(length(nvalsarr),4);
    meanUbyN = zeros(length(nvalsarr),4);
    
    for i = 1:3
        tline = fgetl(all_Emeans);
    end
    
    lennvals = length(nvalsarr); arrcnt = 1; archcnt = 1;
    while ~feof(all_Emeans)
        tline = fgetl(all_Emeans);
        strarr = strsplit(tline);
        meanUbyN(arrcnt,archcnt) = str2double(strarr{6});
        if rem(arrcnt,lennvals) == 0
            arrcnt = 1; archcnt = archcnt + 1;
        else
            arrcnt = arrcnt + 1;
        end
    end
    
    %Cols 3-6 (UbyN value std)
    
    kvar = 3; %See above comment
    for i = 1:4
        varUbyN(:,i) = errUbyN.data(:,kvar);
        kvar = kvar + 1;
    end
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle \Delta U \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        errorbar(nvalsarr/ngraft,meanUbyN(:,i),varUbyN(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'deltaUbyN_witherr','png');
    fclose(all_Emeans)
end

%% Fig 3b. Entropy with errorbar
% Ref: finplots_energy.m

if fig3b == 1
    nvalsarr = [32,48,64,72]; ncutoff = '1.50';
    
    fylename = 'deltaF_all.dat';
    if exist(fylename, 'file') == 2
        
        all_Emeans = fopen('./biascalc/AllEner.dat','r');
        errUbyN = importdata('./biascalc/deludelN.txt');
        
        varUbyN  = zeros(length(nvalsarr),4);
        meanUbyN = zeros(length(nvalsarr),4);
        
        for i = 1:3
            tline = fgetl(all_Emeans);
        end
        
        lennvals = length(nvalsarr); arrcnt = 1; archcnt = 1;
        while ~feof(all_Emeans)
            tline = fgetl(all_Emeans);
            strarr = strsplit(tline);
            meanUbyN(arrcnt,archcnt) = str2double(strarr{6});
            if rem(arrcnt,lennvals) == 0
                arrcnt = 1; archcnt = archcnt + 1;
            else
                arrcnt = arrcnt + 1;
            end
        end
        
        %Cols 3-6 (N value std)
        
        kvar = 3; %See above comment
        for i = 1:4
            varUbyN(:,i) = errUbyN.data(:,kvar);
            kvar = kvar + 1;
        end

        
        dataN = importdata('./biascalc/delN.txt');
        %Cols 3-6 (N value std), Cols 10-13 (N value mean)
        
        varN  = zeros(length(nvalsarr),4);
        meanN = zeros(length(nvalsarr),4);
        
        kvar = 3; kmean = 10; %See above comment
        for i = 1:4
            varN(:,i) = dataN.data(:,kvar);
            meanN(:,i) = dataN.data(:,kmean);
            kvar = kvar + 1; kmean = kmean+1;
        end
        
        
        delF = zeros(length(nvalsarr),4);
        errF = zeros(length(nvalsarr),4);
        k = 1;
        ffree = fopen(fylename,'r');
        header = fgetl(ffree);
        free_energy = zeros(length(nvalsarr)*4,2);
        
        while ~feof(ffree) && k <= length(nvalsarr)*4
            tline = fgetl(ffree);
            strarr = strsplit(tline);
            free_energy(k,1) = str2double(strarr{3});
            free_energy(k,2) = str2double(strarr{4});
            k = k + 1;
        end
        fclose(ffree);
        
        if k ~= length(nvalsarr)*4+1
            error('Unequal number of free energy columns%d\t%d\n',k,length(nvalsarr)*4+1);
        end
        
        k = 1;
        for i = 1:length(nvalsarr)
            for j = 1:4
                delF(i,j) = free_energy(k,1);
                errF(i,j) = free_energy(k,2);
                k = k + 1;
            end
        end
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
        ylabel('$\langle T \Delta S \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            delS = meanUbyN(:,i) - delF(:,i);
            errfromU = varUbyN(:,i);
            errfromF = errF(:,i);
            errforS  = sqrt(errfromU.^2 + errfromF.^2);
            errorbar(nvalsarr/ngraft,delS,varUbyN(:,i),'color',pclr{i},'LineWidth',2, ...
                'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
            
            
        end
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,'deltaS','png');
        
        
    else
        fprintf('%s not found',fylename);
        error('Cannot compute Delta S');
    end
   
end

%% Fig 4a. Density Profiles @n-32
% Ref:paper_methods.m

if fig4a == 1
    
    nfree = 32;
    ncntr = abs(ngraft*nmongraft/2-nmonfree*nfree/2);
    ntotmons = nfree*nmonfree+ngraft*(nmongraft+nbackmons) + ncntr+nsalt*2;
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    h = zeros(4,1);
    for i = 1:4
        
        h(i) = subplot(4,1,i);
        
        if i == 1
            dirstr = 'bl_bl';
        elseif i == 2
            dirstr = 'bl_al';
        elseif i == 3
            dirstr = 'al_bl';
        elseif i == 4
            dirstr = 'al_al';
        else
            disp('No Correct String')
            break;
        end
        
        fid = fopen(sprintf('./results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata     = fld(:,1);
        pl_data(:,1) = fld(:,2)*ngraft*nmongraft/ntotmons;
        pl_data(:,2) = fld(:,3)*nfree*nmonfree/ntotmons;
        
        patch(rdata/lz,pl_data(:,1),gold);
        alpha(0.3);
        patch(rdata/lz,pl_data(:,2),'b')
        alpha(0.5);
        xlim([0 0.3])
        fclose(fid);
        
    end
    
    set(h(1),'xticklabel',[]);
    set(h(2),'xticklabel',[]);
    set(h(3),'xticklabel',[]);
    
    pos=get(h,'position');
    bottom=pos{4}(2);
    top=pos{1}(2)+pos{1}(4);
    plotspace=top-bottom;
    
    pos{4}(4)=plotspace/4;
    pos{3}(4)=plotspace/4;
    pos{2}(4)=plotspace/4;
    pos{1}(4)=plotspace/4;
    
    pos{1}(2)=bottom + 3*plotspace/4;
    pos{2}(2)=bottom + 2*plotspace/4;
    pos{3}(2)=bottom + plotspace/4;
    
    
    set(h(1),'position',pos{1},'FontSize',16);
    set(h(2),'position',pos{2},'FontSize',16);
    set(h(3),'position',pos{3},'FontSize',16);
    set(h(4),'position',pos{4},'FontSize',16);
    box on
    
    yticks(h(4),[])
    if(nfree >= 80)
        yticks(h(1),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(1),[0 0.25 0.5 0.75])
        yticks(h(2),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(2),[0 0.25 0.5 0.75])
        yticks(h(3),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(3),[0 0.25 0.5 0.75])
        yticks(h(4),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(4),[0 0.25 0.5 0.75])
    else
        yticks(h(1),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(1),[0 0.5 1.0 1.5])
        yticks(h(2),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(2),[0 0.5 1.0 1.5])
        yticks(h(3),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(3),[0 0.5 1.0 1.5])
        yticks(h(4),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(4),[0 0.5 1.0 1.5])
        
    end
    
    xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho(z) \times 10^5$','FontSize',20,'Interpreter','Latex')
    box(h(1),'on');
    box(h(2),'on');
    box(h(3),'on');
    legendinfo{1} = '$\rho_{pc}(z)$';
    legendinfo{2} = '$\rho_{pa}(z)$';
    legend(h(1), legendinfo, 'Interpreter','Latex','FontSize',16)
    box(legend,'off')
end

%% Fig 4b. Density Profiles
% Ref:paper_methods.m

if fig4b == 1
    
    nfree = 100;
    ncntr = abs(ngraft*nmongraft/2-nmonfree*nfree/2);
    ntotmons = nfree*nmonfree+ngraft*(nmongraft+nbackmons) + ncntr+nsalt*2;
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    h = zeros(4,1);
    for i = 1:4
        
        h(i) = subplot(4,1,i);
        
        if i == 1
            dirstr = 'bl_bl';
        elseif i == 2
            dirstr = 'bl_al';
        elseif i == 3
            dirstr = 'al_bl';
        elseif i == 4
            dirstr = 'al_al';
        else
            disp('No Correct String')
            break;
        end
        
        fid = fopen(sprintf('./results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata     = fld(:,1);
        pl_data(:,1) = fld(:,2)*ngraft*nmongraft/ntotmons;
        pl_data(:,2) = fld(:,3)*nfree*nmonfree/ntotmons;
        
        patch(rdata/lz,pl_data(:,1),gold);
        alpha(0.3);
        patch(rdata/lz,pl_data(:,2),'b')
        alpha(0.5);
        xlim([0 0.3])
        fclose(fid);
        
    end
    
    set(h(1),'xticklabel',[]);
    set(h(2),'xticklabel',[]);
    set(h(3),'xticklabel',[]);
    
    pos=get(h,'position');
    bottom=pos{4}(2);
    top=pos{1}(2)+pos{1}(4);
    plotspace=top-bottom;
    
    pos{4}(4)=plotspace/4;
    pos{3}(4)=plotspace/4;
    pos{2}(4)=plotspace/4;
    pos{1}(4)=plotspace/4;
    
    pos{1}(2)=bottom + 3*plotspace/4;
    pos{2}(2)=bottom + 2*plotspace/4;
    pos{3}(2)=bottom + plotspace/4;
    
    
    set(h(1),'position',pos{1},'FontSize',16);
    set(h(2),'position',pos{2},'FontSize',16);
    set(h(3),'position',pos{3},'FontSize',16);
    set(h(4),'position',pos{4},'FontSize',16);
    box on
    
    yticks(h(4),[])
    if(nfree >= 80)
        yticks(h(1),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(1),[0 0.25 0.5 0.75])
        yticks(h(2),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(2),[0 0.25 0.5 0.75])
        yticks(h(3),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(3),[0 0.25 0.5 0.75])
        yticks(h(4),[0 2.5e-6 5e-6 7.5e-6])
        yticklabels(h(4),[0 0.25 0.5 0.75])
    else
        yticks(h(1),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(1),[0 0.5 1.0 1.5])
        yticks(h(2),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(2),[0 0.5 1.0 1.5])
        yticks(h(3),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(3),[0 0.5 1.0 1.5])
        yticks(h(4),[0 5e-6 1e-5 1.5e-5])
        yticklabels(h(4),[0 0.5 1.0 1.5])
        
    end
    
    xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho(z) \times 10^5$','FontSize',20,'Interpreter','Latex')
    box(h(1),'on');
    box(h(2),'on');
    box(h(3),'on');
    legendinfo{1} = '$\rho_{pc}(z)$';
    legendinfo{2} = '$\rho_{pa}(z)$';
    legend(h(1), legendinfo, 'Interpreter','Latex','FontSize',16)
    box(legend,'off')
end
%% Fig 5a. Charge density, q(z)
% Ref:densananew.m

if fig5a == 1
    
    nfree = 100;rbin = 1;nbins=lz/rbin;
    ncntr = abs((nfree-ngraft)*nmonfree/2);
    ntotal = (nfreearr + ngraft)*nmonfree + 2*nsalt + ncntr + nbase;
    chargearr = [0;1;0;-1;1;-1];
    
    nval = nfree;
    if nval > ngraft
        cntrcharg = 1;
    elseif nval == ngraft
        cntrcharg = 0;
    elseif nval < ngraft
        cntrcharg = -1;
    end
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$q(z)$','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        if i == 1
            dirstr = 'bl_bl';
        elseif i == 2
            dirstr = 'bl_al';
        elseif i == 3
            dirstr = 'al_bl';
        elseif i == 4
            dirstr = 'al_al';
        else
            disp('No Correct String')
            break;
        end
        fid = fopen(sprintf('./results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata     = fld(:,1);
        pneutral_brush = fld(:,2)*ngraft*nmongraft/2*lz*area/nbins;
        pnegativ_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
        pneutral_free  = fld(:,4)*nfree*nmonfree/2*lz*area/nbins;
        ppositive_free = fld(:,5)*nfree*nmonfree/2*lz*area/nbins;
        pposions  = fld(:,6)*nsalt*lz*area/nbins;
        pnegions  = fld(:,7)*nsalt*lz*area/nbins;
        pcntrions = fld(:,8)*ncntr*lz*area/nbins;
        
        if nval ~= ngraft
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
        else
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions;
        end
        
        plot(rdata/lz,netcharge, 'color', pclr{i}, 'LineWidth', 2)
        
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('n%d_integral_netcharge',nval),'png');
    
end

%% Fig 5b: Qb - Npc/Npa

if fig5b == 1
    
    chargearr = [0;1;0;-1;1;-1];
    rbin = 1;nbins=lz/rbin;
    ncntr = abs((nfreearr*nmonfree-ngraft*nmongraft)/2);
    ntotarr = nfreearr*nmonfree + nmongraft*ngraft + nbase + nsalt + ncntr;
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$f$','FontSize',20,'Interpreter','Latex')
    
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('./fig1b_QnetBound_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','NetCharge: Q_{b}=\Delta(n_g)*\int(\sum(q_j n_j(z)dz, j=all entities)z=0,Lz))');
        fprintf(fnr,'%s\t %s\n','Arch','Q_{b}');
        
        if nfree > ngraft
            cntrcharg = 1;
        elseif nfree == ngraft
            cntrcharg = 0;
        elseif nfree < ngraft
            cntrcharg = -1;
        end
        
        for i = 1:4
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            elseif i == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            
            fid = fopen(sprintf('./results_dens/results_%d_%s/dens.lammpstrj',nfree,dirstr));
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            pneutral_brush = fld(:,2)*ngraft*nmongraft/2*lz*area/nbins;
            pcation_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
            pneutral_free  = fld(:,4)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            panion_brush = fld(:,5)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            pposions  = fld(:,6)*nsalt*lz*area/nbins;
            pnegions  = fld(:,7)*nsalt*lz*area/nbins;
            pcntrions = fld(:,8)*ncntr(nvals)*lz*area/nbins;
            
            if nfree*nmonfree ~= ngraft*nmongraft
                netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pcation_brush + ...
                    chargearr(3)*pneutral_free + chargearr(4)*panion_brush + ...
                    chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
            else
                netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pcation_brush + ...
                    chargearr(3)*pneutral_free + chargearr(4)*panion_brush + ...
                    chargearr(5)*pposions + chargearr(6)*pnegions;
            end
            
            sumq = 0.5*netcharge(1);
            
            fid_g  = fopen(sprintf('./results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
            data_g = textscan(fid_g,'%f%f%f','Headerlines',1);
            
            fld_g   = cell2mat(data_g);
            dens_g  = fld_g(:,2);
            [maxden, imaxden]  = max(dens_g);
            
            trapz(rdata,pcntrions)
            
            % Find edge of brush
            
            for k = imaxden:length(dens_g)
                
                if dens_g(k,1) < 0.05*maxden
                    
                    i_edge = k;
                    break;
                    
                end
                
            end
            
            qofr = zeros(i_edge,1); rqofr = zeros(i_edge,1);
            qofr(1,1) = sumq; rqofr(1,1) = 0.5*rdata(1);
            
            % Integrate charge to the edge of the brush
            for k = 2:i_edge
                
                sumq = sumq + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
                qofr(k,1)  = sumq;
                rqofr(k,1) = 0.5*(rdata(k)+rdata(k+1));
                
            end
            
            fprintf(fnr,'%s\t%g\n',dirstr,sumq);
            fclose(fid_g);
            fclose(fid);
            
            if(nvals == 1)
                
                plot(rqofr,qofr,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{1})
                
            end
            
        end
        
        fclose(fnr);
        
    end
    
    % Plot Q_{b} data
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('./fig1b_QnetBound_%d.txt',nfreearr(i)),'r');
        data = textscan(fnr,'%s%f','Headerlines',2);
        
        fld = cell2mat(data(2));
        data_bb(i,1) = fld(1,1);
        data_ba(i,1) = fld(2,1);
        data_ab(i,1) = fld(3,1);
        data_aa(i,1) = fld(4,1);
        
        fclose(fnr);
        
    end
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$Q_{b}$','FontSize',20,'Interpreter','Latex')
    
    plot(nfreearr/ngraft,data_bb,'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
    plot(nfreearr/ngraft,data_ba,'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
    plot(nfreearr/ngraft,data_ab,'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
    plot(nfreearr/ngraft,data_aa,'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
    
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'Fig1b_QnetBound_%s','png');
    clear legendinfo
    
end
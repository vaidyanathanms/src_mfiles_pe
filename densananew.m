%% To analyze the density profiles

clear;
clc;
close all;
format long;

%% Plot Graft Profiles


nfree = [32;64;80]; lz = 120; 
ngrafts = 64;area = 53^2;rbin = 1;nbins=lz/rbin;
nbase = 10*ngrafts;
nmons = 30;
nsalt = 510;
ncntr = abs((nfree-ngrafts)*nmons/2);
ntotal = (nfree + ngrafts)*nmons + 2*nsalt + ncntr + nbase;
pclr = {'r','b','k','m','g'};
nclr = {'--r','--b','--k','--m','--g'};
cclr = {'r*','b*','k*','m*','g*'};

% For a particular N


for j = 1:length(nfree)
    nval = nfree(j);
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho_{f}(r)$','FontSize',20,'Interpreter','Latex')
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
        sprintf('./n_%d/results_%d_%s/grpdens.lammpstrj',nval,nval,dirstr)
        fid = fopen(sprintf('./n_%d/results_%d_%s/grpdens.lammpstrj',nval,nval,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata   = fld(:,1);
        pegraft = fld(:,2);
        pefree  = fld(:,3);
        plot(rdata/lz,pefree, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('n%d_free',nval),'png');
    
end

% Plot Ion density Profiles


for j = 1:length(nfree)
    nval = nfree(j);
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho^{+}(r)$','FontSize',20,'Interpreter','Latex')
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
        fid = fopen(sprintf('./n_%d/results_%d_%s/dens.lammpstrj',nval,nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata     = fld(:,1);
        pposions  = fld(:,6);
        pnegions  = fld(:,7);
        pcntrions = fld(:,8);
        ax(i) = plot(rdata/lz,pposions, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8);
        %         plot(rdata/lz,pnegions, nclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
%         plot(rdata/lz,pcntrions, cclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('n%d_posions',nval),'png');
    
end


% Plot Overall Charge

chargearr = [0;1;0;-1;1;-1];


for j = 1:length(nfree)
    nval = nfree(j);
    if nval > ngrafts
        cntrcharg = 1;
    elseif nval == ngrafts
        cntrcharg = 0;
    elseif nval < ngrafts
        cntrcharg = -1;
    end
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$Q(r)$','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        if i == 1
            dirstr = 'blbl';
        elseif i == 2
            dirstr = 'blal';
        elseif i == 3
            dirstr = 'albl';
        elseif i == 4
            dirstr = 'alal';
        else
            disp('No Correct String')
            break;
        end
        fid = fopen(sprintf('./n_%d/n%d_results/dens_%s.txt',nval,nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata     = fld(:,1);
        pneutral_brush = fld(:,2)*ngrafts*nmons/2;
        pnegativ_brush = fld(:,3)*ngrafts*nmons/2;
        pneutral_free  = fld(:,4)*nfree(j)*nmons/2;
        ppositive_free = fld(:,5)*nfree(j)*nmons/2;
        pposions  = fld(:,6)*nsalt;
        pnegions  = fld(:,7)*nsalt;
        pcntrions = fld(:,8)*ncntr(j);
        
        if nval ~= ngrafts
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
        else
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions;
        end
        plot(rdata/lz,netcharge, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('n%d_netcharge',nval),'png');
    
end


% Plot Integral Charge

chargearr = [0;1;0;-1;1;-1];


for j = 1:length(nfree)
    nval = nfree(j);
    if nval > ngrafts
        cntrcharg = 1;
    elseif nval == ngrafts
        cntrcharg = 0;
    elseif nval < ngrafts
        cntrcharg = -1;
    end
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$ Q(r) $','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        if i == 1
            dirstr = 'blbl';
        elseif i == 2
            dirstr = 'blal';
        elseif i == 3
            dirstr = 'albl';
        elseif i == 4
            dirstr = 'alal';
        else
            disp('No Correct String')
            break;
        end
        fid = fopen(sprintf('./n_%d/n%d_results/dens_%s.txt',nval,nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata     = fld(:,1);
        pneutral_brush = fld(:,2)*ngrafts*nmons/2*lz*area/nbins;
        pnegativ_brush = fld(:,3)*ngrafts*nmons/2*lz*area/nbins;
        pneutral_free  = fld(:,4)*nfree(j)*nmons/2*lz*area/nbins;
        ppositive_free = fld(:,5)*nfree(j)*nmons/2*lz*area/nbins;
        pposions  = fld(:,6)*nsalt*lz*area/nbins;
        pnegions  = fld(:,7)*nsalt*lz*area/nbins;
        pcntrions = fld(:,8)*ncntr(j)*lz*area/nbins;
        
        if nval ~= ngrafts
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
        else
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions;
        end
        
        intnet = zeros(length(netcharge),1);
        intnet(1,1) = 0.5*netcharge(1);
        rnet(1,1) = 0.5*(rdata(1));
        for k = 2:length(netcharge)
            intnet(k,1) = intnet(k-1,1) + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
            rnet(k,1) = 0.5*(rdata(k-1,1)+rdata(k,1));
        end
        plot(rnet/lz,intnet, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('n%d_integral_netcharge',nval),'png');
    
end


% Plot Integral Charge of polymers alone

chargearr = [0;1;0;-1;1;-1];


for j = 1:length(nfree)
    nval = nfree(j);
    if nval > ngrafts
        cntrcharg = 1;
    elseif nval == ngrafts
        cntrcharg = 0;
    elseif nval < ngrafts
        cntrcharg = -1;
    end
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$ Q(r) $','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        if i == 1
            dirstr = 'blbl';
        elseif i == 2
            dirstr = 'blal';
        elseif i == 3
            dirstr = 'albl';
        elseif i == 4
            dirstr = 'alal';
        else
            disp('No Correct String')
            break;
        end
        fid = fopen(sprintf('./n_%d/n%d_results/dens_%s.txt',nval,nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata     = fld(:,1);
        pneutral_brush = fld(:,2)*ngrafts*nmons/2*lz*area/nbins;
        pnegativ_brush = fld(:,3)*ngrafts*nmons/2*lz*area/nbins;
        pneutral_free  = fld(:,4)*nfree(j)*nmons/2*lz*area/nbins;
        ppositive_free = fld(:,5)*nfree(j)*nmons/2*lz*area/nbins;
        pposions  = fld(:,6)*nsalt*lz*area/nbins;
        pnegions  = fld(:,7)*nsalt*lz*area/nbins;
        pcntrions = fld(:,8)*ncntr(j)*lz*area/nbins;
        
        if nval ~= ngrafts
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
        else
            netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
                chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
                chargearr(5)*pposions + chargearr(6)*pnegions;
        end
        
        intnet = zeros(length(netcharge),1);
        intnet(1,1) = 0.5*netcharge(1);
        rnet(1,1) = 0.5*(rdata(1));
        for k = 2:length(netcharge)
            intnet(k,1) = intnet(k-1,1) + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
            rnet(k,1) = 0.5*(rdata(k-1,1)+rdata(k,1));
        end
        plot(rnet/lz,intnet, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('n%d_integral_netcharge',nval),'png');
    
end



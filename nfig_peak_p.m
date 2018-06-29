%10/24/2017 creation
%get f(k)

% interpreter for text
set(0,'DefaultTextInterpreter','latex');        
% interpreter for legend
set(0,'DefaultLegendInterpreter','latex');
% Figure parameters
nlin2=2; nlin2=1; nlin3=0.5;
nmark=1.5;
nfont1=14; nfont2=12; nfont3=10;
nmark=1.5;
col=[0.5 0.5 1];

lsty{1}='k-'; lsty{2}='r-'; lsty{3}='b-';lsty{4}='c-';
lsty{5}='g-'; lsty{6}='m-'; lsty{7}='y-';lsty{8}='k:';

lsty2{1}='k--'; lsty2{2}='r--'; lsty2{3}='b--';lsty2{4}='c--';
lsty2{5}='g--'; lsty2{6}='m--'; lsty2{7}='y--';lsty2{8}='k--';

num = 1:39;
% kn_res = sqrt(1/3*Eg/Ebind*num);
% kn_res1 = kn_res(2:end);
% kn_res2 = kn_res(1:(end-1));
% kn_err = kn_res1 - kn_res2;
% nkn_res = [0,kn_res2 - kn_err/2,kmax];
nkn_res = [0,find(kn_res<kmax),kmax];
i1 = size(nkn_res);
peak_f = zeros(i1(2)-1, Nm_o, Nt);
for N_t = 1:Nt
    for i = 1:(i1(2)-1)
        range = find((kgrid>nkn_res(i))&(kgrid<nkn_res(i+1)));
        peak_f(i,:, N_t) = max(real(fk((Nm_o*N_t-Nm_o+1):(Nm_o*N_t),range)'));
    end
end

fs_0 = abs(fk(Nm_o*(1:Nt)-Nm-2,1));
fs_1 = squeeze(peak_f(2,Nm-1,:));
fs_2 = squeeze(peak_f(3,Nm-1,:));
fs_r = sum(squeeze(peak_f(4:end,Nm-1,:)))';

peak_p = zeros(i1(2)-1, Nm_o, Nt);
for N_t = 1:Nt
    for i = 1:(i1(2)-1)
        range = find((kgrid>nkn_res(i))&(kgrid<nkn_res(i+1)));
        peak_p(i,:, N_t) = max(abs(pk((Nm_o*N_t-Nm_o+1):(Nm_o*N_t),range))');
    end
end

ps_0 = abs(pk(Nm_o*(1:Nt)-Nm-2,1));
ps_1 = squeeze(peak_p(2,Nm-1,:));
ps_2 = squeeze(peak_p(3,Nm-1,:));
ps_r = sum(squeeze(peak_p(4:end,Nm-1,:)))';

figure
subplot(2,1,1)
h = area(t, abs(E_en)/max(abs(E_en))*max([fs_0', fs_1', fs_2', fs_r']));
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on
plot(t, fs_0, lsty{2}, 'Linewidth', nlin1)
hold on
plot(t, fs_1, lsty{3}, 'Linewidth', nlin1)
hold on
plot(t, fs_2, lsty{4}, 'Linewidth', nlin1)
hold on
plot(t, fs_r, lsty{5}, 'Linewidth', nlin1)
legend('Electrical field', '0th(k = 0)()', '1st', '2nd', 'rest')
legend boxoff
num1 = 1:16;
for i= 1:num1(end)
    plot([(i/num1(end) - 1/2)*t_end, (i/num1(end) - 1/2)*t_end], [0,max([fs_0', fs_1', fs_2', fs_r'])], lsty{6}, 'Linewidth', nlin2)
    hold on
end
title(['d state, dt = ',num2str(params(2),3), ' as, dy = ',num2str(dk), ', Emax = ', num2str(params(6),3), 'MV/cm, $\gamma$ = ', num2str(params(7)), 'meV, run',num2str(num_run)],'FontName','Arial','FontSize',nfont1)
xlabel('t(ps)','FontName','Arial','FontSize',nfont2)
ylabel('$f_d$, peak','FontName','Arial','FontSize',nfont2)

subplot(2,1,2)
h = area(t, abs(E_en)/max(abs(E_en))*max([ps_0', ps_1', ps_2', ps_r']));
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on
plot(t, ps_0, lsty{2}, 'Linewidth', nlin1)
hold on
plot(t, ps_1, lsty{3}, 'Linewidth', nlin1)
hold on
plot(t, ps_2, lsty{4}, 'Linewidth', nlin1)
hold on
plot(t, ps_r, lsty{5}, 'Linewidth', nlin1)
legend('Electrical field', '0th(k = 0)()', '1st', '2nd', 'rest')
legend boxoff
for i= 1:num1(end)
    plot([(i/num1(end) - 1/2)*t_end, (i/num1(end) - 1/2)*t_end], [0,max([ps_0', ps_1', ps_2', ps_r'])], lsty{6}, 'Linewidth', nlin2)
    hold on
end
xlabel('t(ps)','FontName','Arial','FontSize',nfont2)
ylabel('$|p_d|$, peak','FontName','Arial','FontSize',nfont2)

h=gcf;
% set(h, 'PaperOrientationMode', 'auto');
% set(h, 'PaperOrientation', 'landscape');
print(h, '-dpdf', '-fillpage', ['/home/liuhai/dynamics/Cir_Harm_Exp/multi_photon/run',num2str(num_run), '/nfig_peak_d'])
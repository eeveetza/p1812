% testing against ITM
d = union(linspace(0,100,101), linspace(100, 1000, 91));
h = zeros(size(d));
R = zeros(size(d));
Ct = 2*ones(size(d)); % open
z = 4*ones(size(d));  % inland
GHz = 6;
Tpc = 50;
Phire = 70.5;
Phirn = -60;
Phite = 75;
Phitn = -60;
Re = 6371;

Hrg = 1.5;
Htg = 1.5;
Grx = 0;
Gtx = 0;
FlagVP = 1; % 1 - horizontal, 2 - vertical
dct = 500;
dcr = 500;
press = 1013;
temp = 20;

count = 1;
imin = 6;
for ii = imin:length(d)
    dd = d(1:ii);
    hh = h(1:ii);
    zz = z(1:ii);
    RR = R(1:ii);
    CCt = Ct(1:ii);

dpnt = dd(end)-dd(1);
[Phire_ii, Phirn_ii, bt2r, dgc] = great_circle_path(Phire, Phite, Phirn, Phitn, Re, dpnt);
    pdr = 0;
    L0(count) = tl_p1812(GHz, Tpc, dd, hh, RR, CCt, zz, Htg, Hrg, FlagVP, pdr, 'lam_t', Phite, 'phi_t', Phitn, 'lam_r', Phire_ii, 'phi_r', Phirn_ii, 'Gt', Gtx, 'Gr', Grx);
    pdr = 1;
    L1(count) = tl_p1812(GHz, Tpc, dd, hh, RR, CCt, zz, Htg, Hrg, FlagVP, pdr, 'lam_t', Phite, 'phi_t', Phitn, 'lam_r', Phire_ii, 'phi_r', Phirn_ii, 'Gt', Gtx, 'Gr', Grx);
    Lfs(count) = 92.4 + 20*log10(GHz) + 10*log10(dd(end).^2 + (Htg-Hrg).^2/1e6);
    count = count + 1;
end

plot(d(imin:end), L0, 'b', 'LineWidth', 2)
hold on
plot(d(imin:end), L1, 'r', 'LineWidth', 2)
grid on
plot(d(imin:end), Lfs, 'g', 'LineWidth', 2)
set(gca,'FontSize', 14)
legend('P.1812-6','P.1812 (new troposcatter)', 'Free-Space')
xlabel('distance (km)')
ylabel('Lb (dB)')
titlestr = ['f = ' num2str(GHz) ' GHz, Htg = ' num2str(Htg) ' m, Hrg = ' num2str(Hrg) ' m' ]
title(titlestr)
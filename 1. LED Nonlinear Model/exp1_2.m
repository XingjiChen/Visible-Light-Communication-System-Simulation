clear;
R = 1;
data = [0.00001 0.015 0.03 0.055 0.09 0.12 0.15 0.19 0.21 0.35 0.475 0.65 0.8 1];
minSumab = 10000;

for Imax = 1:0.01:3
    for k = 0:0.01:5
        v_LED = [2.77 2.8 2.85 2.9 2.95 3 3.05 3.08 3.1 3.2 3.3 3.4 3.5 3.58];
        f_LED = ((1.197.*(v_LED-2.77).^(2))+0.2537.*(v_LED-2.77)+0.001228)./R;
        i_LED = ((f_LED)./((1+((f_LED./Imax).^(2.*k))).^(1./(2.*k))));
        value = (i_LED-data);
        Sumab = sum(abs(value));
        if minSumab>Sumab
             minSumab = Sumab;
             Imaxover = Imax;
             kover = k;
        end
    end
end

vv_LED = [2.77:0.01:5.5];
ff_LED = ((1.197.*(vv_LED-2.77).^(2))+0.2537.*(vv_LED-2.77)+0.001228)./R;
ii_LED = ((ff_LED)./((1+((ff_LED./Imaxover).^(2.*kover))).^(1./(2.*kover))));
plot(vv_LED,ii_LED,'r')
set(gca,'XLim',[2.77 5.5])
hold on
xlabel('输入电压 v_{LED}  (V)')
ylabel('驱动电流 i_{LED}  (A)')
title('LED非线性伏安模型')

v_data=[2.77 2.8 2.85 2.9 2.95 3 3.05 3.08 3.1 3.2 3.3 3.4 3.5 3.58];
i_data=[0.00001 0.015 0.03 0.055 0.09 0.12 0.15 0.19 0.21 0.35 0.475 0.65 0.8 1];
plot(v_data,i_data,'b')
legend('Curve Fitting','Data Sheet')
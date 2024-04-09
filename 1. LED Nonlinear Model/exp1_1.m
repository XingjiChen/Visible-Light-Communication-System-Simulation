clear;

I_max = 0.5;
R = 1;
k = 2; 
v_LED = [-0.5:0.1:1];
v_len = length(v_LED);
f_LED = v_LED/R;
i_LED = (0).*(v_LED < 0)+((f_LED)./((1+((f_LED./I_max).^(2.*k))).^(1./(2.*k)))).*(v_LED >= 0);
plot(v_LED,i_LED,'r')
set(gca,'XLim',[-0.5 1])
hold on
xlabel('输入电压 v_{LED} (V)')
ylabel('驱动电流 i_{LED} (A)')
title('LED非线性伏安模型')

k = 3; 
f_LED = v_LED/R;
i_LED = (0).*(v_LED < 0)+((f_LED)./((1+((f_LED./I_max).^(2.*k))).^(1./(2.*k)))).*(v_LED >= 0);
plot(v_LED,i_LED,'b')
set(gca,'XLim',[-0.5 1])

k = 50; 
f_LED = v_LED/R;
i_LED = (0).*(v_LED < 0)+((f_LED)./((1+((f_LED./I_max).^(2.*k))).^(1./(2.*k)))).*(v_LED >= 0);
plot(v_LED,i_LED,'g')
set(gca,'XLim',[-0.5 1])
legend('k=2','k=3','k=50')
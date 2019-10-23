%script to look at angular VSF of different Duke scientific beads

%3000 series https://www.thermofisher.com/order/catalog/product/3269A?SID=srch-srp-3269A
D=[0.100, 0.205, 0.296, 0.400, 0.498, 0.600, 0.707, 0.799, 0.903];
dD=[0.003, 0.005, 0.006, 0.009, 0.009, 0.009, 0.009, 0.009, 0.012];

for i=1:length(D)
    [VSF(i,:)]=VSF_beads(D(i), dD(i), 700, 20, [90:1:180])
    plot([90:1:180],VSF(i,:)./max(VSF(i,:)),'-')
    hold on
end
xlabel('angle [degrees]');
ylabel('VSF/max(VSF)')
legend('0.1', '0.2','0.3', '0.5','0.6', '0.7','0.8', '0.9')

    
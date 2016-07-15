function show_endmembers(M2,wl)
M = size(M2,1);
names = cell(1,M);
for i = 1:M
    names{i} = ['endmember ',num2str(i)];
end
figure,plot(repmat(wl',1,size(M2,1)),M2','LineWidth',1);
xlabel('Wavelength (micrometer)');
ylabel('Reflectance');
legend(names);

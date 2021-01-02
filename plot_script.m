result = load('result.mat').psnr_record;
result = result(1:10,:,:);

% 画psnr变化图
mean_psnr = squeeze(mean(result,1));
M = 32:32:192;
sampling_rate = (M*512+128*128)/512/512;
phi = ["-","--"];
spa = ["o","x"];
sol = ["r","g","b"];
figure
for i=1:2
    for j=1:2
        for k=1:3
            plot(sampling_rate,mean_psnr(1:6,(i-1)*6+(j-1)*3+k),phi(i)+spa(j)+sol(k));
            hold on;
        end
    end
end
xlabel('samplint rate')
ylabel('psnr')
saveas(gcf,"psnr_sampling.png")

% 画箱式图
sol_set= [ "BP","OMP","Bayes"];
for i=1:3
    pos = i:3:12;
    figure
    ylabel('psnr')
    boxplot(squeeze(result(:,1,pos)),["Gauss,DCT,"+sol_set(i),"Gauss,DWT,"+sol_set(i),"Hada,DCT,"+sol_set(i),"Hada,DWT,"+sol_set(i)])
    saveas(gcf,sol_set(i)+"_boxplot.png");
end
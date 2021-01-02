repeat = 10;
psnr_record = zeros(10,6,12);

for i=1:repeat
    for j=32:32:192
        fprintf("%d,%d\n",i,j);
        psnr_record(i,j/32,:) = image_reconstruction_CS(j);
    end
end
save('result.mat','psnr_record')

function [ psnr_set ] = image_reconstruction_CS( M )
% 利用压缩感知进行图像重建测试
% 稀疏基：余弦变换 and 小波变换
% 观测矩阵：高斯随机 and 局部哈达玛矩阵
% 重建方法：基追踪 and 正交匹配追踪 and 贝叶斯估计
% M: 观测矩阵的行数
warning('off','all')
img = imread('lena512.bmp');  % 读入图像
img = im2double(img);
[h,w] = size(img);

% 构建高斯观测矩阵
phi_gauss = randn(M,h)./sqrt(M);
% 构建局部哈达玛观测矩阵
rand_row = randperm(h,M);
hadamard_mat = hadamard(h);
phi_hadamard = hadamard_mat(rand_row(1:M),:);

phi_set= {phi_gauss,phi_hadamard};

% 构建稀疏化向量
% DWT
[C, S] = wavedec2(img,2,'haar');
CA = C(1:S(1)*S(2));
C(1:S(1)*S(2))=1e-9;
img_dwt = reshape(C,h,w);
% DCT
img_dct = dct2(img);
CT = img_dct(1:S(1),1:S(2));
img_dct(1:S(1),1:S(2)) = 1e-9;
%imshow(img_dct);
sparse_set={img_dct,img_dwt}; 

% 求解器
sol_set= ["BP","OMP","Bayes"];

psnr_set = [];
for phi_ind = 1:length(phi_set)
    phi = phi_set{phi_ind};
    for sparse_ind = 1:length(sparse_set)
        % y = phi*theta
        img_compressed = phi*sparse_set{sparse_ind};
        for sol_ind = 1:length(sol_set)
            rec_tmp = zeros(h,w);
            for i=1:w
                if sol_set(sol_ind) == "OMP"
                    rec_tmp(:,i) = cs_omp(img_compressed(:,i),phi);
                elseif sol_set(sol_ind) == "BP"
                    rec_tmp(:,i) = cs_bp(img_compressed(:,i),phi);
                elseif sol_set(sol_ind) == "Bayes"
                    rec_tmp(:,i) = cs_bayes(img_compressed(:,i),phi);
                end
            end
            if sparse_ind == 1
                rec_tmp(1:S(1),1:S(2)) = CT;
                img_rec = idct2(rec_tmp);
            else
                rec_tmp = reshape(rec_tmp,1,h*w);
                rec_tmp(1:S(1)*S(2))= CA;
                img_rec = waverec2(rec_tmp,S,'haar');
            end
            %subplot(3,4,4*(sol_ind-1)+phi_ind*2-1+sparse_ind-1);
            % 计算psnr
            [peak_snr,~] = psnr(img_rec,img);
            psnr_set = [psnr_set, peak_snr];
            %fprintf("phi=%d,sparse=%d,%s,%0.4f\n",phi_ind,sparse_ind,sol_set(sol_ind),peak_snr);
            %imwrite(img_rec,sprintf("phi=%d,sparse=%d,%s.png",phi_ind,sparse_ind,sol_set(sol_ind)))
            %imshow(img_rec);
            %title(sprintf("phi=%d,sparse=%d,%s",phi_ind,sparse_ind,sol_set(sol_ind)));
        end
    end
end
end
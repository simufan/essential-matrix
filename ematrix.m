
clc; clear all; close all;

% Corresponding points between two images

temp = dir(['C:/Users/80677/Desktop/8point-algorithm-master/coresponding/','*']);

N = length(temp);   
dis_e1 = zeros(18, 0);
dis_e2 = zeros(18, 0);
rep_e1 = zeros(18, 0);
rep_e2 = zeros(18, 0);


 
for iii = 3: 19
    
 f =  fullfile('C:/Users/80677/Desktop/8point-algorithm-master/coresponding/',temp(iii).name,'/');
 path_png = dir([f,'*.txt']);
path_png1 =  fullfile(f,path_png(1).name);
path_png2 =  fullfile(f,path_png(2).name);
p1=load(path_png1);
p2=load(path_png2);
s = length(p1);
kpt1 = [p1, ones(s,1)];
kpt2 = [p2, ones(s,1)];
kpt1 = kpt1';
kpt2 = kpt2';

 
    
% compute true rotation and true translation
path_mat = dir([f,'*.camera']);
path_mat1 =  fullfile(f,path_mat(1).name);
path_mat2 =  fullfile(f,path_mat(2).name);
data1 = load(path_mat1);
K = data1(1:3, :); %all cameras share the same intrinsic matrix k
R1 = data1(5:7, :);
t1 = data1(8, :);
data2 = load(path_mat2);
R2 = data2(5:7, :);
t2 = data2(8, :);
% 
Rt = inv(R1) * R2;
Tt = t2 - t1;
Tt = Tt';
Tt = Tt/norm(Tt);


%Essential matrix estimation
[E1, x, y1, y2] = TwoViewGeometry.EMatrix (kpt1, kpt2);
E2 = essentialMatrix(kpt1, kpt2);



%four possible options
[R1, t1,~] = decompose(E1, K,kpt1, kpt2);
[R2, t2, ~] = decompose(E2, K,kpt1, kpt2);

kpt1 = kpt1(1:2,:);
kpt2 = kpt2(1:2,:);
P1 = K * [ eye(3),  zeros(3, 1)];
P2 = K * [R1, t1];
[pts3d1, error1] = triangulate(kpt1', kpt2', P1', P2');
P1 = K * [ eye(3),  zeros(3, 1)];
P2 = K * [R2, t2];
[pts3d2, error2] = triangulate(kpt1', kpt2', P1', P2');
fprintf(1, 'error1 error2 = %f  %f\n', [error1,error2]);



% reprojection error
% error1 = reprojection(kpt1, kpt2, K, pts3d1, R1, t1);
% error2 = reprojection(kpt1, kpt2, K, pts3d2, R2, t2);
% fprintf(1, 'error1 error2 = %f  %f\n', [error1,error2]);



% distance error
% kpt1 = inv(K)* kpt1 ;
% kpt2 = inv(k)* kpt2 ;
% [dis_er1] = diserror(E1,kpt1, kpt2);
% [dis_er2] = diserror(E2,kpt1, kpt2);
% fprintf(1, 'dis_er1 dis_er2 = %f  %f\n', [dis_er1,dis_er2]);



end



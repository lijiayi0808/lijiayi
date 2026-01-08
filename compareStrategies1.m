% 该函数是用来比较高斯算法过程中是否达到均衡用的
% A和B为两个等长x行y列的数列
% cS 比较结果，cS为1时，两个序列不相同，cS为0时，两个序列相同
function cS = compareStrategies1(A,B)

% 进行比较
C=(A-B);
[x,y]=size(C);
cSc = zeros(x,y);

for i=1:x
    for ii=1:y
        if abs(C(i,ii))<= 0.01
            cSc(i,ii) = 0;
        else
            cSc(i,ii) = 1;
        end
    end
end

if cSc == 0
    cS = 0;
else
    cS = 1;
end

end




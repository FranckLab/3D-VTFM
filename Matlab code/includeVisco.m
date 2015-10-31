function [Tij] = includeVisco(Fij, Gneq, eda, dt,J)
timeTotal = length(Fij);
FijV = cell(timeTotal,1);
Tij = cell(timeTotal,1);
for i = 1:3
    for j = 1:3
    Tij{1}{i,j} = zeros(size(Fij{1}{1,1}));
    end
end
sizeV = size(Fij{1}{1,1})
FijV{1}{1,1} = ones(sizeV);
FijV{1}{1,2} = zeros(sizeV);
FijV{1}{1,3} = zeros(sizeV);
FijV{1}{2,1} = zeros(sizeV);
FijV{1}{2,2} = ones(sizeV);
FijV{1}{2,3} = zeros(sizeV);
FijV{1}{3,1} = zeros(sizeV);
FijV{1}{3,2} = zeros(sizeV);
FijV{1}{3,3} = ones(sizeV);

for time = 1:timeTotal-1
    % Compute the trial elastic deformation gradient
    FijVinv = funInvMatrix(FijV{time});
    FijE_trial = cellfun(@(x,y) x.*y, Fij{time+1}, FijVinv,'UniformOutput',false);
    Uij = funCalculateUij(Fij{time+1});
    [Rij,~] = funCalculateRij(Fij{time+1}, Uij);
    CijE_trial = funCalcualteCij(FijE_trial);
    EijE_trial_ = funCalculateLN(CijE_trial);
    % Compute trial kinematics
    EijE_trial = cellfun(@(x) x*0.5, EijE_trial_, 'UniformOutput',false);
    trace = EijE_trial{1,1}+EijE_trial{2,2}+EijE_trial{3,3};
    Eij0 = cell(3,3);
    for i = 1:3
        for j = 1:3
            dij = i==j;
            Eij0{i,j} = EijE_trial{i,j}-1/3 * trace*dij;
        end
    end
    % Compute the trial Mandel Stress, which is deviatoric
    MtrE = cellfun(@(x) x*2*Gneq, Eij0, 'UniformOutput',false);
    % Compute the trial equilibrium tensile stress
    TauTr = funCalculateTau(MtrE);
    % Compute direction of viscous flow
    N_v = cellfun(@(x) 1/sqrt(2)*x./TauTr, MtrE, 'UniformOutput',false);
    for jj = 1:9
       N_v_ = N_v{jj};
       if sum(isnan(N_v_(:)))>0
           N_v_(isnan(N_v_)) = 0;
       end
       N_v{jj} = N_v_;
    end
    % Compute the equivalent shear viscous strain rate
    Gamma_v = TauTr/(eda+Gneq*dt);
    % Compute the viscous stretching
    D_v = cellfun(@(x) 1/(sqrt(2)) * Gamma_v.*x, N_v, 'UniformOutput',false);
    %Viscous deformation gradient
    FijV{time+1}= funCalculateExp(dt, D_v);
    FijV{time+1} = cellfun(@(x,y) x.*y, FijV{time+1}, FijV{time}, 'UniformOutput',false);
    %Stress
    MnUpdate = cellfun(@(x,y) x-2*Gneq*dt*y, MtrE, D_v,'UniformOutput',false);
    Tn_ = funUpdateMn(MnUpdate,Rij);
    Tn = cellfun(@(x) x./J{time}, Tn_,'UniformOutput',false);
    Tij{time+1} = Tn;
end
end
%% Functions
function cij = funCalcualteCij(Fij)
cij = num2cell(zeros(3));
for i = 1:3
    for j = 1:3
        for k = 1:3, cij{i,j} = cij{i,j} + Fij{k,i}.*Fij{k,j}; end
    end
end
end


function lnTensorCell = funCalculateLN(Tensor)
sizeTensor = size(Tensor{1});
for k = 1:9
    Tensor{k} = reshape(Tensor{k},1,1,[]);
end

Tensor_ = cell2mat(Tensor);

n = size(Tensor_,3);
D = zeros(n,3);
lnD_ = zeros(3,3);
lnTensor = zeros(size(Tensor_));
for i = 1:n
    [V,D_] = eig(Tensor_(:,:,i));
    D(i,1) = D_(1); D(i,2) = D_(5); D(i,3) = D_(9);
    lnD_(1) = log(D_(1));lnD_(5) = log(D_(5));lnD_(9) = log(D_(9));
    lnTensor(:,:,i) = V*lnD_*V';
end



for i = 1:3
    for j = 1:3
        lnTensorCell{i,j} = reshape(lnTensor(i,j,:), sizeTensor);
    end
end
end

function TauTr = funCalculateTau(Tensor)

sizeTensor = size(Tensor{1});
for k = 1:9
    Tensor{k} = reshape(Tensor{k},1,1,[]);
end

Tensor_ = cell2mat(Tensor);
n = size(Tensor_,3);
InnerProduct = zeros(1,n);
for i = 1:n
    T = Tensor_(:,:,i);
    InnerProduct(:,i) = T(1,1)*T(1,1)+T(1,2)*T(1,2) + T(1,3)*T(1,3) + ...
        T(2,1)*T(2,1) + T(2,2)*T(2,2)+T(2,3)*T(2,3)+...
        T(3,1)*T(3,1) + T(3,2)*T(3,2)+T(3,3)*T(3,3);
end

TauTr_ = 1/(sqrt(2))*sqrt(InnerProduct);
TauTr = reshape(TauTr_, sizeTensor);
end

function expTensorCell = funCalculateExp(dt,Tensor)
sizeTensor = size(Tensor{1});
for k = 1:9
    Tensor{k} = reshape(Tensor{k},1,1,[]);
end

Tensor_ = cell2mat(Tensor);

n = size(Tensor_,3);
D = zeros(n,3);
expD_ = zeros(3,3);
expTensor = zeros(size(Tensor_));
for i = 1:n
    [V,D_] = eig(Tensor_(:,:,i));
    D(i,1) = D_(1); D(i,2) = D_(5); D(i,3) = D_(9);
    expD_(1) = exp(dt*D_(1));expD_(5) = exp(dt*D_(5));expD_(9) = exp(dt*D_(9));
    expTensor(:,:,i) = V*expD_*V';
end



for i = 1:3
    for j = 1:3
        expTensorCell{i,j} = reshape(expTensor(i,j,:), sizeTensor);
    end
end
end

function Uij = funCalculateUij(Fij)
Cij = cell(3,3);
Uij = cell(3,3);
Umat = cell(length(Fij{3,3}(:)),1);

Cij{1,1} = (Fij{1,1}.*Fij{1,1} + Fij{2,1}.*Fij{2,1} + Fij{3,1}.*Fij{3,1});
Cij{1,2} = (Fij{1,1}.*Fij{1,2} + Fij{2,1}.*Fij{2,2} + Fij{3,1}.*Fij{3,2});
Cij{1,3} = (Fij{1,1}.*Fij{1,3} + Fij{2,1}.*Fij{2,3} + Fij{3,1}.*Fij{3,3});
Cij{2,2} = (Fij{1,2}.*Fij{1,2} + Fij{2,2}.*Fij{2,2} + Fij{3,2}.*Fij{3,2});
Cij{2,3} = (Fij{1,2}.*Fij{1,3} + Fij{2,2}.*Fij{2,3} + Fij{3,2}.*Fij{3,3});
Cij{3,3} = (Fij{1,3}.*Fij{1,3} + Fij{2,3}.*Fij{2,3} + Fij{3,3}.*Fij{3,3});

for j=1:length(Fij{3,3}(:));
    Umat{j} = sqrtm([Cij{1,1}(j), Cij{1,2}(j), Cij{1,3}(j); Cij{1,2}(j), Cij{2,2}(j), Cij{2,3}(j); Cij{1,3}(j), Cij{2,3}(j), Cij{3,3}(j)]);
    Uij{1,1}{j} = Umat{j}(1,1);
    Uij{1,2}{j} = Umat{j}(1,2);
    Uij{1,3}{j} = Umat{j}(1,3);
    Uij{2,2}{j} = Umat{j}(2,2);
    Uij{2,3}{j} = Umat{j}(2,3);
    Uij{3,3}{j} = Umat{j}(3,3);
    
end

Uij{1,1} = reshape(cell2mat(Uij{1,1}),size(Fij{1,1}));
Uij{1,2} = reshape(cell2mat(Uij{1,2}),size(Fij{1,1}));
Uij{1,3} = reshape(cell2mat(Uij{1,3}),size(Fij{1,1}));
Uij{2,2} = reshape(cell2mat(Uij{2,2}),size(Fij{1,1}));
Uij{2,3} = reshape(cell2mat(Uij{2,3}),size(Fij{1,1}));
Uij{3,3} = reshape(cell2mat(Uij{3,3}),size(Fij{1,1}));
end

function [Rij,angles] = funCalculateRij(Fij, Uij)
Rij = cell(3,3);
Umat = cell(length(Fij{3,3}(:)),1);
Rmat = cell(length(Fij{3,3}(:)),1);
Fmat = cell(length(Fij{3,3}(:)),1);
angles = cell(3,1);


for j=1:length(Fij{3,3}(:));
    Umat{j} = [Uij{1,1}(j), Uij{1,2}(j), Uij{1,3}(j); Uij{1,2}(j), Uij{2,2}(j), Uij{2,3}(j); Uij{1,3}(j), Uij{2,3}(j), Uij{3,3}(j)];
    Fmat{j} = [Fij{1,1}(j), Fij{1,2}(j), Fij{1,3}(j); Fij{2,1}(j), Fij{2,2}(j), Fij{2,3}(j); Fij{3,1}(j), Fij{3,2}(j), Fij{3,3}(j)];
    Rmat{j} = Fmat{j}/Umat{j};
    Rij{1,1}{j} = Rmat{j}(1,1);
    Rij{1,2}{j} = Rmat{j}(1,2);
    Rij{1,3}{j} = Rmat{j}(1,3);
    Rij{2,1}{j} = Rmat{j}(2,1);
    Rij{2,2}{j} = Rmat{j}(2,2);
    Rij{2,3}{j} = Rmat{j}(2,3);
    Rij{3,1}{j} = Rmat{j}(3,1);
    Rij{3,2}{j} = Rmat{j}(3,2);
    Rij{3,3}{j} = Rmat{j}(3,3);
    
    angles{1}{j}=atan2(Rij{3,2}{j},Rij{3,3}{j})*180/pi; %range: -pi to pi
    angles{2}{j}=atan2(-Rij{3,1}{j},sqrt(Rij{3,2}{j}.^2+Rij{3,3}{j}.^2))*180/pi; %range: -pi/2 to pi/2
    angles{3}{j}=atan2(Rij{2,1}{j},Rij{1,1}{j})*180/pi; %range: -pi to pi
end

Rij{1,1} = reshape(cell2mat(Rij{1,1}),size(Fij{1,1}));
Rij{1,2} = reshape(cell2mat(Rij{1,2}),size(Fij{1,1}));
Rij{1,3} = reshape(cell2mat(Rij{1,3}),size(Fij{1,1}));
Rij{2,1} = reshape(cell2mat(Rij{2,2}),size(Fij{1,1}));
Rij{2,2} = reshape(cell2mat(Rij{2,2}),size(Fij{1,1}));
Rij{2,3} = reshape(cell2mat(Rij{2,3}),size(Fij{1,1}));
Rij{3,1} = reshape(cell2mat(Rij{3,1}),size(Fij{1,1}));
Rij{3,2} = reshape(cell2mat(Rij{3,2}),size(Fij{1,1}));
Rij{3,3} = reshape(cell2mat(Rij{3,3}),size(Fij{1,1}));
angles{1} = reshape(cell2mat(angles{1}),size(Fij{1,1}));
angles{2} = reshape(cell2mat(angles{2}),size(Fij{1,1}));
angles{3} = reshape(cell2mat(angles{3}),size(Fij{1,1}));
end


function Tn = funUpdateMn(MnUpdate,Rij)

sizeTensor = size(MnUpdate{1});
for i = 1:3
    for j = 1:3
        Tensor{i,j} = reshape(MnUpdate{i,j},1,1,[]);
        TensorR{i,j} = reshape(Rij{i,j},1,1,[]);
    end
end
Tensor_ = cell2mat(Tensor);
TensorR_ = cell2mat(TensorR);
n = size(Tensor_,3);

Tn_ = zeros(size(Tensor_));
for i = 1:n
    Tn_(:,:,i) = TensorR_(:,:,i)*Tensor_(:,:,i)*TensorR_(:,:,i)';
end

for i = 1:3
    for j = 1:3
        Tn{i,j} = reshape(Tn_(i,j,:), sizeTensor);
    end
end

end

function InverseA = funInvMatrix(TensorD)

sizeTensor = size(TensorD{1});
for i = 1:3
    for j = 1:3
        Tensor{i,j} = reshape(TensorD{i,j},1,1,[]);
    end
end
ATensor = cell2mat(Tensor);
n = size(ATensor,3);

A_inv = zeros(size(ATensor));
for i = 1:n
    A = ATensor(:,:,i);
    det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) - ...
        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +...
        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3));
    
    det_A_inv = 1/abs(det_A);
    
    A_inv_(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3));
    A_inv_(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3));
    A_inv_(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3));
    A_inv_(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3));
    A_inv_(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3));
    A_inv_(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3));
    A_inv_(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2));
    A_inv_(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2));
    A_inv_(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2));
    
    A_inv(:,:,i) = A_inv_;
end

for i = 1:3
    for j = 1:3
        InverseA{i,j} = reshape(A_inv(i,j,:), sizeTensor);
    end
end
end

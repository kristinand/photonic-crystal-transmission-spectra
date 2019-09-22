% Input conditions
eps1 = 11.9;
eps2 = 5.8;
b1 = 1028; % Thickness
b2 = 1472;
b12 = 1250;
nm = 10^(-9);
c = 3*10^8;
N = 20; % Bilayer number
P = 2; % Period
l = 10000:21000; % wavelength
w = 2*pi*c./(l*nm);
lam0 = 14160; %nm
w0 = 2*pi*c/(lam0*nm);

% Refractive index array
n1 = sqrt(eps1);
n2 = sqrt(eps2);
n= [n1, n2];
nArray = [1, repmat(n,1,N), 1]; % without defect, where 1 is in air

% Thickness array
b = zeros(1,N*P); 
for i = 1:2:N*P
    b(i) = b1;
end

for i = 2:2:N*P
    b(i) = b2;
end

% Wave vector in vacuum
k=2*pi./l;                          

% Fresnel coefficients
r=zeros(1, N*P+1);         
t=zeros(1, N*P+1);

for i=1:1:N*P+1
    r(i)=(nArray(i)-nArray(i+1))/(nArray(i)+nArray(i+1));
    t(i)=2*nArray(i)/(nArray(i)+nArray(i+1));
end

% Transfer and transmit matricies in layer j
M=zeros(2,2,N*P+1);
F=zeros(2,2,length(l),N*P);

for i=1:1:N*P+1
    M(1,1,i)=1/t(i);
    M(1,2,i)=r(i)/t(i);
    M(2,1,i)=r(i)/t(i);
    M(2,2,i)=1/t(i);
end

for i=1:1:N*P
    for j=1:length(l)
    F(1,1,j,i)=exp(1i*nArray(i+1)*k(j)*b(i));
    F(1,2,j,i)=0;
    F(2,1,j,i)=0;
    F(2,2,j,i)=exp(-1i*nArray(i+1)*k(j)*b(i));
    end
end

% Full matrix
T=zeros(2,2,length(l));
for j=1:length(l)
    T(:,:,j)=M(:,:,1);
    for i=1:1:N*P
        T(:,:,j)=T(:,:,j)*F(:,:,j,i)*M(:,:,i+1);
    end
end

M=zeros(length(l));
for j=1:length(l)
    M(j)=abs( (T(1,1,j).*T(2,2,j)-T(1,2,j).*T(2,1,j))./T(2,2,j) ).^2;
end

figure(1)
plot(w/w0,M(:,1),'k','linewidth',1.5);
axis([0.7 1.3 0 1]);
grid on;
yticks([0,0.25,0.5,0.75,1]);
set(gca,'FontSize',16,'FontName', 'Times New Roman');
hold on;
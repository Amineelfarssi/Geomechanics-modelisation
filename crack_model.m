tic;
%%%%%%%%%%%%%%%%%%%% input parameters
x_beg=1/4; x_end=3/4; % fault limits in [0, 1], in units of lambda
time=1;             % total time of simulation, in unit of lambda/c
p=6;                  % 2^p points to discretize the replication period lambda
h=0.5;                % "Courant" parameter c*dt/dx
%%%%%%%%%%%%%%%%%%%% end inpt 

c=1; mu=1; tau0=1; friction=0; lambda=1;  % can be changed, but no point
nx=2^p; dx=lambda/nx; dt=h*dx/c; nt=floor(time/dt); 

i_beg=floor(x_beg/lambda/dx+1)+1; i_end=floor(x_end/lambda/dx)-1;     % fault limits

% save('variables_p6.mat','dt','nt','nx','dx','i_beg','i_end','time')

slip=zeros(nx,1);
 velo=zeros(nx,1);
 func=zeros(nx,1);
slipfft=zeros(nx,1);
funcfft=zeros(nx,1);
histfft=zeros(nx/2+1,nt);
velohist=zeros(nx,nt+1);


i_beg=(nx/4+1)+1; i_end=(3*nx/4+1)-1; % fault limits
a=2*pi*c*dt/lambda; b=-2*mu*c*pi^2/lambda^2; % some constants
% Probably better to use an (external) M-file for kern, for efficiency;
% kern(t)=J1(t)/t, which -> 0.5 when t->0 (J1 1st kind bessel)
kern=@(t)real([t(t==0)+0.5 besselj(1,t(t>eps))/t(t>eps) t(t<eps & t~=0 )/t(t<eps & t~=0)-0.5]);

for n=0:nt % time loop
    % Euler scheme
    for i=i_beg:i_end
        slip(i)=slip(i)+velo(i)*dt;
    end
    slipfft=fft(slip);% DFT of slip and save.
    histfft(1:nx/2+1,n+1)=slipfft(1:nx/2+1);

    % Computation of the functional in the Fourier domain
    funcfft(1:nx)=0;
    for j=1:nx/2+1
        % Time integration for each Fourier mode (0->nx/2, i.e., j: 1->nx/2+1),
        % using the "preintegrated" kernel version (w/  trapezoidal rule).
        for m=1:n
            funcfft(j)=funcfft(j) + ...
            histfft(j,m+1)*dt*( kern(a*(j-1)*(n-m+1))+kern(a*(j-1)*(n-m)) )/2 ;
        end
        funcfft(j)=b*(j-1)^2*funcfft(j);
    end
    for j=nx/2+2:nx;
        funcfft(j)=conj(funcfft(nx+2-j));
    end % symetrization
    func=real(ifft(funcfft)); % back in the spatial domain
    %Solve for the velocity using the BCs on the fault (friction)
    for i=i_beg:i_end
        velo(i)=(2*c)/mu * (tau0+func(i)-friction);
        velohist(i,n+1)=velo(i); % save just for plotting
    end
end

%  save('Test_p6.mat','velohist','velo')
plot(velo);
 hold on
mesh(real(velohist));
xlabel('time')
ylabel('Fault position')
zlabel('Velocity')
figure;
plot(velohist(:,floor(1/10/dt)));
temps=toc 

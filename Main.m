function Ans = Main()



I = sqrt(-1);
L = 60;
Ngrid = 512;
dx = L/Ngrid;
[x,y] = meshgrid(-L/2:dx:L/2);
[px,py] = meshgrid(-pi*Ngrid/L:2*pi/L:pi*Ngrid/L);
[theta,r] = cart2pol(x,y);
m1 = 1;
m2 = -1;
lambda1 = 0.1;
lambda2 = 0.1;
alpha = 1;
lim = 0.00000000001;
num = 50;
st = 0.2;
NumIter = 50;


NgridR = Ngrid*10;
dxR = L/(2*NgridR);
ro = linspace(0,L/2,NgridR);

positionVector1 = [0.05, 0.05, 0.23, 0.4];
positionVector2 = [0.325, 0.05, 0.135, 0.4];
positionVector3 = [0.05, 0.55, 0.23, 0.4];
positionVector4 = [0.325, 0.55, 0.135, 0.4];
positionVector5 = [0.475, 0.05, 0.422, 0.9];



dob = 0*x.^2;
dob(1,1) = 0.08;
dobt = 0*x.^2;
%dobt(1,1) = 1.4;

numticks = 6;
%cticks = linspace(log10(0.25),0,numticks);
%clabels = round(((10.^(cticks))*4-1)*100)/100;
clabels = linspace(0,0.08,numticks);
cticks = clabels;
numtickst = 8;
%ctickst = linspace(log10(1/4.5),0,numtickst);
%clabelst = round(((10.^(ctickst))*4.5-1)*100)/100;
clabelst = linspace(0.6,2,numtickst);
ctickst = clabelst;

%clabelsnum = num2str(clabels);
disp(clabels);
disp(clabelst);

%for alpha=1:-0.1:0.5
    %alpha = 0.6;


    
%[psi1,psi2] = CalcGrSt();
%eta = CalcThetaR(NgridR,dxR,psi1,psi2,ro);


%plot(ro,psi1,ro,eta,'--');
%xlabel('r');
%legend('\psi_n','\theta');
%saveas(gcf,'conf1.jpg');


%----------------------------------------------
%main part

    [Phi1,Phi2] = Initiate();
    Theta = CalcTheta(Phi1,Phi2);
    disp('Ground state calculated');
    
    
    
    [Phi1,Phi2,Theta]=Dynamics(Phi1,Phi2,Theta);

    
%-----------------------------------------------
%disp(Energy(Phi1,Phi2,Theta));
    
%end



%mesh(x,y,abs(Phi1));

%pcolor(x,y,abs(Phi1));
%shading interp;
%colormap jet;


    function [Phi1R,Phi2R] = CalcGrSt()
        
        
        Phi1R = (ro.^(abs(m1))).*exp(-ro);
        Phi2R = (ro.^(abs(m2))).*exp(-ro);
        
        
        du = 1;
        
        a = -dxR^(-2) + ro.^(-1)*dxR^(-1)*0.5;
        a(1) = 0;
        b1 = lambda1 + 2*dxR^(-2) + m1^2*ro.^(-2);
        b2 = lambda2 + 2*dxR^(-2) + m2^2*ro.^(-2);
        c = -dxR^(-2) - ro.^(-1)*dxR^(-1)*0.5;
        c(1) = -2*dxR^(-2);
        
        
        while (abs(du)>lim)
            
            ThetaR = CalcThetaR(NgridR,dxR,Phi1R,Phi2R,ro);
            
            d1 = ThetaR.*Phi1R;
            d2 = ThetaR.*Phi2R;
            
            NormP1 = CalcN(Phi1R);
            NormP2 = CalcN(Phi2R);
            
            Phi1R = TDMAsolver(a,b1,c,d1);
            Phi2R = TDMAsolver(a,b2,c,d2);
            
            Norm1 = CalcN(Phi1R);
            Norm2 = CalcN(Phi2R);
            S1 = (NormP1/Norm1)^(0.7);
            S2 = (NormP2/Norm2)^(0.7);
            
            Phi1R=Phi1R*S1;
            Phi2R=Phi2R*S2;
            
            
            du = max(abs(S1-1),abs(S2-1));
            
        end
        
        
        
        
    end    
    function ThetaR = CalcThetaR(NgridR,dxR,Phi1R,Phi2R,ro)
        
        a = -dxR^(-2) + ro.^(-1)*dxR^(-1)*0.5;
        b = alpha*ro.^0 + 2*dxR^(-2);
        c = -dxR^(-2) - ro.^(-1)*dxR^(-1)*0.5;
        d = Phi1R.^2 + Phi2R.^2;
        c(1) = -2*dxR^(-2);
        ThetaR = TDMAsolver(a,b,c,d);
    
    end    
    function NormR = CalcN(Phi)
        NormR = sum(Phi.*Phi);
    end
    function Dev = CalcDev(Phi,Theta,lambda,m)
        
        n = length(Phi);
        
        Dev = 0;
        
        a = -dxR^(-2) + ro.^(-1)*dxR^(-1)*0.5;
        a(1) = 0;
        b = lambda + 2*dxR^(-2) + m^2*ro.^(-2);
        c = -dxR^(-2) - ro.^(-1)*dxR^(-1)*0.5;
        c(1) = -2*dxR^(-2);
        d = Theta.*Phi;
        
        for in = 2:n-1
            k = a(in)*Phi(in-1)+b(in)*Phi(in)+c(in)*Phi(in+1)-d(in);
            Dev = Dev +k;
        end
            
    
    
    end
    function phi = Transform(PhiR)
        
        phi = exp(-r);
        
        
        for in=1:Ngrid
            for j=1:Ngrid
                [h,n]=min(abs(ro-r(in,j)));
                if (n == 1)
                    phi(in,j) = PhiR(1);
                else
                    if (n == NgridR)
                        phi(in,j) = 0;
                    else
                  
                        %phi(i,j) = d^2*a+d*b+c;
                        phi(in,j) = PhiR(n);
                    end
                end
            end
        end
                
        
        
        
        
    end
    function [f1,f2]=Initiate()
    
        [Phi1R,Phi2R] = CalcGrSt();
        
        ThetaR = CalcThetaR(NgridR,dxR,Phi1R,Phi2R,ro);
        
        disp(CalcDev(Phi1R,ThetaR,lambda1,m1));
        disp(CalcDev(Phi2R,ThetaR,lambda2,m2));
        
        %plot(ro,ThetaR);
        
        f1 = Transform(Phi1R);
        f2 = Transform(Phi2R);
        f1 = f1.*exp(I*m1*theta);
        f2 = f2.*exp(I*m2*theta);
    
    end
    function t = CalcTheta(phi1,phi2)
        
        t = abs(ifft2(ifftshift(((alpha^2+px.^2+py.^2).^(-1)).*fftshift(fft2(phi1.*conj(phi1)+phi2.*conj(phi2))))));
        
        
    end
    function psi = LinearEvol(phi,tau)
        psi = ifft2(ifftshift(exp(-I*(px.^2+py.^2)*tau).*fftshift(fft2(phi))));
    end
    function psi = NonlinearEvol(phi,tau,th)
        psi = exp(I*th*tau).*phi;
    end
    function [phi1,phi2,th] = Evolve(phi1,phi2,th,Num,Step)
    
        th = CalcTheta(phi1,phi2);
        
        for i=1:Num
            phi1 = LinearEvol(phi1,Step/2);
            phi1 = NonlinearEvol(phi1,Step,th);
            phi1 = LinearEvol(phi1,Step/2);
            
            phi2 = LinearEvol(phi2,Step/2);
            phi2 = NonlinearEvol(phi2,Step,th);
            phi2 = LinearEvol(phi2,Step/2);
            
            th = CalcTheta(phi1,phi2);
        end
    
    end
    function buf = Norm2D(phi)
        buf = sum(sum(abs(phi.*conj(phi))))*dx*dx;
    end
    function E = EnergyU(phi1,phi2,th)

        E = abs(sum(sum(-th.*(phi1.*conj(phi1)+phi2.*conj(phi2))))*0.5);
    end
    function E = EnergyK(phi1,phi2,th)
        E = abs(sum(sum(conj(phi1).*Lapl(phi1)+conj(phi2).*Lapl(phi2))));

    end
    function E = Energy(phi1,phi2,th)
        E = EnergyU(phi1,phi2,th)-EnergyK(phi1,phi2,th);
    end
    function delt = Lapl(phi)
        delt = ifft2(ifftshift((px.^2+py.^2).*fftshift(fft2(phi))));
    end
    function sr = AvRad(phi)
        sr = abs(sum(sum(sqrt((x.^2+y.^2)).*phi.*conj(phi)))*dx*dx/Norm2D(phi));
    end
    function sr = AvRadT(th)
        sr = abs(sum(sum(sqrt(x.^2+y.^2).*th)))/sum(sum(th));
    end
    function Px = ImpulseX(phi)
        buf = conj(phi).*ifft2(ifftshift(-I*px.*fftshift(fft2(phi)))) - phi.*ifft2(ifftshift(-I*px.*fftshift(fft2(conj(phi)))));
        Px = abs(sum(sum(buf))*dx*dx*I*0.5);
        
    end
    function Py = ImpulseY(phi)
        buf = conj(phi).*ifft2(ifftshift(-I*py.*fftshift(fft2(phi)))) - phi.*ifft2(ifftshift(-I*py.*fftshift(fft2(conj(phi)))));
        Py = abs(sum(sum(buf))*dx*dx*I*0.5);
    end
    function M = Moment(phi)
        buf1 = conj(phi).*ifft2(ifftshift(-I*px.*fftshift(fft2(phi)))) - phi.*ifft2(ifftshift(-I*px.*fftshift(fft2(conj(phi)))));
        buf2 = conj(phi).*ifft2(ifftshift(-I*py.*fftshift(fft2(phi)))) - phi.*ifft2(ifftshift(-I*py.*fftshift(fft2(conj(phi)))));
        M = abs(sum(sum(x.*buf2-y.*buf1))*I*dx*dx*0.5);
    end
    function [Phi1,Phi2,Theta] = Dynamics(Phi1,Phi2,Theta)
        
        name = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\alpha' num2str(alpha)];
        namep1r = [name '\Phi1R'];
        namep1a = [name '\Phi1A'];
        namep2r = [name '\Phi2R'];
        namep2a = [name '\Phi2A'];
        namet = [name '\Theta'];
        nameo = [name '\overall'];
        
        
        mkdir(name);
        mkdir(namep1r);
        mkdir(namep1a);
        mkdir(namep2r);
        mkdir(namep2a);
        mkdir(namet);
        mkdir(nameo);
        
        
        
        
        
        Rstart = AvRad(Phi1);
        in = 0;
        
        for in=1:300
        %while (Rstart/AvRad(Phi1)>0.7)
            %in = in + 1;
            
         
            [Phi1,Phi2,Theta] = Evolve(Phi1,Phi2,Theta,num,st);
            EK(in) = EnergyK(Phi1,Phi2,Theta);
            EU(in) = EnergyU(Phi1,Phi2,Theta);
            T(in) = num*st*in;
            R(in) = AvRad(Phi1);
            Rt(in) = AvRadT(Theta);
            M(in) = Moment(Phi1);
            M2(in) = Moment(Phi2);
            E(in) = Energy(Phi1,Phi2,Theta);
            N1(in) = Norm2D(Phi1);
            N2(in) = Norm2D(Phi2);
            I1x(in) = ImpulseX(Phi1);
            I1y(in) = ImpulseY(Phi1);
            I2x(in) = ImpulseX(Phi2);
            I2y(in) = ImpulseY(Phi2);
            
            namef = ['\' num2str(T(in)) '.jpg'];
            
            figure
            pcolor(x,y,(abs(Phi1).^2));
            shading interp;
            colormap jet;
            filename = [namep1r namef];
            set(gcf, 'PaperPosition', [0 0 15 15]);
            saveas(gcf,filename);
    
            pcolor(x,y,angle(Phi1));
            shading interp;
            colormap jet;
            filename = [namep1a namef];
            set(gcf, 'PaperPosition', [0 0 15 15]);
            saveas(gcf,filename);
            
            pcolor(x,y,(abs(Phi2).^2));
            shading interp;
            colormap jet;
            filename = [namep2r namef];
            set(gcf, 'PaperPosition', [0 0 15 15]);
            saveas(gcf,filename);
    
            pcolor(x,y,angle(Phi2));
            shading interp;
            colormap jet;
            filename = [namep2a namef];
            set(gcf, 'PaperPosition', [0 0 15 15]);
            saveas(gcf,filename);
            
            
    
            pcolor(x,y,Theta);
            shading interp;
            colormap jet;
            filename = [namet namef];
            set(gcf, 'PaperPosition', [0 0 15 15]);
            saveas(gcf,filename);
            close(gcf);
            
            figure
            subplot('Position',positionVector3);
            %pcolor(x,y,log10((abs(Phi1).^2+dob)/1.8));
            pcolor(x,y,(abs(Phi1).^2)+dob);
            %pcolor(x,y,(abs(Phi1).^2));
            axis([-15;15;-15;15]);
            %xlabel('x');
            %ylabel('y');
            shading interp;
            colormap jet;
            %colorbar('location','westoutside','YTick',cticks,'YTickLabel',clabels);
            %colorbar();
            colorbar('location','westoutside');
            title('Phi1^2')
            set(gca,'xtick',[-15 0 15]);
            %set(gca,'xticklabel',[]);
            set(gca,'ytick',[-15 0 15]);
            %set(gca,'yticklabel',[]);


            subplot('Position',positionVector4);
            %pcolor(x,y,log10((abs(Phi2).^2+dob)/1.8));
            pcolor(x,y,(abs(Phi2).^2)+dob);
            %pcolor(x,y,(abs(Phi2).^2));
            axis([-15;15;-15;15]);
            %xlabel('x');
            %ylabel('y');
            shading interp;
            colormap jet;
            title('Phi2^2')
            set(gca,'xtick',[-15 0 15]);
            %set(gca,'xticklabel',[]);
            set(gca,'ytick',[-15 0 15]);
            %set(gca,'yticklabel',[]);

            %Theta(1,1) = 2;
            subplot('Position',positionVector5);
            %pcolor(x,y,log10((Theta+dobt)/3));
            pcolor(x,y,Theta);
            axis([-15;15;-15;15]);
            %xlabel('x');
            %ylabel('y');
            shading interp;
            colormap jet;
            %colorbar('location','westoutside','YTick',ctickst,'YTickLabel',clabelst);
            colorbar('location','westoutside');
            title_name = ['Theta,   z=' num2str(T(in))];
            title(title_name)
            set(gca,'xtick',[-15 0 15]);
            %set(gca,'xticklabel',[]);
            set(gca,'ytick',[-15 0 15]);
            %set(gca,'yticklabel',[]);


            subplot('Position',positionVector1);
            pcolor(x,y,angle(Phi1));
            shading interp;
            colormap jet;
            colorbar('location','westoutside');
            title('angle(Phi1)')
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);


            subplot('Position',positionVector2);
            pcolor(x,y,angle(Phi2));
            shading interp;
            colormap jet;
            title('angle(Phi2)')
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);


            filename = [nameo namef];
            
            set(gcf, 'PaperPosition', [0 0 45 15]);
            
            
            saveas(gcf,filename);
            
            close(gcf);
    
            disp(T(in));
            
    
    
        end
        
        
        namegr = [name '\graph'];
        mkdir(namegr);
        
        figure
        fig = plot(T,EK);
        filename = [namegr '\EK.jpg'];
        saveas(fig,filename);

        fig = plot(T,EU);
        filename = [namegr '\EU.jpg'];
        saveas(fig,filename);

        fig = plot(T,R);
        filename = [namegr '\R.jpg'];
        saveas(fig,filename);

        fig = plot(T,Rt);
        filename = [namegr '\Rt.jpg'];
        saveas(fig,filename);

        fig = plot(T,M);
        filename = [namegr '\M1.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,E);
        filename = [namegr '\E.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,M2);
        filename = [namegr '\M2.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,N2);
        filename = [namegr '\N2.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,N1);
        filename = [namegr '\N1.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,I1x);
        filename = [namegr '\I1x.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,I1y);
        filename = [namegr '\I1y.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,I2x);
        filename = [namegr '\I2x.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,I2y);
        filename = [namegr '\I2y.jpg'];
        saveas(fig,filename);
        
        close(gcf);
    end

    function [Phi1,Phi2,Theta] = Dynamics2(Phi1,Phi2,Theta)
        for in=1:NumIter
            [Phi1,Phi2,Theta] = Evolve(Phi1,Phi2,Theta,num,st);
            EK(in) = EnergyK(Phi1,Phi2,Theta);
            EU(in) = EnergyU(Phi1,Phi2,Theta);
            T(in) = num*st*in;
            R(in) = AvRad(Phi1);
            Rt(in) = AvRadT(Theta);
            M(in) = Moment(Phi1);
            E(in) = Energy(Phi1,Phi2,Theta);
    
            fig=pcolor(x,y,abs(Phi1));
            shading interp;
            colormap jet;
            filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\Phi1R\im' num2str(T(in)) '.jpg'];
            saveas(fig,filename);
    
            fig=pcolor(x,y,angle(Phi1));
            shading interp;
            colormap jet;
            filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\Phi1A\im' num2str(T(in)) '.jpg'];
            saveas(fig,filename);
    
            fig=pcolor(x,y,Theta);
            shading interp;
            colormap jet;
            filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\Theta\im' num2str(T(in)) '.jpg'];
            saveas(fig,filename);
    
    
            disp(T(in));
    
    
        end


        fig = plot(T,EK);
        filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\graphs\EK.jpg'];
        saveas(fig,filename);

        fig = plot(T,EU);
        filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\graphs\EU.jpg'];
        saveas(fig,filename);

        fig = plot(T,R);
        filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\graphs\R.jpg'];
        saveas(fig,filename);

        fig = plot(T,Rt);
        filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\graphs\Rt.jpg'];
        saveas(fig,filename);

        fig = plot(T,M);
        filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\graphs\M.jpg'];
        saveas(fig,filename);
        
        fig = plot(T,E);
        filename = ['C:\Users\MichaelS\Documents\MATLAB\Bachelor\graphs\E.jpg'];
        saveas(fig,filename);
    end
end


%% Fig.7
clear all
Nt=16;
Nr=16;
Nt_Beam=12;
Nr_Beam=12;
P_bar=2;
N3=2;
N_bar=10;
alphaScale=5;
alphaScale_N=15;

Nrf=4;
G=32;
L=12;
Qp=10;
lambda=0.001;
d=lambda/2;
sigma_alpha2=1;
sigma_AS=3.5/180*pi;
sigman2=1;
epsilon=1e-6;
M=2*Nr_Beam*N_bar ;
M_bar=2*Nr_Beam*Nt_Beam;

T_v=[1,5,10,15,20,25];% the number of iterations
Nx=100;

%%
for ttt=1:6
    sum1=0;
    sum2=0;
    Tmax=T_v(ttt);
    
    for f=1:Nx
        a_A=rand(L,1)*2*pi;
        b_A=sigma_AS/sqrt(2);
        a_D=rand(L,1)*2*pi;
        b_D=sigma_AS/sqrt(2);
        
        AoA=zeros(L,Qp);
        AoD=zeros(L,Qp);
        
        aR_vector=zeros(Nr,L*Qp);
        aT_vector=zeros(Nt,L*Qp);
        
        H=zeros(Nr,Nt);
        
        for kk=1:Qp
            alpha_vector(:,kk)=sqrt(sigma_alpha2/2)*(randn(L,1)+1i*randn(L,1));%/sqrt(L*Qp);
            
            AoA(:,kk)=laplace(a_A,b_A,L,1);
            
            AoD(:,kk)=laplace(a_D,b_D,L,1);
        end
        
        for pp=1:L
            for j=1:Qp
                for k=1:Nr
                    aR_vector(k,1)=(1/sqrt(Nr))*exp((-1i)*2*pi/lambda*(k-1)*d*cos(AoA(pp,j)));
                end
                for x=1:Nt
                    aT_vector(x,1)=(1/sqrt(Nt))*exp((-1i)*2*pi/lambda*(x-1)*d*cos(AoD(pp,j)));
                end
                h_h=alpha_vector(pp,j)* aR_vector*aT_vector';
                H=H+h_h;
            end
        end
        H=H*sqrt(Nr*Nt/(L*Qp));
        vec_H=reshape(H,Nr*Nt,1);
        
        F_rf=dftmtx(Nt)/sqrt(Nt);
        W_rf=dftmtx(Nr)/sqrt(Nr);
        
        F_BB=eye(Nt);
        W_BB=eye(Nr);
        F_bb=F_BB(:,randperm(size(F_BB,2),Nt_Beam));
        W_bb=W_BB(:,randperm(size(W_BB,2),Nr_Beam));
        
        F=F_rf*F_bb;
        W=W_rf*W_bb;
        W_bar=mat2cell(W,ones(Nr/Nr,1)*Nr,ones(Nr_Beam/(Nr_Beam/(Nr_Beam/Nrf)),1)*Nr_Beam/(Nr_Beam/Nrf));
        a=blkdiag(W_bar{:});
        b=a';
        
        n_vector=zeros(Nr_Beam/Nrf*Nr,N_bar);
        for p=1:N_bar
            n_vector(:,p)=(rand(Nr_Beam/Nrf*Nr,1)+1i*rand(Nr_Beam/Nrf*Nr,1))/sqrt(2)*sqrt(sigman2);
        end
        
        Noise=b*n_vector;
        n_bar=reshape(Noise,Nr_Beam*N_bar,1);
        
        Q=kron((F_bb.')*(F_rf.'),(W_bb')*(W_rf'));
        
        P=10;
        X12=blkdiag(sqrt(P)*eye(Nrf),sqrt(P)*eye(Nrf));
        X3=zeros(Nrf,2);
        for xxx=1:2
            X3(:,xxx)=sqrt(2/(P/Nrf))*(randn(Nrf,1)+1i*randn(Nrf,1));
        end
        X=blkdiag(X12,X3);
        
        Y_bar_12=W'*H*F*X+Noise;
        y_bar_12=reshape(Y_bar_12,Nr_Beam*N_bar,1);
        
        Y_bar=W'*H*F*X(:,1:P_bar*Nrf)+Noise(:,1:P_bar*Nrf);
        y_bar=reshape(Y_bar,Nt_Beam*P_bar*Nrf,1);
        [Y_sort,index]=sort(abs(y_bar(:)),'descend');
        [row,col]=ind2sub(size(Y_bar),index);
        col_L=zeros(Nrf,1);
        row_L=zeros(Nrf,1);
        
        col_L(1)=col(1);
        row_L(1)=row(1);
        pp=2;
        mm=1;
        II=zeros(Nrf,1);
        II(1)=1;
        while(mm<=Nrf-1)
            for ii=pp:(Nr_Beam*Nt_Beam)
                RR=abs(col(ii)-col_L);
                GG=abs(row(ii)-row_L);
                if (min(RR)~=0)&&(min(GG)~=0)
                    mm=mm+1;
                    col_L(mm)=col(ii);
                    row_L(mm)=row(ii);
                    II(mm)=ii;
                    break;
                end
            end
            pp=ii+1;
        end
        W_rf_c=zeros(Nr,Nrf);
        F_rf_c=zeros(Nt,Nrf);
        for aaa=1:Nrf
            W_rf_c(:,aaa)=W_rf(:,find(W_bb(:,row_L(aaa))));
            F_rf_c(:,aaa)=F_rf(:,find(F_bb(:,col_L(aaa))));
        end
        F_bb_c=eye(Nrf);
        
        W_rf_cc=W_rf_c;
        F_rf_cc=F_rf_c;
        
        
        ML=zeros(Nrf^2,1);
        mmm=1;
        for rrr=1:Nrf
            for ccc=1:Nrf
                ML(mmm)=(col_L(ccc)-1)*Nt_Beam+row_L(rrr);
                mmm=mmm+1;
            end
            
        end
        
        %%
        H_wave=W'*H*F;
        h_wave=reshape(H_wave,Nr_Beam*Nt_Beam,1);
        xvec=[real(h_wave);imag(h_wave)];
        
        Y_wave=H_wave*X+Noise;
        y_wave=reshape(Y_wave,Nr_Beam*N_bar,1);
        y=[real(y_wave);imag(y_wave)];
        
        A_wave=kron((X.'),eye(Nr_Beam));
        A=[real(A_wave),(-imag(A_wave));imag(A_wave),real(A_wave)];
        A_conj=conj(A);
        
        %% GAMP-Laplace[13]
        t=1;
        x_nold=zeros(M_bar,1);
        mu_xold=2*alphaScale^2*ones(M_bar,1);
        x_n=zeros(M_bar,1);
        mu_x=zeros(M_bar,1);
        s_mold=zeros(M,1);
        
        mu_p=zeros(M,1);
        p_m=zeros(M,1);
        mu_z=zeros(M,1);
        z_m=zeros(M,1);
        mu_s=zeros(M,1);
        s_m=zeros(M,1);
        mu_r=zeros(M_bar,1);
        r_n=zeros(M_bar,1);
        
        alphaNega=zeros(M_bar,1);
        alphaPosi=zeros(M_bar,1);
        gammaNega=zeros(M_bar,1);
        gammaPosi=zeros(M_bar,1);
        psi=zeros(M_bar,1);
        
        while(t<=T_v(ttt))
            mu_p=zeros(M,1);
            p_m=zeros(M,1);
            for m=1:M
                for u=1:M_bar
                    aa=(abs(A(m,u)))^2*mu_xold(u);
                    mu_p(m)=mu_p(m)+aa;
                end
                for uu=1:M_bar
                    bb=A(m,uu)*x_nold(uu);
                    p_m(m)=p_m(m)+bb;
                end
                p_m(m)=p_m(m)-mu_p(m)*s_mold(m);
                mu_z(m)=mu_p(m)*sigman2/(mu_p(m)+sigman2);
                z_m(m)=(mu_p(m)*y(m)+sigman2*p_m(m))/(mu_p(m)+sigman2);
                mu_s(m)=(1-mu_z(m)/mu_p(m))/mu_p(m);
                s_m(m)=(z_m(m)-p_m(m))/mu_p(m);
            end
            
            for n=1:M_bar
                cc_sum=0;
                for r=1:M
                    cc=(abs(A(r,n)))^2*mu_s(r);
                    cc_sum=cc_sum+cc;
                end
                mu_r(n)=cc_sum^(-1);
                dd_sum=0;
                for rr=1:M
                    dd=A_conj(rr,n)*s_m(rr);
                    dd_sum=dd_sum+dd;
                end
                r_n(n)=x_nold(n)+mu_r(n)*dd_sum;
                
                alphaNega(n)=-r_n(n)/alphaScale-mu_r(n)/(2*alphaScale^2);
                alphaPosi(n)=r_n(n)/alphaScale-mu_r(n)/(2*alphaScale^2);
                gammaNega(n)=r_n(n)+mu_r(n)/alphaScale;
                gammaPosi(n)=r_n(n)-mu_r(n)/alphaScale;
                
                psi(n)=(exp(-alphaNega(n))*qfunc(gammaNega(n)/sqrt(mu_r(n)))...
                    +exp(-alphaPosi(n))*qfunc(-gammaPosi(n)/sqrt(mu_r(n))))...
                    /(2*alphaScale);
                
                x_n(n)=(exp(-alphaNega(n))*gammaNega(n)*qfunc(gammaNega(n)/sqrt(mu_r(n)))...
                    +exp(-alphaPosi(n))*gammaPosi(n)*qfunc(-gammaPosi(n)/sqrt(mu_r(n))))...
                    /(2*alphaScale*psi(n));
                mu_x(n)=(((gammaPosi(n))^2+mu_r(n))*exp(-alphaPosi(n))*qfunc(-gammaPosi(n)/sqrt(mu_r(n)))...
                    +((gammaNega(n))^2+mu_r(n))*exp(-alphaNega(n))*qfunc(gammaNega(n)/sqrt(mu_r(n)))...
                    -2*(mu_r(n))^2*exp(-(r_n(n))^2/(2*mu_r(n)))/(alphaScale*sqrt(2*pi*mu_r(n))))/(2*alphaScale*psi(n))...
                    -(x_n(n))^2;
            end
            
            if norm(x_n-x_nold,2)<=epsilon*norm(x_nold,2)
                break;
            else
                t=t+1;
                s_mold=s_m;
                x_nold=x_n;
                mu_xold=mu_x;
            end
        end
        
        %% Proposed
        N=zeros(M_bar,1);
        NL=sort(ML);
        LL=zeros(2*length(NL),1);
        for xx=1:length(NL)
            LL(xx)=NL(xx);
            LL(xx+length(NL))=NL(xx)+M_bar/2;
        end
        N(LL)=1;
        
        t_N=1;
        x_nold_N=zeros(M_bar,1);
        mu_xold_N=2*alphaScale_N^2*ones(M_bar,1);
        x_n_N=zeros(M_bar,1);
        mu_x_N=zeros(M_bar,1);
        s_mold_N=zeros(M,1);
        
        mu_p_N=zeros(M,1);
        p_m_N=zeros(M,1);
        mu_z_N=zeros(M,1);
        z_m_N=zeros(M,1);
        mu_s_N=zeros(M,1);
        s_m_N=zeros(M,1);
        mu_r_N=zeros(M_bar,1);
        r_n_N=zeros(M_bar,1);
        
        alphaNega_N=zeros(M_bar,1);
        alphaPosi_N=zeros(M_bar,1);
        gammaNega_N=zeros(M_bar,1);
        gammaPosi_N=zeros(M_bar,1);
        psi_N=zeros(M_bar,1);
        
        while(t_N<T_v(ttt))
            mu_p_N=zeros(M,1);
            p_m_N=zeros(M,1);
            for m=1:M
                for u=1:M_bar
                    if N(u)~=0
                        aa=(abs(A(m,u)))^2*mu_xold_N(u);
                        mu_p_N(m)=mu_p_N(m)+aa;
                    end
                end
                for uu=1:M_bar
                    if N(uu)~=0
                        bb=A(m,uu)*x_nold_N(uu);
                        p_m_N(m)=p_m_N(m)+bb;
                    end
                end
                p_m_N(m)=p_m_N(m)-mu_p_N(m)*s_mold_N(m);
                mu_z_N(m)=mu_p_N(m)*sigman2/(mu_p_N(m)+sigman2);
                z_m_N(m)=(mu_p_N(m)*y(m)+sigman2*p_m_N(m))/(mu_p_N(m)+sigman2);
                if mu_p_N(m)~=0
                    mu_s_N(m)=(1-mu_z_N(m)/mu_p_N(m))/mu_p_N(m);
                    s_m_N(m)=(z_m_N(m)-p_m_N(m))/mu_p_N(m);
                end
            end
            
            for n=1:M_bar
                cc_sum=0;
                for r=1:M
                    cc=(abs(A(r,n)))^2*mu_s_N(r);
                    cc_sum=cc_sum+cc;
                end
                mu_r_N(n)=cc_sum^(-1);
                dd_sum=0;
                for rr=1:M
                    dd=A_conj(rr,n)*s_m_N(rr);
                    dd_sum=dd_sum+dd;
                end
                r_n_N(n)=x_nold_N(n)+mu_r_N(n)*dd_sum;
            end
            
            for vv=1:M_bar
                if N(vv)~=0
                    alphaNega_N(vv)=-r_n_N(vv)/alphaScale_N-mu_r_N(vv)/(2*alphaScale_N^2);
                    alphaPosi_N(vv)=r_n_N(vv)/alphaScale_N-mu_r_N(vv)/(2*alphaScale_N^2);
                    gammaNega_N(vv)=r_n_N(vv)+mu_r_N(vv)/alphaScale_N;
                    gammaPosi_N(vv)=r_n_N(vv)-mu_r_N(vv)/alphaScale_N;
                    
                    psi_N(vv)=(exp(-alphaNega_N(vv))*qfunc(gammaNega_N(vv)/sqrt(mu_r_N(vv)))...
                        +exp(-alphaPosi_N(vv))*qfunc(-gammaPosi_N(vv)/sqrt(mu_r_N(vv))))...
                        /(2*alphaScale_N);
                    
                    x_n_N(vv)=(exp(-alphaNega_N(vv))*gammaNega_N(vv)*qfunc(gammaNega_N(vv)/sqrt(mu_r_N(vv)))...
                        +exp(-alphaPosi_N(vv))*gammaPosi_N(vv)*qfunc(-gammaPosi_N(vv)/sqrt(mu_r_N(vv))))...
                        /(2*alphaScale_N*psi_N(vv));
                    mu_x_N(vv)=(((gammaPosi_N(vv))^2+mu_r_N(vv))*exp(-alphaPosi_N(vv))*qfunc(-gammaPosi_N(vv)/sqrt(mu_r_N(vv)))...
                        +((gammaNega_N(vv))^2+mu_r_N(vv))*exp(-alphaNega_N(vv))*qfunc(gammaNega_N(vv)/sqrt(mu_r_N(vv)))...
                        -2*(mu_r_N(vv))^2*exp(-(r_n_N(vv))^2/(2*mu_r_N(vv)))/(alphaScale_N*sqrt(2*pi*mu_r_N(vv))))/(2*alphaScale_N*psi_N(vv))...
                        -(x_n_N(vv))^2;
                end
            end
            
            if norm(x_n_N-x_nold_N,2)<=epsilon*norm(x_nold_N,2)
                break;
            else
                t_N=t_N+1;
                s_mold_N=s_m_N;
                x_nold_N=x_n_N;
                mu_xold_N=mu_x_N;
            end
        end
        
        %%
        for bbb=1:Nr_Beam*Nt_Beam
            h_wave_c(bbb)=x_n(bbb)+1i*x_n(bbb+Nr_Beam*Nt_Beam);
            h_wave_c_N(bbb)=x_n_N(bbb)+1i*x_n_N(bbb+Nr_Beam*Nt_Beam);
        end
        
        H_wave_c=reshape(h_wave_c,Nr_Beam,Nt_Beam);
        H_wave_c_N=reshape(h_wave_c_N,Nr_Beam,Nt_Beam);
        H_GAMP=W*H_wave_c*F';
        H_GAMP_N=W*H_wave_c_N*F';
        
        H_c=W_rf_cc'*H*F_rf_cc;
        W_bb_c=inv(H_c);
        W_c=W_rf_cc*W_bb_c';
        Nosie_c=W_c'*n_vector(1:Nt,1:Nrf);
        Y_c=sqrt(P)*W_c'*H*F_rf_cc*F_bb_c+Nosie_c;
        vec_H_c=reshape(H_c,Nrf^2,1);
        
        H_cGAMP=W_rf_cc'*H_GAMP*F_rf_cc;
        W_bb_cGAMP=inv(H_cGAMP);
        W_cGAMP=W_rf_cc*W_bb_cGAMP';
        Nosie_GAMP=W_cGAMP'*n_vector(1:Nt,1:Nrf);
        Y_cGAMP=sqrt(P)*W_bb_cGAMP*H_c+Nosie_GAMP;
        
        H_cGAMP_N=W_rf_cc'*H_GAMP_N*F_rf_cc;
        W_bb_cGAMP_N=inv(H_cGAMP_N);
        W_cGAMP_N=W_rf_cc*W_bb_cGAMP_N';
        Nosie_GAMP_N=W_cGAMP_N'*n_vector(1:Nt,1:Nrf);
        Y_cGAMP_N=sqrt(P)*W_bb_cGAMP_N*H_c+Nosie_GAMP_N;
        
        J=H_c-H_cGAMP_N;%
        sum1=sum1+trace(J'*J)/(norm(vec_H_c))^2;% Proposed
        
        U=H_c-H_cGAMP;
        sum2=sum2+trace(U'*U)/(norm(vec_H_c))^2;% GAMP-Laplace[13]
        
    end
    NMSE_1(ttt)=10*log10(sum1/Nx);
    NMSE_2(ttt)=10*log10(sum2/Nx);
end
figure(7)
plot(T_v,NMSE_1,'-bo');
hold on
plot(T_v,NMSE_2,'-r.');
xlabel('The number of iterations');
ylabel('NMSE [dB]');
legend({'Proposed','GAMP-Laplace[13]'},'Location','northeast');


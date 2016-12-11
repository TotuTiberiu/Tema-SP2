%Am luat D ca fiind 10 deoarece am calculat formula coeficientiilor Fourier
%manual si am avut nevoie de o durata a semnalului care sa fie divizibila
%cu perioad. De asemena daca as fi luat-o 24, durata totala a semnalului ar
%fi fost peste perioada semnalului. Am incercat o implementare directa a
%coeficientilor dar nu a mers asa ca am realizat calculul pana la final.
%Conform teoriei semnalelor orice semnal se poate descompune ca o suma de
%sinusuri si cosinusuri( care se pot transforma in sinusuri). Acest lucru
%se poate realiza prin intermediul seriilor Fouriei. Astfel dupa
%reconstructie semnalul va arata aproape identic cu cel initial doar ca
%varfurile ascutite ale acestuia for disparea deoarece functiile sinus si
%cosinus ce ar putea avea frecventele atat de ridicate incat sa genereze o
%forma ascutita sunt la frecvente foarte inalte iar amplitudinile acesotra
%ar fi foarte mici. Totodata, functia sinus si cosinus nu poate avea
%variatii instante astfel, la limitele unui palier ale unui semnal
%dreptunghiular vor aparea niste spike-uri datorate suprapunerii mai multor
%sinusoide ce nu se pot redresa instantaneu.


%Construim semnalul



P=40;

D=10;

t=-120:0.1:120;

s=zeros(1,length(t));

m=1/10;

j=0;

for i=40:0.1:160
    j=j+1;
    c=floor(i/P);
    i=i-c*P;
    
    if(j==1)
        s(1,j+1200)=m*(t(j+1200));
    else if (i<=10)
            s(1,j+1200)=m*(t(j+1200))-m*40*(c-1);
            %Am adunca 40*(c-1) pentru a modela si celelalte perioade
        else if (i>=10&&i<=20)
                s(1,j+1200)=-m*(t(j+1200))+2+m*40*(c-1);
            end
        end
    end
end


j=0;

for i=40:0.1:160
    j=j+1;
    c=floor(i/P);
    i=i-c*P;
    
    if(j==1)
        s(1,j)=m*(t(j+1200));
    else if (i<=10)
            s(1,j)=m*(t(j+1200))-m*40*(c-1);
            %Am adunca 40*(c-1) pentru a modela si celelalte perioade
        else if (i>=10&&i<=20)
                s(1,j)=-m*(t(j+1200))+2+m*40*(c-1);
            end
        end
    end
end


N=50;

p=-9:1:9;

Ak=zeros(length(2*p));
Ak(10)=0.25;

for k=1:2:9
    
        %Ak(k+10)=4*0.5/(k*k*pi*pi);
        %Ak(10-k)=Ak(k+10);
        %Formula calculata dupa calculul lui Ck
        Ak(k+10)=2*sqrt(2*(-1)^k+4*(-1)^(k-1)+1)/(k*k*pi*pi);
        Ak(10-k)=Ak(k+10);
end

w0=2*pi/P;

r=linspace(-120,120,241);
f=0*r;


for k=-N+1:1:N-1
    
   if (k==0)
       continue;
   end
   
   C_k=(2*(-1j)^k-(-1)^k-1)/(k*k*pi*pi);
   f_k=C_k*exp(1j*k*w0*r);
   f=f+f_k;
   
    
end

%Adunam componenta continua
f=f+0.25;

figure (1)

plot(t,s),grid,xlabel('Timp'),ylabel('s[t]'),title('Graficul lui s in functie de timp');

figure (2)

stem(p,Ak),grid,xlabel('k'),ylabel('Ak'),title('Spectrul de amplitudini ale lui s');

figure (3)

plot(r,f),grid,xlabel('f[t]'),ylabel('Timp'),title('Graficul dependentei functiei s de timp refacut cu ajutorul seriei Fourier exponentiala');


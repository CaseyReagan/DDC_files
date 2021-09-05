clear all
close all
clc

% L=256;
% ai=sign(randn(1,L));
% aq=sign(randn(1,L));
% U=8;
% a=ai+j*aq;
% b=upsample(a,U);
% fb=rcosfir(0.35,3,U,1,'sqrt');
% 
% c=filter(fb,1,b);
% 
% d=resample(c,1866,80);
% 
% e=d.*exp(j*2*pi*46.666/186.6*(1:length(d)));
% 
% f=real(e);

load datasource.mat

%load sig_source.dat

%qf = round(f*10000);


%%
fb1=fir1(128,8/186.6);
fb2=fir1(128,2/186.6);
fb3=fir1(256,2/186.6);

l = 129;     %length
pl = 16384;   %pipe_length
res = 0;

f = f * 10000;

t = length(f);

        ip11 = f(1:pl+l);
        ip12 = f(pl+1:2*pl+l);
%        ipx11 = [zeros(1,pl),ip12(1:pl)];
%        ipx12 = [ip11(1:pl),ip12(1:pl)];       %same with f

        ip21 = ip11 .* cos(2*pi*(46.666/186.6)*(1:pl+l));
        ip22 = ip12 .* cos(2*pi*(46.666/186.6)*(pl+1:2*pl+l));

        ip1 = [ip21(1:pl),ip22(1:pl)];
%        ipx21 = [zeros(1,pl),ip22(1:pl)];
%         ipx22 = [ip21(1:pl),zeros(1,pl)];


%             for ax=1:pl
%                 ix21(ax) = ip11(ax) * cos(2*pi*(46.69/186.6)*(ax));
%             end
%             
%             for ax=1:pl
%                 ix22(ax) = ip12(ax) * cos(2*pi*(46.69/186.6)*(ax + pl));
%             end
%             
%             ipx23 = [zeros(1,pl),ix22];
        
        
        qp11 = f(1:pl+l);
        qp12 = f(pl+1:2*pl+l);

        qp21 = qp11 .* sin(2*pi*(46.666/186.6)*(1:pl+l));
        qp22 = qp12 .* sin(2*pi*(46.666/186.6)*(pl+1:2*pl+l));

        qp1 = [qp21(1:pl),qp22(1:pl)];
        

i1 = f .* cos(2*pi*(46.666/186.6)*(1:t));
q1 = f .* sin(2*pi*(46.666/186.6)*(1:t));

i2 = filter(fb1, 1, i1);
q2 = filter(fb1, 1, q1);

        ip31 = filter(fb1,1,ip21);
        ip33 = filter(fb1,1,ip22);
        ip32 = [ip31(pl+1:end),ip33(l+1:end)];

        qp31 = filter(fb1,1,qp21);
        qp33 = filter(fb1,1,qp22);
        qp32 = [qp31(pl+1:end),qp33(l+1:end)];
        
        ipx31 = [ip31(1:pl+l),ip32(l+1:pl+l)];     % same with i2

%%
y1 = 0;
y2 = 0;
y3 = 0;

x1 = 0;
x2 = 0;
x3 = 0;

comb_stage = 3;
%%

for a=1:length(i2)
    y1 = i2(a) + y1;
    y2 = y1 + y2;
    y3 = y2 + y3;
    y4(a) = y3;
    
    x1 = q2(a) + x1;
    x2 = x1 + x2;
    x3 = x2 + x3;
    x4(a) = x3;
end

        y1 = 0;
        y2 = 0;
        y3 = 0;

        x1 = 0;
        x2 = 0;
        x3 = 0;

         for a=1:length(ip31)
%             if a == pl
%                 y111 = y1;
%             end
            
            y1 = ip31(a) + y1;
            y2 = y1 + y2;
            y3 = y2 + y3;
            ip41(a) = y3;

            x1 = qp31(a) + x1;
            x2 = x1 + x2;
            x3 = x2 + x3;
            qp41(a) = x3;
            
            if a == pl
                y11 = y1;
                y22 = y2;
                y33 = y3;

                x11 = x1;
                x22 = x2;
                x33 = x3;
            end
        end
        
        y1 = y11;
        y2 = y22;
        y3 = y33;

        x1 = x11;
        x2 = x22;
        x3 = x33;

        for a=1:length(ip32)
            y1 = ip32(a) + y1;
            y2 = y1 + y2;
            y3 = y2 + y3;
            ip42(a) = y3;

            x1 = qp32(a) + x1;
            x2 = x1 + x2;
            x3 = x2 + x3;
            qp42(a) = x3;
        end

%         ip2 = [ip41(1:pl),ip42(1:pl)];
%         qp2 = [qp41(1:pl),qp42(1:pl)];
        
%         ipx41 = [zeros(1,pl),ip42(1:pl)];
%         ipx42 = [zeros(1,pl),qp42(1:pl)];
%%
% t2 = ceil(length(y4)/16);
% 
% y5 = zeros(1,t2);
% x5 = zeros(1,t2);
% 
% for b=1:t2
%     y5(b) = y4(1+16*(b-1));
%     x5(b) = x4(1+16*(b-1));
% end

y5 = y4(1:16:end);
x5 = x4(1:16:end);

% td1 =16*floor(length(ip41)/16);
% td2 = 16 - rem(length(pl),16);
% 
% ip51 = ip41(1:16:td1);
% ip52 = ip42(td2:16:)

        ip51 = ip41(1:16:16512);
        ip511 = ip41(1:16:16384);
        ip52 = ip42(16:16:16511);
        ip522 = ip42(1:16:16384);

        qp51 = qp41(1:16:16512);
        qp511 = ip41(1:16:16384);
        qp52 = qp42(16:16:16511);
        qp522 = ip42(1:16:16384);
        
%        ip3 = [ip511,ip522];

%%

x6 = 0;
x7 = 0;
x8 = 0;

y6 = 0;
y7 = 0;
y8 = 0;

ym1 = 0;
ym2 = 0;
ym3 = 0;

xm1 = 0;
xm2 = 0;
xm3 = 0;

for c=1:length(y5)
%    for cs=1:comb_stage
        ym1 = y5(c) - y6;
        y6 = y5(c);
        ym2 = ym1 - y7;
        y7 = ym1;
        ym3 = ym2 - y8;
        y8 = ym2;
        y(c) = ym3;
        
        xm1 = x5(c) - x6;
        x6 = x5(c);
        xm2 = xm1 - x7;
        x7 = xm1;
        xm3 = xm2 - x8;
        x8 = xm2;
        x(c) = xm3;
end

% for c=1:length(y5)
%     y1 = y5(c) - y1;
%     y(c) = y1;
%     
%     x1 = x5(c) - x1;
%     x(c) = x1;
% end

%%
td1 = 16/186.6;
td2 = 1/8;

t1 = 1:td1:length(y)*td1;
t2 = 1:td2:length(y)*td1;

tl = td1;
d_t = t2(2) - t1(2);    % td2 - td1

for d=1:length(t2)
    
    if d == 1
        u = 0;
        j = 1;
    else if d == 2
            u = d_t;
        end
    end
    
    if d > 2
        u = u + d_t;
        if u > tl
            u = u - tl;
            j = j + 1;
        end
    end
    yy(d) = u / tl * (y(j+1) - y(j)) + y(j);
    j = j + 1;
end

for e=1:length(t2)
    
    if e == 1
        u = 0;
        j = 1;
    else if e == 2
            u = d_t;
        end
    end
    
    if e > 2
        u = u + d_t;
        if u > tl
            u = u - tl;
            j = j + 1;
        end
    end
    xx(e) = u / tl * (x(j+1) - x(j)) + x(j);
    j = j + 1;
end

fir = rcosfir(0.35,3,8,1,'sqrt');

% yy2 = filter(fir, 1, yy);
% xx2 = filter(fir, 1, xx);

zz=filter(fir,1,xx+sqrt(-1)*yy);
% zz=yy2+sqrt(-1)*xx2;

